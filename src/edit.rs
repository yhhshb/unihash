//! Deletion and defragmentation.
//!
//! ## Deletion
//!
//! Deleting the k-mers at window starts `lo..=hi` of a piece `P` (length `L`)
//! cuts `P` into
//!
//! * a left part `P[0, lo + k - 1)` (k-mers `0..lo`), which keeps `P`'s id —
//!   only its length shrinks, so its index entries stay valid — and
//! * a right part `P[hi + 1, L)` (k-mers `hi + 1..`), a new piece *viewing the
//!   same bases* (no copy), whose index items are re-filed under its own id.
//!
//! Each part keeps the `k - 1` bases that flank the deleted run, so every
//! k-mer outside the run survives. Bases covered by neither part become dead
//! space in the base array until [`KmerTable::defrag`].
//!
//! Index maintenance is exact: the index items of a window (layer-1 entries,
//! layer-2 entries, or its layer-3 location) are a function of the window's
//! bases and of which minimizers are hot, so the items of the windows around
//! the cut are recomputed and only those no surviving window needs are
//! removed.
//!
//! A k-mer repeated inside one inserted run is stored twice (plan decision
//! 5), so deleting one occurrence is not enough: [`KmerTable::delete_read`]
//! re-checks its input with a streaming query and deletes leftover copies.
//!
//! ## Defragmentation
//!
//! [`KmerTable::defrag`] merges pieces whose ends overlap by exactly `k - 1`
//! bases (in either orientation), rewrites the base array without dead
//! space, renumbers pieces, and rebuilds the indexes from the pieces. Merging
//! two such pieces adds no k-mer (the junction windows are the pieces' own
//! end windows), so any pairing of matching ends is correct. It is never
//! called automatically: see [`KmerTable::fragmentation`].

use crate::codec::{canonical, canonicalize, encode, revcomp_kmer, RollingKmer, INVALID};
use crate::minimizer::{kmer_minimizer_mask, order_value, WindowMin};
use crate::piece_store::{PieceId, PieceStore};
use crate::table::{Hit, KmerTable, QueryTrace};
use std::fmt;

/// Cumulative deletion and defragmentation counters.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct EditStats {
    /// `delete_read` / `delete` calls.
    pub delete_calls: u64,
    /// Stored k-mer occurrences removed.
    pub kmers_deleted: u64,
    /// Contiguous runs cut out of pieces.
    pub cuts: u64,
    /// Cuts that left both a left and a right part.
    pub splits: u64,
    /// Index items re-filed under a right part's new id.
    pub items_moved: u64,
    /// Extra passes of `delete_read` caused by duplicated k-mers.
    pub recheck_passes: u64,
    /// Bases no piece covers any more (reset by `defrag`).
    pub dead_bases: u64,
    /// Bases shared by the two parts of a cut (counted twice in live bases;
    /// reset by `defrag`).
    pub shared_bases: u64,
    /// Cuts since the last `defrag`.
    pub cuts_since_defrag: u64,
    /// Pieces created by insertion since the last `defrag`.
    pub inserted_pieces_since_defrag: u64,
    /// Number of `defrag` calls.
    pub defrags: u64,
    /// Number of `defrag_incremental` calls.
    pub incremental_defrags: u64,
    /// Piece joins made by all `defrag` calls.
    pub joins: u64,
}

impl fmt::Display for EditStats {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "  deletions                    {} calls, {} k-mers, {} cuts ({} splits), {} items moved, {} re-check passes",
            self.delete_calls,
            self.kmers_deleted,
            self.cuts,
            self.splits,
            self.items_moved,
            self.recheck_passes
        )?;
        write!(
            f,
            "  fragmentation                {} dead bases, {} shared bases, {} cuts and {} new pieces since defrag ({} full + {} incremental defrags, {} joins)",
            self.dead_bases,
            self.shared_bases,
            self.cuts_since_defrag,
            self.inserted_pieces_since_defrag,
            self.defrags,
            self.incremental_defrags,
            self.joins
        )
    }
}

/// Snapshot of how far the table is from its defragmented layout, to decide
/// when to call [`KmerTable::defrag`] (see [`Fragmentation::overhead`]).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Fragmentation {
    /// Non-empty pieces.
    pub pieces: u64,
    /// Pieces emptied by deletion (still holding an id).
    pub empty_pieces: u64,
    /// Non-empty pieces right after the last `defrag` (or 0).
    pub pieces_at_last_defrag: u64,
    /// Bases in the array, live bases (Σ lengths), dead bases.
    pub array_bases: u64,
    pub live_bases: u64,
    pub dead_bases: u64,
    pub shared_bases: u64,
    pub cuts_since_defrag: u64,
    pub inserted_pieces_since_defrag: u64,
    /// Index block bytes beyond the bytes per entry of the last rebuild
    /// (freed slots, partly empty overflow blocks), both index layers.
    pub index_slack_bytes: u64,
    /// The part of the index slack `defrag_incremental` returns: the slack
    /// of each index layer once it reaches [`INDEX_SHRINK_SLACK`] of that
    /// layer (it then re-sizes the layer).
    pub shrinkable_slack_bytes: u64,
    /// Heap bytes of pieces + layer-1/2 indexes + layer 3 (by length).
    pub table_bytes: u64,
}

impl Fragmentation {
    /// Bytes held by dead bases (2 bits each): what `defrag_incremental`
    /// reclaims by sliding the base array.
    pub fn dead_bytes(&self) -> f64 {
        self.dead_bases as f64 / 4.0
    }

    /// Bytes held by pieces beyond the count left by the last full `defrag`
    /// — for each, a piece record, its k − 1 overlap bases and about one
    /// junction index entry (≈ 13 B) — plus bases shared by the two parts of
    /// a cut and emptied piece records. Merging pieces returns them: a full
    /// `defrag` merges all, `defrag_incremental` only chains of short pieces.
    /// On a table that also grows by new sequence, new pieces that extend
    /// nothing are counted too, so this is an upper bound.
    pub fn piece_bytes(&self, k: usize) -> f64 {
        let extra = self.pieces.saturating_sub(self.pieces_at_last_defrag) as f64;
        extra * (8.0 + (k - 1) as f64 / 4.0 + 13.0)
            + self.shared_bases as f64 / 4.0
            + self.empty_pieces as f64 * 8.0
    }

    /// Estimated bytes a full `defrag` would reclaim: dead bytes, piece
    /// bytes and the index slack (block bytes beyond the bytes per entry
    /// right after the last index rebuild: slots freed by deletions, partly
    /// empty overflow blocks).
    pub fn reclaimable_bytes(&self, k: usize) -> f64 {
        self.dead_bytes() + self.piece_bytes(k) + self.index_slack_bytes as f64
    }

    /// [`reclaimable_bytes`](Self::reclaimable_bytes) relative to the table
    /// size: compare with a threshold before calling `defrag`.
    pub fn overhead(&self, k: usize) -> f64 {
        self.reclaimable_bytes(k) / self.table_bytes.max(1) as f64
    }

    /// What `defrag_incremental` reclaims for sure, relative to the table
    /// size: dead bytes and the index slack it re-sizes away
    /// ([`shrinkable_slack_bytes`](Self::shrinkable_slack_bytes)).
    pub fn incremental_overhead(&self) -> f64 {
        (self.dead_bytes() + self.shrinkable_slack_bytes as f64) / self.table_bytes.max(1) as f64
    }

    /// Piece bytes relative to the table size: the part left to a full
    /// `defrag` (except chains of short pieces).
    pub fn piece_overhead(&self, k: usize) -> f64 {
        self.piece_bytes(k) / self.table_bytes.max(1) as f64
    }

    /// Index slack relative to the table size.
    pub fn index_slack(&self) -> f64 {
        self.index_slack_bytes as f64 / self.table_bytes.max(1) as f64
    }
}

/// Result of one [`KmerTable::defrag`] call.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct DefragReport {
    pub pieces_before: u64,
    pub pieces_after: u64,
    pub joins: u64,
    pub array_bases_before: u64,
    pub array_bases_after: u64,
    pub entries_before: u64,
    pub entries_after: u64,
    pub table_bytes_before: u64,
    pub table_bytes_after: u64,
}

impl fmt::Display for DefragReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "defrag: pieces {} → {} ({} joins), bases {} → {}, entries {} → {}, table {} → {} B",
            self.pieces_before,
            self.pieces_after,
            self.joins,
            self.array_bases_before,
            self.array_bases_after,
            self.entries_before,
            self.entries_after,
            self.table_bytes_before,
            self.table_bytes_after
        )
    }
}

/// Index item of a window: `(layer, key, position)`. Layers 1 and 2: the
/// minimizer order value and the ℓ- (ℓ2-) mer offset in the piece. Layer 3:
/// the canonical k-mer and the window start.
type Item = (u8, u64, u32);

/// A contiguous run of found k-mers in one piece (deletion pass).
#[derive(Clone, Copy, Debug)]
struct Run {
    pid: PieceId,
    first: usize,
    last: usize,
    /// The piece holds the query's k-mers in the query's orientation.
    same: bool,
}

const NONE: u64 = u64::MAX;

/// `defrag_incremental` re-sizes an index when its slack is at least this
/// fraction of the index bytes.
pub const INDEX_SHRINK_SLACK: f64 = 1.0 / 16.0;

fn index_slack_worth_shrinking(ix: &crate::index::MinimizerIndex) -> bool {
    ix.slack_bytes() as f64 >= INDEX_SHRINK_SLACK * ix.heap_bytes().0 as f64 && ix.slack_bytes() > 0
}

/// Merge limit of `defrag_incremental` (0 = no merging). A join saves ≈ 28 B
/// (piece record, k − 1 overlap bases, a junction entry) but costs 8 index
/// lookups and re-filing every window of the chain; a full `defrag`, which
/// processes every base anyway, joins more cheaply per reclaimed byte.
pub const INCREMENTAL_MERGE_MEAN_LEN: usize = 0;

impl KmerTable {
    fn assert_no_open_read(&self, what: &str) {
        assert!(
            !self.stream.as_ref().is_some_and(|s| s.in_read),
            "{what}: a read is open (call end_read first)"
        );
    }

    /// Where a k-mer (either orientation) is stored, if present.
    pub fn locate(&self, kmer: u64) -> Option<Hit> {
        if kmer & !crate::codec::kmer_mask(self.cfg.k) != 0 {
            return None;
        }
        self.locate_canonical(canonicalize(kmer, self.cfg.k))
    }

    fn locate_canonical(&self, x: u64) -> Option<Hit> {
        let (m, ties) = kmer_minimizer_mask(x, self.cfg.k, self.cfg.l);
        self.find_with_minimizer(x, m, ties, &mut QueryTrace::default())
    }

    // ------------------------------------------------------------------
    // Deletion
    // ------------------------------------------------------------------

    /// Delete one k-mer (either orientation), every stored occurrence of it.
    /// Returns whether it was present.
    ///
    /// # Panics
    /// If a read is open (`insert_chunk` without `end_read`).
    pub fn delete(&mut self, kmer: u64) -> bool {
        self.assert_no_open_read("delete");
        self.edits.delete_calls += 1;
        let Some(mut h) = self.locate(kmer) else {
            return false;
        };
        let x = canonicalize(h.kmer, self.cfg.k);
        loop {
            self.delete_windows(h.piece, h.pos as usize, h.pos as usize);
            match self.locate_canonical(x) {
                Some(n) => h = n,
                None => return true,
            }
        }
    }

    /// Delete every k-mer of `read` (split at non-ACGT characters like
    /// insertion). Absent k-mers are ignored. Successive k-mers found at
    /// successive positions of a piece are removed as one run, so a stretch
    /// costs one cut per piece it touches. Returns the number of stored
    /// k-mer occurrences removed.
    ///
    /// # Panics
    /// If a read is open (`insert_chunk` without `end_read`).
    pub fn delete_read(&mut self, read: &[u8]) -> u64 {
        self.assert_no_open_read("delete_read");
        self.edits.delete_calls += 1;
        let before = self.edits.kmers_deleted;
        loop {
            self.delete_pass(read);
            // Copies of a deleted k-mer elsewhere in its piece (intra-run
            // repeats) survive the pass; look for any k-mer still present.
            let mut any = false;
            let mut sq = self.query_stream();
            sq.push(read, |_, p| any |= p);
            if !any {
                break;
            }
            self.edits.recheck_passes += 1;
        }
        self.edits.kmers_deleted - before
    }

    fn delete_pass(&mut self, read: &[u8]) {
        let k = self.cfg.k;
        let mut roll = RollingKmer::new(k);
        let mut run: Option<Run> = None;
        for &c in read {
            let b = encode(c);
            if b == INVALID {
                if let Some(r) = run.take() {
                    self.cut(r);
                }
                roll.reset();
                continue;
            }
            if !roll.push(b) {
                continue;
            }
            if let Some(r) = run.as_mut() {
                // Does the next k-mer sit right after the run (one base)?
                let (g0, len) = self.store.view(r.pid);
                let ok = if r.same {
                    r.last + k < len as usize
                        && self.store.code_at_global(g0 + (r.last + k) as u64) == b
                } else {
                    r.last >= 1 && self.store.code_at_global(g0 + r.last as u64 - 1) == 3 - b
                };
                if ok {
                    if r.same {
                        r.last += 1;
                    } else {
                        r.last -= 1;
                    }
                    continue;
                }
                let r = run.take().expect("run");
                self.cut(r);
            }
            let (fwd, rc) = (roll.fwd(), roll.rc());
            if let Some(h) = self.locate_canonical(canonical(fwd, rc)) {
                run = Some(Run {
                    pid: h.piece,
                    first: h.pos as usize,
                    last: h.pos as usize,
                    same: h.kmer == fwd,
                });
            }
        }
        if let Some(r) = run {
            self.cut(r);
        }
    }

    fn cut(&mut self, r: Run) {
        let (lo, hi) = if r.first <= r.last {
            (r.first, r.last)
        } else {
            (r.last, r.first)
        };
        self.delete_windows(r.pid, lo, hi);
    }

    /// Index items of the windows `s_from..=s_to` of piece `pid` (appended to
    /// `out`, unsorted). `vals` is scratch.
    fn window_items(
        &self,
        pid: PieceId,
        s_from: usize,
        s_to: usize,
        vals: &mut Vec<u64>,
        out: &mut Vec<Item>,
    ) {
        let (k, l, w, l2) = (self.cfg.k, self.cfg.l, self.w, self.l2);
        if s_from > s_to {
            return;
        }
        let g0 = self.store.start(pid);
        vals.clear();
        let mut roll = RollingKmer::new(l);
        for g in g0 + s_from as u64..g0 + (s_to + k) as u64 {
            if roll.push(self.store.code_at_global(g)) {
                vals.push(order_value(roll.canonical()));
            }
        }
        let mut wm = WindowMin::new(w);
        let mut last: Option<(u64, bool)> = None;
        for (i, &v) in vals.iter().enumerate() {
            wm.push(v);
            let Some((m, ties)) = wm.current() else {
                continue;
            };
            let s = s_from + i + 1 - w;
            let hot = l2 > 0
                && match last {
                    Some((lm, h)) if lm == m => h,
                    _ => {
                        let h = self.index.is_hot(m);
                        last = Some((m, h));
                        h
                    }
                };
            if hot {
                let f = self.store.kmer_at(pid, s as u32, k);
                let (m2, t2) = kmer_minimizer_mask(f, k, l2);
                if self.index2.is_hot(m2) {
                    out.push((3, canonicalize(f, k), s as u32));
                } else {
                    let mut bits = t2;
                    while bits != 0 {
                        out.push((2, m2, (s + bits.trailing_zeros() as usize) as u32));
                        bits &= bits - 1;
                    }
                }
            } else {
                let mut bits = ties;
                while bits != 0 {
                    out.push((1, m, (s + bits.trailing_zeros() as usize) as u32));
                    bits &= bits - 1;
                }
            }
        }
    }

    /// Delete the k-mers at window starts `lo..=hi` of piece `pid` (see the
    /// module documentation).
    pub(crate) fn delete_windows(&mut self, pid: PieceId, lo: usize, hi: usize) {
        let (k, l) = (self.cfg.k, self.cfg.l);
        let span = k - l; // farthest an entry lies from its window start
        let (g0, len) = self.store.view(pid);
        let len = len as usize;
        assert!(len >= k && lo <= hi && hi <= len - k, "bad deletion range");
        let nw = len - k + 1;

        let mut vals = Vec::new();
        // Items of the windows that lose their piece (cut run and its left
        // neighbourhood), of the right part, and those the left part keeps.
        let mut cut = Vec::new();
        self.window_items(pid, lo.saturating_sub(span), hi, &mut vals, &mut cut);
        let mut moved = Vec::new();
        if hi + 1 < nw {
            self.window_items(pid, hi + 1, nw - 1, &mut vals, &mut moved);
        }
        let mut keep = Vec::new();
        if lo > 0 {
            self.window_items(
                pid,
                lo.saturating_sub(2 * span),
                lo - 1,
                &mut vals,
                &mut keep,
            );
        }
        for v in [&mut cut, &mut moved, &mut keep] {
            v.sort_unstable();
            v.dedup();
        }

        // 0. Room for the items both parts need (copies), reserved while every
        //    entry still matches the old piece layout (growth recomputes
        //    minimizers from the pieces).
        let (mut n1, mut n2) = (0u64, 0u64);
        for it in moved.iter().filter(|it| keep.binary_search(it).is_ok()) {
            n1 += (it.0 == 1) as u64;
            n2 += (it.0 == 2) as u64;
        }
        if n1 + n2 > 0 {
            let store = &self.store;
            let l2 = self.l2;
            self.index.reserve(n1, |p, q| {
                order_value(canonicalize(store.kmer_at(p, q, l), l))
            });
            self.index2.reserve(n2, |p, q| {
                order_value(canonicalize(store.kmer_at(p, q, l2), l2))
            });
        }

        // 1. Re-point the pieces: the left part keeps the id.
        let a_len = if lo > 0 { lo + k - 1 } else { 0 };
        let b_off = hi + 1;
        let b_len = if hi + 1 < nw { len - b_off } else { 0 };
        let right = if a_len > 0 {
            self.store.set_view(pid, g0, a_len as u32);
            (b_len > 0).then(|| self.store.add_view(g0 + b_off as u64, b_len as u32))
        } else {
            self.store.set_view(pid, g0 + b_off as u64, b_len as u32);
            (b_len > 0).then_some(pid)
        };

        // 2. Items only the deleted windows needed: remove.
        for &it in &cut {
            if keep.binary_search(&it).is_ok() || moved.binary_search(&it).is_ok() {
                continue;
            }
            let (layer, key, q) = it;
            let n = match layer {
                1 => self.index.remove(key, pid, q),
                2 => self.index2.remove(key, pid, q),
                _ => self.layer3.remove(key, (pid, q)) as usize,
            };
            debug_assert!(n > 0, "index item {it:?} missing");
        }

        // 3. The right part's items: re-file under its id in place (the
        //    bucket depends only on the key), or add a copy when the left
        //    part needs the item too (bases shared by both parts).
        if right == Some(pid) {
            // The right part keeps the id with shifted offsets: an in-place
            // re-key could make a moved entry look like one not yet moved
            // (same id and offset, tag collision), so remove all, then add.
            let off = b_off as u32;
            for &(layer, key, q) in &moved {
                match layer {
                    1 => {
                        self.index.remove(key, pid, q);
                    }
                    2 => {
                        self.index2.remove(key, pid, q);
                    }
                    _ => {
                        self.layer3.remove(key, (pid, q));
                    }
                }
            }
            for &(layer, key, q) in &moved {
                match layer {
                    1 => {
                        self.index.insert(key, pid, q - off);
                    }
                    2 => {
                        self.index2.insert(key, pid, q - off);
                    }
                    _ => self.layer3.add(key, (pid, q - off)),
                }
            }
            self.edits.items_moved += moved.len() as u64;
        } else if let Some(b) = right {
            let off = b_off as u32;
            for &it in &moved {
                let (layer, key, q) = it;
                let copy = keep.binary_search(&it).is_ok();
                match layer {
                    1 if copy => {
                        self.index.insert(key, b, q - off);
                    }
                    1 => {
                        let ok = self.index.rekey(key, pid, q, b, q - off);
                        debug_assert!(ok, "index item {it:?} missing");
                    }
                    2 if copy => {
                        self.index2.insert(key, b, q - off);
                    }
                    2 => {
                        // Layer 2 may hold duplicate entries: move them all.
                        while self.index2.rekey(key, pid, q, b, q - off) {}
                    }
                    _ => {
                        if !copy {
                            self.layer3.remove(key, (pid, q));
                        }
                        self.layer3.add(key, (b, q - off));
                    }
                }
            }
            self.edits.items_moved += moved.len() as u64;
        }

        // Accounting.
        let both = a_len > 0 && b_len > 0;
        let overlap = if both { a_len.saturating_sub(b_off) } else { 0 };
        let union = a_len + b_len - overlap;
        self.edits.dead_bases += (len - union) as u64;
        self.edits.shared_bases += overlap as u64;
        let n = (hi - lo + 1) as u64;
        self.stored_kmers -= n;
        self.edits.kmers_deleted += n;
        self.edits.cuts += 1;
        self.edits.cuts_since_defrag += 1;
        if a_len > 0 {
            self.mark_dirty(pid);
        }
        if let Some(b) = right {
            self.mark_dirty(b);
        }
        self.edits.splits += both as u64;
    }

    // ------------------------------------------------------------------
    // Defragmentation
    // ------------------------------------------------------------------

    /// How far the table is from its defragmented layout.
    pub fn fragmentation(&self) -> Fragmentation {
        let e = &self.edits;
        let st = self.store.len() as u64 - self.store.empty_pieces();
        Fragmentation {
            pieces: st,
            empty_pieces: self.store.empty_pieces(),
            pieces_at_last_defrag: self.pieces_at_defrag,
            array_bases: self.store.total_bases(),
            live_bases: self.store.live_bases(),
            dead_bases: e.dead_bases,
            shared_bases: e.shared_bases,
            cuts_since_defrag: e.cuts_since_defrag,
            inserted_pieces_since_defrag: e.inserted_pieces_since_defrag,
            index_slack_bytes: self.index.slack_bytes() + self.index2.slack_bytes(),
            shrinkable_slack_bytes: [&self.index, &self.index2]
                .iter()
                .filter(|ix| index_slack_worth_shrinking(ix))
                .map(|ix| ix.slack_bytes())
                .sum(),
            table_bytes: self.table_bytes() as u64,
        }
    }

    /// Heap bytes (by length) of the pieces, both indexes and layer 3.
    pub fn table_bytes(&self) -> usize {
        self.store.heap_bytes().0
            + self.index.heap_bytes().0
            + self.index2.heap_bytes().0
            + self.layer3.heap_bytes()
            + self.dirty.len() * 8
    }

    /// Merge pieces that overlap by `k - 1` bases, drop dead bases and
    /// emptied pieces, and rebuild the indexes. O(stored bases); never
    /// called automatically. Answers are unchanged.
    ///
    /// Peak memory: the indexes are freed first, so the peak is the larger
    /// of the base array plus ≈ 56 B per piece (pairing ends), the old and
    /// new base arrays (copying), and the new table being re-indexed.
    ///
    /// # Panics
    /// If a read is open (`insert_chunk` without `end_read`).
    pub fn defrag(&mut self) -> DefragReport {
        self.assert_no_open_read("defrag");
        let k = self.cfg.k;
        let km1 = k - 1;
        let cap = self.cfg.max_piece_len as u64;
        let n = self.store.len();
        let mut rep = DefragReport {
            pieces_before: n as u64 - self.store.empty_pieces(),
            array_bases_before: self.store.total_bases(),
            entries_before: self.index.num_entries() + self.index2.num_entries(),
            table_bytes_before: self.table_bytes() as u64,
            ..Default::default()
        };

        // 0. Free the indexes and layer 3 first: they are rebuilt from the
        //    new pieces, and nothing below needs them. The temporaries of
        //    steps 1-2 (≈ 56 B per piece) then never coexist with the index.
        let n1 = self.index.num_entries();
        let n2 = self.index2.num_entries();
        self.index.release();
        self.index2.release();
        self.layer3.clear();

        // 1. Ends of every piece: canonical (k-1)-mer, side, orientation.
        //    Two ends can be joined iff their keys are equal and
        //    (sides differ) == (both or neither are in canonical orientation).
        let mut ends: Vec<(u64, u32, u8)> = Vec::with_capacity(2 * rep.pieces_before as usize);
        for pid in 0..n as PieceId {
            let (g0, len) = self.store.view(pid);
            if (len as usize) < k {
                continue;
            }
            let pre = self.store.kmer_at_global(g0, km1);
            let suf = self.store.kmer_at_global(g0 + len as u64 - km1 as u64, km1);
            for (side, z) in [(0u8, pre), (1u8, suf)] {
                let c = canonicalize(z, km1);
                if c == revcomp_kmer(c, km1) {
                    continue; // palindromic overlap: orientation ambiguous
                }
                ends.push((c, pid, side | (((z == c) as u8) << 1)));
            }
        }
        ends.sort_unstable();

        // 2. Pair ends greedily; union-find keeps chains acyclic and within
        //    max_piece_len. partner[2 * pid + side] = 2 * pid' + side'.
        let mut partner = vec![NONE; 2 * n];
        let mut parent: Vec<u32> = (0..n as u32).collect();
        let mut comp_len: Vec<u32> = (0..n as PieceId).map(|p| self.store.length(p)).collect();
        fn find(parent: &mut [u32], mut x: u32) -> u32 {
            while parent[x as usize] != x {
                let p = parent[x as usize];
                parent[x as usize] = parent[p as usize];
                x = p;
            }
            x
        }
        let mut i = 0;
        while i < ends.len() {
            let mut j = i;
            while j < ends.len() && ends[j].0 == ends[i].0 {
                j += 1;
            }
            for a in i..j {
                let (_, pa, ea) = ends[a];
                let ia = 2 * pa as usize + (ea & 1) as usize;
                if partner[ia] != NONE {
                    continue;
                }
                // Scan a bounded number of candidates (repeat-rich keys can
                // gather many ends; any compatible pairing is correct).
                for &(_, pb, eb) in &ends[a + 1..j.min(a + 65)] {
                    let ib = 2 * pb as usize + (eb & 1) as usize;
                    if partner[ib] != NONE || pa == pb {
                        continue;
                    }
                    let sides_differ = (ea & 1) != (eb & 1);
                    let flags_equal = (ea >> 1) == (eb >> 1);
                    if sides_differ != flags_equal {
                        continue;
                    }
                    let (ra, rb) = (find(&mut parent, pa), find(&mut parent, pb));
                    if ra == rb {
                        continue;
                    }
                    let joined =
                        comp_len[ra as usize] as u64 + comp_len[rb as usize] as u64 - km1 as u64;
                    if joined > cap {
                        continue;
                    }
                    parent[rb as usize] = ra;
                    comp_len[ra as usize] = joined as u32;
                    partner[ia] = ib as u64;
                    partner[ib] = ia as u64;
                    rep.joins += 1;
                    break;
                }
            }
            i = j;
        }
        drop(ends);
        drop(parent);
        drop(comp_len);

        // 3. Write the chains into a new base array.
        let old = std::mem::take(&mut self.store);
        let mut store = PieceStore::with_capacity(old.live_bases(), rep.pieces_before as usize);
        let mut visited = vec![false; n];
        let copy = |store: &mut PieceStore, pid: PieceId, rev: bool, skip: usize| {
            let (g0, len) = old.view(pid);
            let len = len as u64;
            if rev {
                for i in (0..len - skip as u64).rev() {
                    store.push_code(3 - old.code_at_global(g0 + i));
                }
            } else {
                for i in skip as u64..len {
                    store.push_code(old.code_at_global(g0 + i));
                }
            }
        };
        for start in 0..n {
            if visited[start] || (old.length(start as PieceId) as usize) < k {
                continue;
            }
            // Start chains at a free end (union-find rules out cycles).
            let enter = if partner[2 * start] == NONE {
                0
            } else if partner[2 * start + 1] == NONE {
                1
            } else {
                continue;
            };
            let g = store.end();
            let (mut pid, mut side) = (start, enter);
            let mut first = true;
            loop {
                visited[pid] = true;
                copy(
                    &mut store,
                    pid as PieceId,
                    side == 1,
                    if first { 0 } else { km1 },
                );
                first = false;
                let next = partner[2 * pid + (1 - side)];
                if next == NONE {
                    break;
                }
                pid = (next / 2) as usize;
                side = (next % 2) as usize;
            }
            store.close_piece(g);
        }
        debug_assert!((0..n).all(|p| visited[p] || (old.length(p as PieceId) as usize) < k));
        drop(partner);
        drop(visited);
        drop(old);
        self.store = store;

        // 4. Re-index every piece (same code path as insertion), into heads
        //    sized for ≈ 80 % load (allocated only now, after the old base
        //    array is gone).
        let buckets =
            |e: u64| ((e as f64 / (0.8 * crate::index::BLOCK_SLOTS as f64)).ceil() as usize).max(1);
        self.index.clear_entries(buckets(n1));
        self.index2.clear_entries(buckets(n2));
        self.stored_kmers = 0;
        let mut st = self.stream.take().expect("insertion state present");
        for pid in 0..self.store.len() as PieceId {
            let len = self.store.length(pid) as usize;
            self.stored_kmers += (len - k + 1) as u64;
            self.lmer_values_of_piece(pid, &mut st.vals);
            self.index_piece(pid, &mut st);
        }
        st.shrink_to_fit();
        self.stream = Some(st);
        self.index.compact();
        self.index2.compact();
        self.index.mark_baseline();
        self.index2.mark_baseline();

        let e = &mut self.edits;
        e.defrags += 1;
        e.joins += rep.joins;
        e.dead_bases = 0;
        e.shared_bases = 0;
        e.cuts_since_defrag = 0;
        e.inserted_pieces_since_defrag = 0;
        self.dirty = Vec::new();
        rep.pieces_after = self.store.len() as u64;
        self.pieces_at_defrag = rep.pieces_after;
        rep.array_bases_after = self.store.total_bases();
        rep.entries_after = self.index.num_entries() + self.index2.num_entries();
        rep.table_bytes_after = self.table_bytes() as u64;
        rep
    }

    /// Order values of every ℓ-mer of piece `pid`, read from the base array.
    fn lmer_values_of_piece(&self, pid: PieceId, out: &mut Vec<u64>) {
        let l = self.cfg.l;
        let (g0, len) = self.store.view(pid);
        out.clear();
        let mut roll = RollingKmer::new(l);
        for g in g0..g0 + len as u64 {
            if roll.push(self.store.code_at_global(g)) {
                out.push(order_value(roll.canonical()));
            }
        }
    }
}

impl KmerTable {
    /// Diagnostic, O(stored bases): check that the three index layers hold
    /// exactly the items the pieces require — every window findable, no
    /// entry left behind by deletions. Returns a description of the first
    /// discrepancy.
    pub fn verify_index(&self) -> Result<(), String> {
        use std::collections::HashSet;
        let (k, l, l2) = (self.cfg.k, self.cfg.l, self.l2);
        let mut want: HashSet<(u8, u64, PieceId, u32)> = HashSet::new();
        let mut vals = Vec::new();
        let mut items = Vec::new();
        for pid in 0..self.store.len() as PieceId {
            let len = self.store.length(pid) as usize;
            if len < k {
                continue;
            }
            items.clear();
            self.window_items(pid, 0, len - k, &mut vals, &mut items);
            want.extend(items.iter().map(|&(a, b, q)| (a, b, pid, q)));
        }
        let mut have: HashSet<(u8, u64, PieceId, u32)> = HashSet::new();
        let mut n1 = 0u64;
        for e in self.index.iter_entries() {
            if e.q as usize + l > self.store.length(e.piece_id) as usize {
                return Err(format!("layer-1 entry {e:?} points past its piece"));
            }
            let m = order_value(canonicalize(self.store.kmer_at(e.piece_id, e.q, l), l));
            have.insert((1, m, e.piece_id, e.q));
            n1 += 1;
        }
        if n1 != self.index.num_entries() {
            return Err(format!(
                "layer-1 count {} != {n1} entries",
                self.index.num_entries()
            ));
        }
        if l2 > 0 {
            for e in self.index2.iter_entries() {
                if e.q as usize + l2 > self.store.length(e.piece_id) as usize {
                    return Err(format!("layer-2 entry {e:?} points past its piece"));
                }
                let m = order_value(canonicalize(self.store.kmer_at(e.piece_id, e.q, l2), l2));
                have.insert((2, m, e.piece_id, e.q));
            }
        }
        for (x, (pid, s)) in self.layer3.iter() {
            if s as usize + k > self.store.length(pid) as usize
                || canonicalize(self.store.kmer_at(pid, s, k), k) != x
            {
                return Err(format!("layer-3 location of {x:#x} is stale"));
            }
            have.insert((3, x, pid, s));
        }
        if let Some(m) = want.iter().find(|i| !have.contains(i)) {
            return Err(format!("missing index item {m:?}"));
        }
        if let Some(m) = have.iter().find(|i| !want.contains(i)) {
            return Err(format!("unneeded index item {m:?}"));
        }
        Ok(())
    }
}

// ----------------------------------------------------------------------
// Incremental defragmentation
// ----------------------------------------------------------------------

impl KmerTable {
    #[inline]
    pub(crate) fn mark_dirty(&mut self, pid: PieceId) {
        let w = pid as usize >> 6;
        if w >= self.dirty.len() {
            // Grow by ≥ 12.5 % (1 bit per piece).
            let n = (w + 1).max(self.dirty.len() + self.dirty.len() / 8);
            self.dirty.resize(n, 0);
        }
        self.dirty[w] |= 1 << (pid & 63);
    }

    /// Where the walk can continue from end `side` (0 = left, 1 = right) of
    /// piece `pid`, leaving the piece through that end: the piece whose end
    /// overlaps it by exactly k − 1 bases, and the end it is entered by.
    ///
    /// The outgoing (k−1)-mer `z` can only continue with the k-mers `z·c`;
    /// each is stored once, so a lookup tells where it is. The pieces join
    /// iff it is the first k-mer of a piece read forward (entered by its
    /// left end) or the last k-mer of a piece read backwards (entered by its
    /// right end). A k-mer found inside a piece means `z` is already followed
    /// there: no join through that base.
    fn next_piece(&self, pid: PieceId, side: u8) -> Option<(PieceId, u8)> {
        let k = self.cfg.k;
        let km1 = k - 1;
        let (g0, len) = self.store.view(pid);
        let z = if side == 1 {
            self.store.kmer_at_global(g0 + len as u64 - km1 as u64, km1)
        } else {
            revcomp_kmer(self.store.kmer_at_global(g0, km1), km1)
        };
        for c in 0..4u64 {
            let y = (z << 2) | c;
            let yr = revcomp_kmer(y, k);
            if y == yr {
                continue; // palindrome: orientation ambiguous
            }
            let Some(h) = self.locate_canonical(y.min(yr)) else {
                continue;
            };
            if h.piece == pid {
                continue;
            }
            let lq = self.store.length(h.piece) as usize;
            if h.kmer == y && h.pos == 0 {
                return Some((h.piece, 0));
            }
            if h.kmer == yr && h.pos as usize + k == lq {
                return Some((h.piece, 1));
            }
        }
        None
    }

    /// Reclaim what deletions freed without rebuilding the index:
    ///
    /// * slide the base array over the dead bases (pieces keep their ids;
    ///   only their starts move): one pass of 32-base word moves, ≈ 12 ms on
    ///   chr21 against 2 s for a full `defrag`;
    /// * drop unlinked overflow blocks, and re-size an index layer whose
    ///   slack reached [`INDEX_SHRINK_SLACK`] of it (O(entries) rebuild, no
    ///   minimizer windows recomputed).
    ///
    /// Pieces are not merged by default (see
    /// [`defrag_incremental_with`](Self::defrag_incremental_with)): the
    /// extra pieces and shared bases left by cuts
    /// ([`Fragmentation::piece_overhead`]) are for a full
    /// [`defrag`](Self::defrag).
    ///
    /// # Panics
    /// If a read is open (`insert_chunk` without `end_read`).
    pub fn defrag_incremental(&mut self) -> DefragReport {
        self.defrag_incremental_with(INCREMENTAL_MERGE_MEAN_LEN)
    }

    /// [`defrag_incremental`](Self::defrag_incremental), and before it, join
    /// the pieces created or cut since the last defragmentation with the
    /// pieces whose ends overlap them by k − 1 bases (found through the
    /// index: 8 lookups per piece), as long as the chain's pieces average at
    /// most `max_mean_len` bases (0: no merging; `usize::MAX`: any chain
    /// within `max_piece_len`). Each merged chain is rewritten at the end of
    /// the base array and its index items re-filed. Measured on chr21, a join
    /// cost ≈ 30 µs for ≈ 28 B saved, about 5× more per reclaimed byte than
    /// a full `defrag`, which is why merging is off by default.
    pub fn defrag_incremental_with(&mut self, max_mean_len: usize) -> DefragReport {
        self.assert_no_open_read("defrag_incremental");
        let k = self.cfg.k;
        let km1 = k - 1;
        let cap = self.cfg.max_piece_len as usize;
        let merge_mean = max_mean_len;
        let mut rep = DefragReport {
            pieces_before: self.store.len() as u64 - self.store.empty_pieces(),
            array_bases_before: self.store.total_bases(),
            entries_before: self.index.num_entries() + self.index2.num_entries(),
            table_bytes_before: self.table_bytes() as u64,
            ..Default::default()
        };
        // With merging off, the changed-piece marks are kept for a later
        // merging call (or cleared by a full defrag).
        let dirty = if merge_mean == 0 {
            Vec::new()
        } else {
            std::mem::take(&mut self.dirty)
        };
        let n_ids = self.store.len();
        let mut consumed = vec![0u64; n_ids.div_ceil(64)];
        let is_set = |v: &[u64], p: PieceId| {
            v.get(p as usize >> 6)
                .is_some_and(|w| w >> (p & 63) & 1 == 1)
        };
        let mut chain: Vec<(PieceId, bool)> = Vec::new();
        for (wi, &word) in dirty.iter().enumerate() {
            let mut bits = word;
            while bits != 0 {
                let pid = (wi * 64 + bits.trailing_zeros() as usize) as PieceId;
                bits &= bits - 1;
                if pid as usize >= n_ids || is_set(&consumed, pid) {
                    continue;
                }
                let len0 = self.store.length(pid) as usize;
                // Long pieces do not seed chains (they can still be joined
                // from a short neighbour): this keeps the 8 lookups per
                // piece for the pieces worth merging.
                if len0 < k || len0 > merge_mean.saturating_mul(2) {
                    continue;
                }
                // Grow a chain through both ends: (piece, reversed).
                chain.clear();
                chain.push((pid, false));
                let mut total = len0;
                for dir in [1u8, 0u8] {
                    loop {
                        // Leave the current end piece through its outer side.
                        let (p, rev) = if dir == 1 {
                            *chain.last().unwrap()
                        } else {
                            chain[0]
                        };
                        let exit = if (dir == 1) != rev { 1 } else { 0 };
                        let Some((np, enter)) = self.next_piece(p, exit) else {
                            break;
                        };
                        let nl = self.store.length(np) as usize;
                        let grown = total + nl - km1;
                        if is_set(&consumed, np)
                            || chain.iter().any(|&(q, _)| q == np)
                            || grown > cap
                            || grown > merge_mean.saturating_mul(chain.len() + 1)
                        {
                            break;
                        }
                        total += nl - km1;
                        // Entered by its left end: read forward in the walk
                        // direction; by its right end: reversed.
                        let walk_rev = enter == 1;
                        if dir == 1 {
                            chain.push((np, walk_rev));
                        } else {
                            chain.insert(0, (np, !walk_rev));
                        }
                    }
                }
                for &(p, _) in &chain {
                    let w = p as usize >> 6;
                    if w < consumed.len() {
                        consumed[w] |= 1 << (p & 63);
                    }
                }
                if chain.len() > 1 {
                    rep.joins += chain.len() as u64 - 1;
                    self.merge_chain(&chain, total);
                }
            }
        }
        // Slide the base array over dead bases; drop index orphans; give
        // back index head slots freed by deletions (O(entries) rebuild into
        // fewer buckets) once they are 1/16 of the index.
        let shared = if self.edits.dead_bases > 0 || rep.joins > 0 {
            self.store.compact_bases().1
        } else {
            self.edits.shared_bases
        };
        for ix in [&mut self.index, &mut self.index2] {
            if ix.orphans() > 0 {
                ix.compact();
            }
        }
        let (l, l2) = (self.cfg.l, self.l2);
        let store = &self.store;
        if index_slack_worth_shrinking(&self.index) {
            self.index.shrink_to_load(0.8, |p, q| {
                order_value(canonicalize(store.kmer_at(p, q, l), l))
            });
        }
        if index_slack_worth_shrinking(&self.index2) {
            self.index2.shrink_to_load(0.8, |p, q| {
                order_value(canonicalize(store.kmer_at(p, q, l2), l2))
            });
        }

        let e = &mut self.edits;
        e.incremental_defrags += 1;
        e.joins += rep.joins;
        e.dead_bases = 0;
        e.shared_bases = shared;
        rep.pieces_after = self.store.len() as u64 - self.store.empty_pieces();
        rep.array_bases_after = self.store.total_bases();
        rep.entries_after = self.index.num_entries() + self.index2.num_entries();
        rep.table_bytes_after = self.table_bytes() as u64;
        rep
    }

    /// Merge `chain` (pieces in walk order, `true` = read reversed; each
    /// overlapping the next by k − 1 bases) into one piece of `total` bases,
    /// written at the end of the base array under the first piece's id; the
    /// other ids are freed. The index items of every window are re-filed at
    /// the merged position (removed, then added once: items of the k − 1
    /// overlap are shared by consecutive pieces).
    fn merge_chain(&mut self, chain: &[(PieceId, bool)], total: usize) {
        let (k, l, l2) = (self.cfg.k, self.cfg.l, self.l2);
        let km1 = k - 1;
        let mut vals = Vec::new();
        let mut items = Vec::new();
        let mut new_items: Vec<Item> = Vec::new();
        let mut old_items: Vec<(Item, PieceId)> = Vec::new();
        let mut d = 0usize;
        for &(p, rev) in chain {
            let lp = self.store.length(p) as usize;
            items.clear();
            self.window_items(p, 0, lp - k, &mut vals, &mut items);
            for &(layer, key, q) in &items {
                let span = match layer {
                    1 => l,
                    2 => l2,
                    _ => k,
                };
                let qn = if rev {
                    d + lp - span - q as usize
                } else {
                    d + q as usize
                };
                new_items.push((layer, key, qn as u32));
                old_items.push(((layer, key, q), p));
            }
            d += lp - km1;
        }
        debug_assert_eq!(d + km1, total);
        new_items.sort_unstable();
        new_items.dedup();
        old_items.sort_unstable();
        old_items.dedup();

        // 1. Remove the chain's items (all copies).
        for &((layer, key, q), p) in &old_items {
            match layer {
                1 => {
                    self.index.remove(key, p, q);
                }
                2 => {
                    self.index2.remove(key, p, q);
                }
                _ => {
                    self.layer3.remove(key, (p, q));
                }
            }
        }
        // 2. Write the merged bases and re-point the ids.
        let g = self.store.end();
        for (i, &(p, rev)) in chain.iter().enumerate() {
            let (g0, lp) = self.store.view(p);
            let skip = if i == 0 { 0 } else { km1 as u64 };
            let lp = lp as u64;
            if rev {
                for j in (0..lp - skip).rev() {
                    let b = self.store.code_at_global(g0 + j);
                    self.store.push_code(3 - b);
                }
            } else {
                for j in skip..lp {
                    let b = self.store.code_at_global(g0 + j);
                    self.store.push_code(b);
                }
            }
        }
        let first = chain[0].0;
        for &(p, _) in &chain[1..] {
            let s = self.store.start(p);
            self.store.set_view(p, s, 0);
        }
        self.store.set_view(first, g, total as u32);
        // 3. File the merged items under the first id.
        let n1 = new_items.iter().filter(|it| it.0 == 1).count() as u64;
        let n2 = new_items.iter().filter(|it| it.0 == 2).count() as u64;
        let store = &self.store;
        self.index.reserve(n1, |p, q| {
            order_value(canonicalize(store.kmer_at(p, q, l), l))
        });
        self.index2.reserve(n2, |p, q| {
            order_value(canonicalize(store.kmer_at(p, q, l2), l2))
        });
        for &(layer, key, q) in &new_items {
            match layer {
                1 => {
                    self.index.insert(key, first, q);
                }
                2 => {
                    self.index2.insert(key, first, q);
                }
                _ => self.layer3.add(key, (first, q)),
            }
        }
        // The merged piece has exactly the chain's windows, so the number
        // of stored k-mers is unchanged.
    }
}
