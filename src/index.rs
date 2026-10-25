//! Minimizer index: flat, fixed-capacity buckets.
//!
//! An entry records one minimizer occurrence `(piece_id, q)` plus an 8-bit
//! tag. The minimizer value itself is not stored: a tag false-match is
//! rejected later by the exact k-mer comparison, so the tag only filters
//! work, never correctness.
//!
//! ## Layout
//!
//! All entries live in one flat `Vec<Block>`. A block is 80 bytes and holds
//! up to [`BLOCK_SLOTS`] = 8 entries in struct-of-arrays form (8 tags, 8
//! piece ids, 8 offsets), its fill count, and the index of an overflow block.
//! Blocks `0..n_buckets` are the bucket heads; overflow blocks are appended
//! after them. There is no per-bucket allocation, so there is no per-bucket
//! header: an entry costs 9 bytes plus its share of the block (80 / 8 = 10
//! bytes per slot).
//!
//! ## Placement
//!
//! `h = hash(m)`; `tag = low 8 bits of h`. Each minimizer has two candidate
//! buckets `b1 = range(h)` and `b2 = range(hash(h))` (multiply-shift range
//! reduction, so the bucket count need not be a power of two). A new entry
//! goes to the head with more free slots ("two choices"), which keeps heads
//! evenly filled. When both heads are full, the entry goes to the bucket's
//! first overflow block; a full overflow block gets a new one pushed in front
//! of it, so insertion is O(1) and only the first overflow block can be
//! partly filled. A lookup scans the chains of both `b1` and `b2`.
//!
//! ## Growth
//!
//! With fixed-capacity heads the bucket count must follow the number of
//! entries, otherwise memory is wasted (too many heads) or chains get long
//! (too few). When entries would exceed [`MAX_LOAD`] of the head slots, the
//! index grows by [`GROWTH`] and every entry is re-placed. Entries do not
//! store `m`, so the caller supplies a function recomputing `m` from
//! `(piece_id, q)` (the ℓ-mer at `q` in the piece).
//!
//! The old blocks are never kept alongside the new ones: their entries are
//! first written out as 8-byte `(piece_id, q)` records — to a temporary file
//! when a spill directory is configured, otherwise to a compact in-memory
//! list — the old blocks are freed, the new blocks allocated, and the records
//! streamed back in. Peak memory during growth is therefore the new index
//! plus at most 8 bytes per entry (nothing, with a spill file).
//!
//! ## Hot minimizers
//!
//! A minimizer can be declared *hot* ([`drain_minimizer`]): its entries are
//! removed and returned to the caller (which re-indexes those k-mers in a
//! fallback layer), and the minimizer is remembered so that lookups can be
//! routed. A flag bit in the head block of `b1` makes the common "not hot"
//! answer cost no extra memory access.
//!
//! [`drain_minimizer`]: MinimizerIndex::drain_minimizer

use crate::minimizer::mix64;
use crate::piece_store::PieceId;
use std::collections::HashSet;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::PathBuf;

/// Entries per block.
pub const BLOCK_SLOTS: usize = 8;
/// Fraction of head slots that may be used before the index grows.
pub const MAX_LOAD: f64 = 0.85;
/// Bucket-count growth factor.
pub const GROWTH: f64 = 1.5;

const NIL: u32 = u32::MAX;
const HOT_FLAG: u8 = 1;
const INDEX_SEED: u64 = 0x2545_f491_4f6c_dd1d;
const SECOND_SEED: u64 = 0x6a09_e667_f3bc_c909;

/// Fixed-capacity bucket block (80 bytes).
#[derive(Clone, Copy, Debug)]
#[repr(C)]
struct Block {
    tags: [u8; BLOCK_SLOTS],
    len: u8,
    /// Head blocks: bit 0 set if a hot minimizer has this bucket as `b1`.
    flags: u8,
    _pad: [u8; 2],
    next: u32,
    pids: [u32; BLOCK_SLOTS],
    qs: [u32; BLOCK_SLOTS],
}

const EMPTY: Block = Block {
    tags: [0; BLOCK_SLOTS],
    len: 0,
    flags: 0,
    _pad: [0; 2],
    next: NIL,
    pids: [0; BLOCK_SLOTS],
    qs: [0; BLOCK_SLOTS],
};

/// Size of one block in bytes.
pub const BLOCK_BYTES: usize = std::mem::size_of::<Block>();

impl Block {
    #[inline(always)]
    fn push(&mut self, tag: u8, pid: PieceId, q: u32) {
        let i = self.len as usize;
        debug_assert!(i < BLOCK_SLOTS);
        self.tags[i] = tag;
        self.pids[i] = pid;
        self.qs[i] = q;
        self.len += 1;
    }
}

/// One minimizer occurrence as seen through the public API.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Entry {
    pub tag: u8,
    pub piece_id: PieceId,
    pub q: u32,
}

/// Counters filled by a single lookup.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ScanTrace {
    /// Entries visited in the scanned chains.
    pub entries_scanned: u64,
    /// Entries whose tag matched (callback invocations).
    pub tag_matches: u64,
}

/// Where a minimizer's entries can be: two buckets and a tag.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Slot {
    pub b1: usize,
    pub b2: usize,
    pub tag: u8,
}

#[inline(always)]
fn range(h: u64, n: usize) -> usize {
    ((h as u128 * n as u128) >> 64) as usize
}

/// Bytes written to and read back from spill storage during growth.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SpillStats {
    pub rebuilds: u32,
    pub bytes_to_disk: u64,
    pub bytes_in_memory: u64,
}

#[derive(Clone, Debug)]
pub struct MinimizerIndex {
    /// Heads `0..n_buckets`, then overflow blocks.
    blocks: Vec<Block>,
    n_buckets: usize,
    n_entries: u64,
    grows: u32,
    /// Overflow blocks no longer reachable (after draining); dropped by
    /// [`compact`](Self::compact).
    orphans: usize,
    /// Hot minimizers (their entries were drained to a fallback layer).
    hot: HashSet<u64>,
    /// Where growth spills old entries (`None`: compact in-memory list).
    spill_dir: Option<PathBuf>,
    spill: SpillStats,
    /// Bytes per entry right after the last rebuild (growth or `defrag`);
    /// 0 = never rebuilt. See [`slack_bytes`](Self::slack_bytes).
    baseline_bpe: f64,
}

impl MinimizerIndex {
    /// Create an index with `2^bucket_bits` buckets initially
    /// (`bucket_bits ≤ 32`).
    pub fn new(bucket_bits: u32) -> Self {
        assert!(bucket_bits <= 32, "bucket_bits must be ≤ 32");
        Self::with_buckets(1usize << bucket_bits)
    }

    /// Create an index with `n` (≥ 1) buckets initially.
    pub fn with_buckets(n: usize) -> Self {
        let n = n.max(1);
        assert!(n < NIL as usize, "too many buckets");
        Self {
            blocks: vec![EMPTY; n],
            n_buckets: n,
            n_entries: 0,
            grows: 0,
            orphans: 0,
            hot: HashSet::new(),
            spill_dir: None,
            spill: SpillStats::default(),
            baseline_bpe: 0.0,
        }
    }

    /// Spill old entries to temporary files in `dir` during growth.
    pub fn with_spill_dir(mut self, dir: Option<PathBuf>) -> Self {
        self.spill_dir = dir;
        self
    }

    pub fn spill_stats(&self) -> SpillStats {
        self.spill
    }

    /// Buckets and tag of a minimizer value. Exposed for tests.
    #[inline(always)]
    pub fn slot(&self, m: u64) -> Slot {
        let h = mix64(m ^ INDEX_SEED);
        Slot {
            b1: range(h, self.n_buckets),
            b2: range(mix64(h ^ SECOND_SEED), self.n_buckets),
            tag: h as u8,
        }
    }

    /// Grow first if `additional` more entries would exceed the load limit.
    /// `remin(piece_id, q)` must return the minimizer value of an entry.
    pub fn reserve<F: Fn(PieceId, u32) -> u64>(&mut self, additional: u64, remin: F) {
        let limit = (MAX_LOAD * (BLOCK_SLOTS * self.n_buckets) as f64) as u64;
        if self.n_entries + additional <= limit {
            return;
        }
        let mut n = self.n_buckets;
        while ((self.n_entries + additional) as f64) > MAX_LOAD * (BLOCK_SLOTS * n) as f64 {
            n = ((n as f64 * GROWTH).ceil() as usize).max(n + 1);
        }
        self.rebuild(n, remin);
    }

    /// Re-place every entry into `n` buckets without holding the old and new
    /// blocks at the same time (see the module documentation).
    fn rebuild<F: Fn(PieceId, u32) -> u64>(&mut self, n: usize, remin: F) {
        let n_entries = self.n_entries;
        let old_overflow = self.overflow_blocks().saturating_sub(self.orphans);
        // 1. Write (piece_id, q) of every reachable entry, then free the blocks.
        let mut sink = Spill::new(self.spill_dir.as_ref(), n_entries);
        for b in 0..self.n_buckets {
            for blk in self.chain(b) {
                for i in 0..blk.len as usize {
                    sink.push(blk.pids[i], blk.qs[i]);
                }
            }
        }
        self.blocks = Vec::new();
        // 2. Allocate the new blocks and stream the entries back in.
        self.blocks = Vec::with_capacity(n + old_overflow / 2);
        self.blocks.resize(n, EMPTY);
        self.n_buckets = n;
        self.n_entries = 0;
        self.orphans = 0;
        self.grows += 1;
        self.spill.rebuilds += 1;
        let (disk, mem) = sink.finish(|pid, q| {
            let m = remin(pid, q);
            self.insert(m, pid, q);
        });
        self.spill.bytes_to_disk += disk;
        self.spill.bytes_in_memory += mem;
        debug_assert_eq!(self.n_entries, n_entries);
        // 3. Hot flags depend on the bucket count.
        let hot: Vec<u64> = self.hot.iter().copied().collect();
        for m in hot {
            let b1 = self.slot(m).b1;
            self.blocks[b1].flags |= HOT_FLAG;
        }
        self.compact();
        self.mark_baseline();
    }

    /// Re-place every entry into fewer buckets when the heads are below
    /// `load` (e.g. after many deletions), with the same spilled rebuild as
    /// growth; O(entries). Returns whether the index was rebuilt.
    pub fn shrink_to_load<F: Fn(PieceId, u32) -> u64>(&mut self, load: f64, remin: F) -> bool {
        let n = ((self.n_entries as f64 / (load * BLOCK_SLOTS as f64)).ceil() as usize).max(1);
        if n >= self.n_buckets {
            return false;
        }
        self.rebuild(n, remin);
        true
    }

    /// Record the current bytes per entry as the compact reference for
    /// [`slack_bytes`](Self::slack_bytes) (done after growth and `defrag`).
    pub fn mark_baseline(&mut self) {
        self.baseline_bpe = if self.n_entries == 0 {
            0.0
        } else {
            (self.blocks.len() * BLOCK_BYTES) as f64 / self.n_entries as f64
        };
    }

    /// Block bytes beyond what the current entries used per entry right after
    /// the last rebuild: slots freed by deletions, partly empty overflow
    /// blocks, unreclaimed orphans. 0 if the index was never rebuilt.
    pub fn slack_bytes(&self) -> u64 {
        if self.baseline_bpe == 0.0 {
            return 0;
        }
        let used = (self.blocks.len() * BLOCK_BYTES) as f64;
        (used - self.n_entries as f64 * self.baseline_bpe).max(0.0) as u64
    }

    /// Lay out every overflow chain contiguously, in chain order, right after
    /// the heads (bucket by bucket). Overflow blocks are otherwise allocated
    /// in insertion order, so a long chain (hot minimizer) would be scattered
    /// over the whole array and scanning it would chase pointers. O(blocks).
    pub fn compact(&mut self) {
        let n = self.n_buckets;
        if self.blocks.len() == n {
            self.orphans = 0;
            return;
        }
        // Heads stay in place; only the overflow area is rebuilt, so the
        // temporary copy is the size of the overflow area, not the index.
        let mut ov: Vec<Block> = Vec::with_capacity(self.blocks.len() - n - self.orphans);
        let mut links: Vec<(usize, u32)> = Vec::new(); // (head, first overflow)
        for b in 0..n {
            let first = self.blocks[b].next;
            if first == NIL {
                continue;
            }
            let start = (n + ov.len()) as u32;
            links.push((b, start));
            let mut cur = first;
            while cur != NIL {
                let mut blk = self.blocks[cur as usize];
                cur = blk.next;
                blk.next = if cur == NIL {
                    NIL
                } else {
                    (n + ov.len() + 1) as u32
                };
                ov.push(blk);
            }
        }
        debug_assert_eq!(ov.len(), self.blocks.len() - n - self.orphans);
        self.blocks.truncate(n);
        for (b, start) in links {
            self.blocks[b].next = start;
        }
        self.blocks.extend_from_slice(&ov);
        self.orphans = 0;
    }

    /// Add one entry (does not grow; call [`reserve`](Self::reserve) first).
    /// Returns `true` if the entry went to an overflow block (both heads were
    /// full), which is when the caller may want to check whether `m` is hot.
    pub fn insert(&mut self, m: u64, piece_id: PieceId, q: u32) -> bool {
        let overflowed;
        let s = self.slot(m);
        let (l1, l2) = (self.blocks[s.b1].len, self.blocks[s.b2].len);
        let head = if l2 < l1 { s.b2 } else { s.b1 };
        if (self.blocks[head].len as usize) < BLOCK_SLOTS {
            self.blocks[head].push(s.tag, piece_id, q);
            overflowed = false;
        } else {
            overflowed = true;
            // Both heads full: use a first overflow block with room, else
            // push a new block at the front of b1's chain.
            let pick = [s.b1, s.b2].into_iter().find(|&b| {
                let o = self.blocks[b].next;
                o != NIL && (self.blocks[o as usize].len as usize) < BLOCK_SLOTS
            });
            match pick {
                Some(b) => {
                    let o = self.blocks[b].next as usize;
                    self.blocks[o].push(s.tag, piece_id, q);
                }
                None => {
                    let b = s.b1;
                    if self.blocks.len() == self.blocks.capacity() {
                        // Grow the overflow area by ~12.5 % rather than 2×.
                        self.blocks.reserve_exact((self.blocks.len() / 8).max(16));
                    }
                    let o = self.blocks.len();
                    assert!(o < NIL as usize, "minimizer index full");
                    let mut nb = EMPTY;
                    nb.next = self.blocks[b].next;
                    nb.push(s.tag, piece_id, q);
                    self.blocks.push(nb);
                    self.blocks[b].next = o as u32;
                }
            }
        }
        self.n_entries += 1;
        overflowed
    }

    /// Remove every entry `(piece_id, q)` filed under minimizer `m`; returns
    /// how many were removed. Survivors of a block are kept packed; an
    /// overflow block that becomes empty is unlinked (and reclaimed by the
    /// next [`compact`](Self::compact)).
    pub fn remove(&mut self, m: u64, piece_id: PieceId, q: u32) -> usize {
        let s = self.slot(m);
        let mut n = self.remove_in_chain(s.b1, s.tag, piece_id, q);
        if s.b2 != s.b1 {
            n += self.remove_in_chain(s.b2, s.tag, piece_id, q);
        }
        self.n_entries -= n as u64;
        if self.orphans * 8 > self.blocks.len() {
            self.compact();
        }
        n
    }

    fn remove_in_chain(&mut self, head: usize, tag: u8, pid: PieceId, q: u32) -> usize {
        let mut n = 0;
        let mut prev = NIL;
        let mut cur = head as u32;
        while cur != NIL {
            let blk = &mut self.blocks[cur as usize];
            let mut i = 0;
            while i < blk.len as usize {
                if blk.tags[i] == tag && blk.pids[i] == pid && blk.qs[i] == q {
                    let last = blk.len as usize - 1;
                    blk.tags[i] = blk.tags[last];
                    blk.pids[i] = blk.pids[last];
                    blk.qs[i] = blk.qs[last];
                    blk.len -= 1;
                    n += 1;
                } else {
                    i += 1;
                }
            }
            let next = blk.next;
            if blk.len == 0 && prev != NIL {
                // Unlink an emptied overflow block (heads always stay).
                self.blocks[prev as usize].next = next;
                self.blocks[cur as usize] = EMPTY;
                self.orphans += 1;
            } else {
                prev = cur;
            }
            cur = next;
        }
        n
    }

    /// Change the entry `(piece_id, q)` filed under `m` into
    /// `(new_piece, new_q)` in place (the slot depends only on `m`). Returns
    /// whether the entry was found.
    pub fn rekey(
        &mut self,
        m: u64,
        piece_id: PieceId,
        q: u32,
        new_piece: PieceId,
        new_q: u32,
    ) -> bool {
        let s = self.slot(m);
        for head in [s.b1, s.b2] {
            let mut cur = head as u32;
            while cur != NIL {
                let blk = &mut self.blocks[cur as usize];
                for i in 0..blk.len as usize {
                    if blk.tags[i] == s.tag && blk.pids[i] == piece_id && blk.qs[i] == q {
                        blk.pids[i] = new_piece;
                        blk.qs[i] = new_q;
                        return true;
                    }
                }
                cur = blk.next;
            }
        }
        false
    }

    /// Free every block (entries are lost; hot minimizers are kept). The
    /// index must be refilled with [`clear_entries`](Self::clear_entries)
    /// before any other use.
    pub fn release(&mut self) {
        self.blocks = Vec::new();
        self.n_entries = 0;
        self.orphans = 0;
    }

    /// Drop every entry but keep the hot minimizers, and restart with `n`
    /// buckets (used to rebuild the index from scratch). Frees the old blocks
    /// before allocating the new heads.
    pub fn clear_entries(&mut self, n: usize) {
        let n = n.max(1);
        assert!(n < NIL as usize, "too many buckets");
        self.n_buckets = n;
        self.blocks = Vec::new();
        self.blocks = vec![EMPTY; n];
        self.n_entries = 0;
        self.orphans = 0;
        let hot: Vec<u64> = self.hot.iter().copied().collect();
        for m in hot {
            let b1 = self.slot(m).b1;
            self.blocks[b1].flags |= HOT_FLAG;
        }
    }

    /// Whether `m` has been declared hot. Costs one flag test on `b1`'s head
    /// (the block a lookup reads first anyway) unless that flag is set.
    #[inline(always)]
    pub fn is_hot(&self, m: u64) -> bool {
        let b1 = self.slot(m).b1;
        self.blocks[b1].flags & HOT_FLAG != 0 && self.hot.contains(&m)
    }

    /// Overflow blocks unlinked by deletions and not yet reclaimed.
    pub fn orphans(&self) -> usize {
        self.orphans
    }

    /// Number of hot minimizers.
    pub fn num_hot(&self) -> usize {
        self.hot.len()
    }

    /// The hot minimizers.
    pub fn hot_minimizers(&self) -> impl Iterator<Item = &u64> + '_ {
        self.hot.iter()
    }

    /// Number of entries whose minimizer is `m` (exact: tag matches are
    /// confirmed with `remin`).
    pub fn count_minimizer<F: Fn(PieceId, u32) -> u64>(&self, m: u64, remin: F) -> usize {
        let mut n = 0;
        self.for_each_candidate(m, &mut ScanTrace::default(), |p, q| {
            n += (remin(p, q) == m) as usize;
            false
        });
        n
    }

    /// Declare `m` hot: remove all its entries, remember it, and return the
    /// removed `(piece_id, q)` pairs.
    pub fn drain_minimizer<F: Fn(PieceId, u32) -> u64>(
        &mut self,
        m: u64,
        remin: F,
    ) -> Vec<(PieceId, u32)> {
        let s = self.slot(m);
        let mut drained = Vec::new();
        self.drain_chain(s.b1, s.tag, m, &remin, &mut drained);
        if s.b2 != s.b1 {
            self.drain_chain(s.b2, s.tag, m, &remin, &mut drained);
        }
        self.n_entries -= drained.len() as u64;
        self.hot.insert(m);
        self.blocks[s.b1].flags |= HOT_FLAG;
        // Reclaim unreachable blocks once they are a noticeable fraction.
        if self.orphans * 8 > self.blocks.len() {
            self.compact();
        }
        drained
    }

    /// Remove `m`'s entries from one chain and repack the survivors so that
    /// the head is filled first and only the first overflow block is partial.
    fn drain_chain<F: Fn(PieceId, u32) -> u64>(
        &mut self,
        head: usize,
        tag: u8,
        m: u64,
        remin: &F,
        drained: &mut Vec<(PieceId, u32)>,
    ) {
        let mut kept: Vec<(u8, PieceId, u32)> = Vec::new();
        let mut ov: Vec<usize> = Vec::new();
        let mut cur = head;
        loop {
            let blk = &self.blocks[cur];
            for i in 0..blk.len as usize {
                let e = (blk.tags[i], blk.pids[i], blk.qs[i]);
                if e.0 == tag && remin(e.1, e.2) == m {
                    drained.push((e.1, e.2));
                } else {
                    kept.push(e);
                }
            }
            if blk.next == NIL {
                break;
            }
            cur = blk.next as usize;
            ov.push(cur);
        }
        if kept.len() + drained.len() == 0 {
            return;
        }
        // Head: first BLOCK_SLOTS survivors.
        let n_head = kept.len().min(BLOCK_SLOTS);
        let rest = kept.len() - n_head;
        let need = rest.div_ceil(BLOCK_SLOTS);
        let first_len = if need == 0 {
            0
        } else {
            rest - BLOCK_SLOTS * (need - 1)
        };
        let mut it = kept.into_iter();
        let flags = self.blocks[head].flags;
        let mut h = EMPTY;
        h.flags = flags;
        for (t, p, q) in it.by_ref().take(n_head) {
            h.push(t, p, q);
        }
        h.next = if need > 0 { ov[0] as u32 } else { NIL };
        self.blocks[head] = h;
        for j in 0..need {
            let mut b = EMPTY;
            let take = if j == 0 { first_len } else { BLOCK_SLOTS };
            for (t, p, q) in it.by_ref().take(take) {
                b.push(t, p, q);
            }
            b.next = if j + 1 < need { ov[j + 1] as u32 } else { NIL };
            self.blocks[ov[j]] = b;
        }
        for &o in &ov[need..] {
            self.blocks[o] = EMPTY;
            self.orphans += 1;
        }
    }

    /// Visit every entry in `m`'s two chains whose tag matches `m`'s tag.
    /// The callback returns `true` to stop the scan early.
    /// Returns whether the scan was stopped by the callback.
    #[inline]
    pub fn for_each_candidate<F>(&self, m: u64, trace: &mut ScanTrace, mut f: F) -> bool
    where
        F: FnMut(PieceId, u32) -> bool,
    {
        let s = self.slot(m);
        if self.scan_chain(s.b1, s.tag, trace, &mut f) {
            return true;
        }
        s.b2 != s.b1 && self.scan_chain(s.b2, s.tag, trace, &mut f)
    }

    #[inline(always)]
    fn scan_chain<F>(&self, head: usize, tag: u8, trace: &mut ScanTrace, f: &mut F) -> bool
    where
        F: FnMut(PieceId, u32) -> bool,
    {
        let mut b = head;
        loop {
            let blk = &self.blocks[b];
            let n = blk.len as usize;
            trace.entries_scanned += n as u64;
            for i in 0..n {
                if blk.tags[i] == tag {
                    trace.tag_matches += 1;
                    if f(blk.pids[i], blk.qs[i]) {
                        return true;
                    }
                }
            }
            if blk.next == NIL {
                return false;
            }
            b = blk.next as usize;
        }
    }

    /// Append every tag-matching entry of `m`'s chains to `out` as
    /// `(piece_id, q)` (used to cache the candidates of a minimizer).
    pub fn collect_candidates(&self, m: u64, trace: &mut ScanTrace, out: &mut Vec<(PieceId, u32)>) {
        self.for_each_candidate(m, trace, |p, q| {
            out.push((p, q));
            false
        });
    }

    /// Current number of buckets.
    pub fn num_buckets(&self) -> usize {
        self.n_buckets
    }

    /// Overflow blocks allocated (including not yet reclaimed orphans).
    pub fn overflow_blocks(&self) -> usize {
        self.blocks.len() - self.n_buckets
    }

    /// Number of times the index has grown.
    pub fn grows(&self) -> u32 {
        self.grows
    }

    pub fn num_entries(&self) -> u64 {
        self.n_entries
    }

    /// Entries stored in bucket `b`'s chain.
    pub fn bucket_len(&self, b: usize) -> usize {
        self.chain(b).map(|blk| blk.len as usize).sum()
    }

    fn chain(&self, head: usize) -> impl Iterator<Item = &Block> + '_ {
        let mut cur = head as u32;
        std::iter::from_fn(move || {
            if cur == NIL {
                return None;
            }
            let blk = &self.blocks[cur as usize];
            cur = blk.next;
            Some(blk)
        })
    }

    /// Every entry with the bucket whose chain holds it.
    pub fn iter_entries_with_bucket(&self) -> impl Iterator<Item = (usize, Entry)> + '_ {
        (0..self.n_buckets).flat_map(move |b| {
            self.chain(b).flat_map(move |blk| {
                (0..blk.len as usize).map(move |i| {
                    (
                        b,
                        Entry {
                            tag: blk.tags[i],
                            piece_id: blk.pids[i],
                            q: blk.qs[i],
                        },
                    )
                })
            })
        })
    }

    /// Every entry, bucket by bucket.
    pub fn iter_entries(&self) -> impl Iterator<Item = Entry> + '_ {
        self.iter_entries_with_bucket().map(|(_, e)| e)
    }

    /// Heap bytes in use (by length) and reserved (by capacity), including
    /// the set of hot minimizers.
    pub fn heap_bytes(&self) -> (usize, usize) {
        let hot = hashset_u64_bytes(self.hot.capacity());
        (
            self.blocks.len() * BLOCK_BYTES + hot,
            self.blocks.capacity() * BLOCK_BYTES + hot,
        )
    }

    /// Compact the overflow chains and release unused capacity.
    pub fn shrink_to_fit(&mut self) {
        self.compact();
        self.blocks.shrink_to_fit();
    }
}

/// Heap bytes of a std `HashSet<u64>` with the given capacity (SwissTable:
/// an 8-byte slot plus one control byte per bucket, plus 16 control bytes).
pub fn hashset_u64_bytes(cap: usize) -> usize {
    if cap == 0 {
        return 0;
    }
    let buckets = if cap < 8 {
        if cap < 4 {
            4
        } else {
            8
        }
    } else {
        (cap * 8 / 7).next_power_of_two()
    };
    buckets * 9 + 16
}

/// Temporary storage for `(piece_id, q)` records during a rebuild.
enum Spill {
    Memory(Vec<(PieceId, u32)>),
    Disk {
        path: PathBuf,
        out: BufWriter<std::fs::File>,
        n: u64,
    },
}

static SPILL_COUNTER: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

impl Spill {
    fn new(dir: Option<&PathBuf>, n: u64) -> Self {
        if let Some(dir) = dir {
            let id = SPILL_COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            let path = dir.join(format!("unihash-spill-{}-{id}.bin", std::process::id()));
            match std::fs::File::create(&path) {
                Ok(f) => {
                    return Spill::Disk {
                        path,
                        out: BufWriter::with_capacity(1 << 20, f),
                        n: 0,
                    }
                }
                Err(e) => eprintln!(
                    "unihash: cannot create spill file {} ({e}); spilling to memory",
                    path.display()
                ),
            }
        }
        Spill::Memory(Vec::with_capacity(n as usize))
    }

    fn push(&mut self, pid: PieceId, q: u32) {
        match self {
            Spill::Memory(v) => v.push((pid, q)),
            Spill::Disk { out, n, .. } => {
                let mut rec = [0u8; 8];
                rec[..4].copy_from_slice(&pid.to_le_bytes());
                rec[4..].copy_from_slice(&q.to_le_bytes());
                out.write_all(&rec).expect("unihash: spill write failed");
                *n += 1;
            }
        }
    }

    /// Replay the records into `f`; returns (bytes via disk, bytes via memory).
    fn finish<F: FnMut(PieceId, u32)>(self, mut f: F) -> (u64, u64) {
        match self {
            Spill::Memory(v) => {
                let bytes = (v.len() * 8) as u64;
                for (p, q) in v {
                    f(p, q);
                }
                (0, bytes)
            }
            Spill::Disk { path, out, n } => {
                let file = out.into_inner().expect("unihash: spill flush failed");
                drop(file);
                let mut r = BufReader::with_capacity(
                    1 << 20,
                    std::fs::File::open(&path).expect("unihash: spill reopen failed"),
                );
                let mut rec = [0u8; 8];
                for _ in 0..n {
                    r.read_exact(&mut rec).expect("unihash: spill read failed");
                    f(
                        u32::from_le_bytes(rec[..4].try_into().unwrap()),
                        u32::from_le_bytes(rec[4..].try_into().unwrap()),
                    );
                }
                let _ = std::fs::remove_file(&path);
                (n * 8, 0)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn block_is_80_bytes() {
        assert_eq!(BLOCK_BYTES, 80);
    }
}
