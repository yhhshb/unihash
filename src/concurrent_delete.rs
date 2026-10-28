//! Deletion in a [`ConcurrentKmerTable`], concurrent with insertion,
//! queries and other deletions.
//!
//! A cut works like [`KmerTable::delete_read`]'s (the left part keeps the
//! piece id, the right part is a new view on the same bases, every surviving
//! k-mer keeps its `k − 1` flanks), but in an order that keeps every
//! surviving k-mer findable at every moment:
//!
//! 1. the right part is published under a **new** id and its index items are
//!    added (its k-mers are now reachable through both pieces);
//! 2. the cut piece's length shrinks — the moment the deleted k-mers
//!    disappear;
//! 3. index items only the cut piece needed are removed (until then they
//!    point past its end or at windows that still exist, which a lookup
//!    verifies against the bases: never a wrong answer).
//!
//! Cuts of one piece are serialized by striped piece locks, and a cut is
//! validated under its lock (the piece still starts where the k-mers were
//! found and still covers them); a stale cut is retried from a fresh lookup.
//! A piece still growing in a writer's open read is never cut: the deletion
//! waits until that read's run closes. A right part is flagged the same way
//! until all its items are added, so no other cut works on it half-filed.
//! Each cut holds the table's gate
//! (shared), during which no minimizer becomes hot, so the index items it
//! computes stay valid; it takes at most one shard lock at a time, after the
//! piece lock (writers never take piece locks): no deadlock.

use super::*;
use crate::codec::{canonical, canonical_kmers, canonicalize, encode, RollingKmer, INVALID};

/// Index item: (layer, key, offset in the piece).
type Item = (u8, u64, u32);

enum Cut {
    Done,
    /// The piece changed since the k-mers were found: look them up again.
    Stale,
    /// The piece is a writer's open run: wait for it to close.
    Open,
}

/// Successive k-mers of a query found at successive windows of one piece.
#[derive(Clone, Copy)]
struct DelRun {
    pid: PieceId,
    /// Piece start and length when the run was found.
    g0: u64,
    len: usize,
    /// First and last window (decreasing if the piece holds the reverse
    /// complement of the query).
    first: usize,
    last: usize,
    same: bool,
    /// Query offset of the first k-mer, and number of k-mers.
    r0: usize,
    n: usize,
}

/// Wait before retrying a cut whose piece was open for the `n`-th time:
/// spin, then yield, then sleep (do not take the CPU from the thread that
/// will close the piece, possibly the other hardware thread of this core).
fn backoff(n: u32, c: &mut ContentionStats) {
    let t0 = Instant::now();
    if n <= 4 {
        for _ in 0..64 << n {
            std::hint::spin_loop();
        }
    } else if n <= 8 {
        std::thread::yield_now();
    } else {
        std::thread::sleep(std::time::Duration::from_micros(50));
    }
    c.open_wait_ns += t0.elapsed().as_nanos() as u64;
}

/// Counters of one deletion call, merged into the table when it returns.
#[derive(Default)]
struct Del {
    edits: EditStats,
    cont: ContentionStats,
    changed: Vec<PieceId>,
}

impl ConcurrentKmerTable {
    /// Delete one k-mer (either orientation). Returns whether this call
    /// removed it. Safe to call while other threads insert, query or delete.
    ///
    /// If the k-mer is in a read that a [`Writer`] has not finished
    /// (`insert_chunk` without `end_read`) and its run is still open, this
    /// waits until it closes: do not delete from a thread whose own writer
    /// has such a read open.
    pub fn delete(&self, kmer: u64) -> bool {
        if kmer & !self.mask != 0 {
            return false;
        }
        let mut d = Del::default();
        d.edits.delete_calls += 1;
        let found = self.delete_one(codec::canonicalize(kmer, self.cfg.k), &mut d);
        self.merge_del(d);
        found
    }

    /// Delete every k-mer of `read` (split at non-ACGT characters like
    /// insertion; absent k-mers are ignored). Successive k-mers found at
    /// successive positions of a piece are cut as one run. Returns the
    /// number of stored k-mer occurrences this call removed. Safe to call
    /// while other threads insert, query or delete; see
    /// [`delete`](Self::delete) about open reads.
    pub fn delete_read(&self, read: &[u8]) -> u64 {
        let mut d = Del::default();
        d.edits.delete_calls += 1;
        // Ranges of `read` still to process (re-queued when a cut is stale,
        // or open, with the number of times it found its piece open).
        let mut pending = vec![(0, read.len(), 0)];
        while let Some((a, b, opens)) = pending.pop() {
            if opens > 0 {
                backoff(opens, &mut d.cont);
            }
            self.delete_pass(read, a, b, opens, &mut d, &mut pending);
        }
        if self.may_have_duplicates {
            // Copies of a k-mer elsewhere in its piece (tables built with
            // deferred run indexing) survive the pass.
            let mut again = false;
            for x in canonical_kmers(read, self.cfg.k) {
                if self
                    .locate_canonical(x, &mut QueryTrace::default())
                    .is_some()
                {
                    again = true;
                    self.delete_one(x, &mut d);
                }
            }
            d.edits.recheck_passes += again as u64;
        }
        let n = d.edits.kmers_deleted;
        self.merge_del(d);
        n
    }

    /// Delete canonical `x` (every copy if the table may hold duplicates).
    fn delete_one(&self, x: u64, d: &mut Del) -> bool {
        let k = self.cfg.k;
        let mut found = false;
        let mut opens = 0;
        while let Some(h) = self.locate_canonical(x, &mut QueryTrace::default()) {
            let g0 = self.store.view(h.piece).0;
            // The id may have been reused since the lookup: check the bases.
            if canonicalize(self.store.kmer_at_global(g0 + h.pos as u64, k), k) != x {
                continue;
            }
            let p = h.pos as usize;
            match self.cut(h.piece, g0, p, p, d) {
                Cut::Done => {
                    found = true;
                    if !self.may_have_duplicates {
                        break;
                    }
                }
                Cut::Stale => d.cont.stale_retries += 1,
                Cut::Open => {
                    opens += 1;
                    d.cont.open_waits += 1;
                    backoff(opens, &mut d.cont);
                }
            }
        }
        found
    }

    /// One deletion pass over `read[a..b]`; ranges whose cut must be retried
    /// are pushed to `pending`.
    fn delete_pass(
        &self,
        read: &[u8],
        a: usize,
        b: usize,
        opens: u32,
        d: &mut Del,
        pending: &mut Vec<(usize, usize, u32)>,
    ) {
        let k = self.cfg.k;
        let mut roll = RollingKmer::new(k);
        let mut run: Option<DelRun> = None;
        let finish = |r: DelRun, d: &mut Del, pending: &mut Vec<(usize, usize, u32)>| {
            let (lo, hi) = (r.first.min(r.last), r.first.max(r.last));
            match self.cut(r.pid, r.g0, lo, hi, d) {
                Cut::Done => {}
                Cut::Stale => {
                    d.cont.stale_retries += 1;
                    pending.push((r.r0, r.r0 + r.n + k - 1, 0));
                }
                Cut::Open => {
                    d.cont.open_waits += 1;
                    pending.push((r.r0, r.r0 + r.n + k - 1, opens + 1));
                }
            }
        };
        for (i, &c) in read.iter().enumerate().take(b).skip(a) {
            let bc = encode(c);
            if bc == INVALID {
                if let Some(r) = run.take() {
                    finish(r, d, pending);
                }
                roll.reset();
                continue;
            }
            if !roll.push(bc) {
                continue;
            }
            if let Some(r) = run.as_mut() {
                // Does the next k-mer sit right after the run in the piece?
                let ok = if r.same {
                    r.last + k < r.len
                        && self.store.code_at_global(r.g0 + (r.last + k) as u64) == bc
                } else {
                    r.last >= 1 && self.store.code_at_global(r.g0 + r.last as u64 - 1) == 3 - bc
                };
                if ok {
                    r.last = if r.same { r.last + 1 } else { r.last - 1 };
                    r.n += 1;
                    continue;
                }
                let r = run.take().expect("run");
                finish(r, d, pending);
            }
            let (fwd, rc) = (roll.fwd(), roll.rc());
            let x = canonical(fwd, rc);
            if let Some(h) = self.locate_canonical(x, &mut QueryTrace::default()) {
                let (g0, len) = self.store.view(h.piece);
                let f = self.store.kmer_at_global(g0 + h.pos as u64, k);
                if (h.pos as usize) + k <= len as usize && canonicalize(f, k) == x {
                    let p = h.pos as usize;
                    run = Some(DelRun {
                        pid: h.piece,
                        g0,
                        len: len as usize,
                        first: p,
                        last: p,
                        same: f == fwd,
                        r0: i + 1 - k,
                        n: 1,
                    });
                } else {
                    // Changed under us: retry this k-mer.
                    d.cont.stale_retries += 1;
                    pending.push((i + 1 - k, i + 1, 0));
                }
            }
        }
        if let Some(r) = run {
            finish(r, d, pending);
        }
    }

    /// Delete windows `lo..=hi` of piece `pid`, found at start `g0` (see the
    /// module documentation).
    fn cut(&self, pid: PieceId, g0: u64, lo: usize, hi: usize, d: &mut Del) -> Cut {
        let (k, l) = (self.cfg.k, self.cfg.l);
        let span = k - l;
        let _gate = rlock(&self.gate);
        let piece_lock = &self.piece_locks[shard_of(pid as u64, PIECE_LOCKS)];
        let _piece = match piece_lock.try_lock() {
            Ok(g) => g,
            Err(std::sync::TryLockError::Poisoned(p)) => p.into_inner(),
            Err(std::sync::TryLockError::WouldBlock) => {
                d.cont.piece_contended += 1;
                let t0 = Instant::now();
                let g = lock(piece_lock);
                d.cont.piece_wait_ns += t0.elapsed().as_nanos() as u64;
                g
            }
        };
        let (g, len) = self.store.view(pid);
        let len = len as usize;
        if g != g0 || hi + k > len || lo > hi {
            return Cut::Stale;
        }
        if self.store.is_open(pid) {
            return Cut::Open;
        }
        let nw = len - k + 1;

        // Items of the windows that lose their piece (cut run and its left
        // neighbourhood), of the right part, and those the left part keeps.
        let (mut cut, mut moved, mut keep) = (Vec::new(), Vec::new(), Vec::new());
        {
            let hot1 = rlock(&self.hot1);
            let hot2 = rlock(&self.hot2);
            let items = |from: usize, to: usize, out: &mut Vec<Item>| {
                self.window_items(g0, from, to, &hot1, &hot2, out);
                out.sort_unstable();
                out.dedup();
            };
            items(lo.saturating_sub(span), hi, &mut cut);
            if hi + 1 < nw {
                items(hi + 1, nw - 1, &mut moved);
            }
            if lo > 0 {
                items(lo.saturating_sub(2 * span), lo - 1, &mut keep);
            }
        }

        // 1. The right part under a new id, with its items.
        let a_len = if lo > 0 { lo + k - 1 } else { 0 };
        let b_off = hi + 1;
        let b_len = if hi + 1 < nw { len - b_off } else { 0 };
        // While its items are being added, the right part is flagged open:
        // a concurrent cut that finds it through the first items waits.
        let right = (b_len > 0).then(|| {
            let b = self.store.alloc_pid();
            self.store.set_open(b);
            self.store.set_view(b, g0 + b_off as u64, b_len as u32);
            b
        });
        if let Some(b) = right {
            let off = b_off as u32;
            let adds: Vec<Item> = moved
                .iter()
                .map(|&(ly, key, q)| (ly, key, q - off))
                .collect();
            self.apply(&adds, b, true, d);
            self.store.clear_open(b);
            d.edits.items_moved += moved.len() as u64;
        }

        // 2. The deleted k-mers disappear.
        self.store.set_view(pid, g0, a_len as u32);

        // 3. Items only the old layout needed.
        let mut removes: Vec<Item> = cut
            .iter()
            .chain(&moved)
            .filter(|it| keep.binary_search(it).is_err())
            .copied()
            .collect();
        removes.sort_unstable();
        removes.dedup();
        self.apply(&removes, pid, false, d);
        if a_len == 0 {
            self.store.free_pid(pid);
        } else {
            d.changed.push(pid);
        }
        d.changed.extend(right);

        // Accounting (as KmerTable::delete_windows).
        let both = a_len > 0 && b_len > 0;
        let overlap = if both { a_len.saturating_sub(b_off) } else { 0 };
        let union = a_len + b_len - overlap;
        let e = &mut d.edits;
        e.dead_bases += (len - union) as u64;
        e.shared_bases += overlap as u64;
        e.kmers_deleted += (hi - lo + 1) as u64;
        e.cuts += 1;
        e.cuts_since_defrag += 1;
        e.splits += both as u64;
        Cut::Done
    }

    /// Add (or remove) index items for piece `pid`, one shard lock at a
    /// time.
    fn apply(&self, items: &[Item], pid: PieceId, add: bool, d: &mut Del) {
        let (n1, n2) = (self.l1.len(), self.l2s.len());
        let shard = |&(ly, key, _): &Item| match ly {
            1 => shard_of(key, n1),
            2 => shard_of(key, n2),
            _ => shard_of(key, LAYER3_SHARDS),
        };
        let mut order: Vec<(u8, usize, u64, u32)> = items
            .iter()
            .map(|it| (it.0, shard(it), it.1, it.2))
            .collect();
        order.sort_unstable();
        let mut i = 0;
        while i < order.len() {
            let (ly, sh) = (order[i].0, order[i].1);
            let j = i + order[i..]
                .iter()
                .take_while(|o| o.0 == ly && o.1 == sh)
                .count();
            let group = &order[i..j];
            i = j;
            match ly {
                1 | 2 => {
                    let (lock, l) = if ly == 1 {
                        (&self.l1[sh], self.cfg.l)
                    } else {
                        (&self.l2s[sh], self.l2)
                    };
                    let mut g = acquire(lock, &mut d.cont);
                    if add {
                        reserve_timed(&mut g, group.len() as u64, self.remin(l), &mut d.cont);
                        for &(_, _, key, q) in group {
                            g.insert(key, pid, q);
                        }
                    } else {
                        for &(_, _, key, q) in group {
                            let n = g.remove(key, pid, q);
                            debug_assert!(
                                n > 0 || ly == 2,
                                "index item ({ly}, {key}, {q}) missing"
                            );
                        }
                    }
                }
                _ => {
                    let mut g = acquire(&self.l3[sh], &mut d.cont);
                    for &(_, _, key, q) in group {
                        if add {
                            g.add(key, (pid, q));
                        } else {
                            g.remove(key, (pid, q));
                        }
                    }
                }
            }
        }
    }

    /// Index items of windows `s_from..=s_to` of the piece starting at
    /// global offset `g0` (see `KmerTable::window_items`).
    fn window_items(
        &self,
        g0: u64,
        s_from: usize,
        s_to: usize,
        hot1: &HashSet<u64>,
        hot2: &HashSet<u64>,
        out: &mut Vec<Item>,
    ) {
        let (k, l, w, l2) = (self.cfg.k, self.cfg.l, self.w, self.l2);
        if s_from > s_to {
            return;
        }
        let mut vals = Vec::with_capacity(s_to - s_from + k);
        let mut roll = RollingKmer::new(l);
        for g in g0 + s_from as u64..g0 + (s_to + k) as u64 {
            if roll.push(self.store.code_at_global(g)) {
                vals.push(order_value(roll.canonical()));
            }
        }
        let mut wm = WindowMin::new(w);
        for (i, &v) in vals.iter().enumerate() {
            wm.push(v);
            let Some((m, ties)) = wm.current() else {
                continue;
            };
            let s = s_from + i + 1 - w;
            if l2 > 0 && hot1.contains(&m) {
                let f = self.store.kmer_at_global(g0 + s as u64, k);
                let (m2, t2) = kmer_minimizer_mask(f, k, l2);
                if hot2.contains(&m2) {
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

    fn merge_del(&self, d: Del) {
        let mut c = lock(&self.carry);
        let (e, o) = (&mut c.edits, &d.edits);
        e.delete_calls += o.delete_calls;
        e.kmers_deleted += o.kmers_deleted;
        e.cuts += o.cuts;
        e.splits += o.splits;
        e.items_moved += o.items_moved;
        e.recheck_passes += o.recheck_passes;
        e.dead_bases += o.dead_bases;
        e.shared_bases += o.shared_bases;
        e.cuts_since_defrag += o.cuts_since_defrag;
        // Writers merge their counts when they finish: wrap until then.
        c.stored_kmers = c.stored_kmers.wrapping_sub(o.kmers_deleted);
        c.contention.merge(&d.cont);
        c.changed.extend(d.changed);
    }

    // ------------------------------------------------------------------
    // Exclusive maintenance
    // ------------------------------------------------------------------

    /// Run `f` on the content as a [`KmerTable`] (converted and back, keeping
    /// the shard count), e.g. for [`KmerTable::defrag`]. O(table).
    pub fn with_table<R>(&mut self, f: impl FnOnce(&mut KmerTable) -> R) -> R {
        let shards = self.l1.len();
        let empty = Self::with_shards(self.cfg.clone(), 1).expect("validated configuration");
        let mut t = std::mem::replace(self, empty).into_table();
        let r = f(&mut t);
        *self = Self::from_table_with_shards(t, shards);
        r
    }

    /// [`KmerTable::defrag`] (merge overlapping pieces, drop dead space).
    pub fn defrag(&mut self) -> crate::edit::DefragReport {
        self.with_table(KmerTable::defrag)
    }

    /// [`KmerTable::defrag_incremental`].
    pub fn defrag_incremental(&mut self) -> crate::edit::DefragReport {
        self.with_table(KmerTable::defrag_incremental)
    }

    /// [`KmerTable::fragmentation`] estimate from the deletion counters
    /// (without conversion): dead bases, shared bases, cuts since the last
    /// defragmentation.
    pub fn edit_stats(&self) -> EditStats {
        lock(&self.carry).edits
    }
}
