//! Defragmentation of a [`ConcurrentKmerTable`] in place (exclusive access),
//! without converting to a [`KmerTable`] and back.
//!
//! * [`defrag`](ConcurrentKmerTable::defrag): the same chains as
//!   [`KmerTable::defrag`] (pieces whose ends overlap by `k − 1` bases, see
//!   `plan_chains`) are written into a new base array, pieces renumbered in
//!   chain order, and every piece re-indexed into the shards by all cores.
//!   The hot minimizers stay hot.
//! * [`defrag_incremental`](ConcurrentKmerTable::defrag_incremental): the
//!   base array is rewritten without dead bases (piece ids and index entries
//!   unchanged, only the views move), unreachable overflow blocks are
//!   dropped, and shards whose deletions freed enough head slots are
//!   re-sized.

use super::delete::Del;
use super::*;
use crate::edit::{plan_chains, DefragReport, INDEX_SHRINK_SLACK};

/// Writes bases one after another into a fresh [`SharedStore`] (single
/// thread, before the store is shared).
struct SeqWriter<'s> {
    store: &'s SharedStore,
    /// Next global offset, and the word being filled.
    g: u64,
    acc: u64,
}

impl<'s> SeqWriter<'s> {
    fn new(store: &'s SharedStore) -> Self {
        Self {
            store,
            g: 0,
            acc: 0,
        }
    }

    #[inline(always)]
    fn push(&mut self, code: u8) {
        self.acc |= (code as u64 & 3) << (62 - 2 * (self.g & 31));
        self.g += 1;
        if self.g & 31 == 0 {
            self.store
                .data
                .slot((self.g >> 5) - 1)
                .store(self.acc, Ordering::Relaxed);
            self.acc = 0;
        }
    }

    /// Append 32 bases packed in `w` (first base in the top bits).
    #[inline(always)]
    fn push32(&mut self, w: u64) {
        let r = (self.g & 31) as u32;
        let word = |me: &Self, i: u64, v: u64| me.store.data.slot(i).store(v, Ordering::Relaxed);
        if r == 0 {
            word(self, self.g >> 5, w);
        } else {
            word(self, self.g >> 5, self.acc | w >> (2 * r));
            self.acc = w << (64 - 2 * r);
        }
        self.g += 32;
    }

    /// Copy `len` bases starting at global offset `from` of `src`, forward
    /// or reverse-complemented (read from the end), 32 at a time.
    fn copy(&mut self, src: &SharedStore, from: u64, len: u64, rev: bool) {
        let mut i = 0;
        while len - i >= 32 {
            let w = if rev {
                codec::revcomp_kmer(src.kmer_at_global(from + len - i - 32, 32), 32)
            } else {
                src.kmer_at_global(from + i, 32)
            };
            self.push32(w);
            i += 32;
        }
        while i < len {
            let c = if rev {
                3 - src.code_at_global(from + len - 1 - i)
            } else {
                src.code_at_global(from + i)
            };
            self.push(c);
            i += 1;
        }
    }

    /// Store the last partial word; returns the number of bases written.
    fn finish(self) -> u64 {
        if self.g & 31 != 0 {
            self.store
                .data
                .slot(self.g >> 5)
                .store(self.acc, Ordering::Relaxed);
        }
        self.g
    }
}

impl ConcurrentKmerTable {
    /// Heap bytes of pieces, index shards and layer 3 (by capacity), for
    /// defragmentation reports.
    fn table_heap_bytes(&self) -> u64 {
        self.stats().heap_bytes as u64
    }

    fn layer_entries(&self) -> (u64, u64) {
        let count = |s: &[Pad<RwLock<MinimizerIndex>>]| -> u64 {
            s.iter().map(|x| rlock(x).num_entries()).sum()
        };
        (count(&self.l1), count(&self.l2s))
    }

    /// Bases between the first and the last live base ever handed out.
    fn array_bases(&self) -> u64 {
        let n = self.store.n_pids.load(Ordering::Relaxed);
        (0..n)
            .map(|p| {
                let (s, l) = self.store.view(p as PieceId);
                if l > 0 {
                    s + l as u64
                } else {
                    0
                }
            })
            .max()
            .unwrap_or(0)
    }

    fn report_before(&self) -> DefragReport {
        let n = self.store.n_pids.load(Ordering::Relaxed);
        let (e1, e2) = self.layer_entries();
        DefragReport {
            pieces_before: (0..n)
                .filter(|&p| self.store.view(p as PieceId).1 > 0)
                .count() as u64,
            array_bases_before: self.array_bases(),
            entries_before: e1 + e2,
            table_bytes_before: self.table_heap_bytes(),
            ..Default::default()
        }
    }

    fn report_after(&self, mut rep: DefragReport) -> DefragReport {
        let n = self.store.n_pids.load(Ordering::Relaxed);
        let (e1, e2) = self.layer_entries();
        rep.pieces_after = (0..n)
            .filter(|&p| self.store.view(p as PieceId).1 > 0)
            .count() as u64;
        rep.array_bases_after = self.array_bases();
        rep.entries_after = e1 + e2;
        rep.table_bytes_after = self.table_heap_bytes();
        rep
    }

    /// [`KmerTable::defrag`] on the concurrent table, in place: merge
    /// pieces whose ends overlap by `k − 1` bases, drop dead space and
    /// emptied pieces, renumber pieces, and re-index them (all cores). Gives
    /// the same pieces as `into_table`, `defrag`, `from_table`, without the
    /// conversions. Answers are unchanged.
    pub fn defrag(&mut self) -> DefragReport {
        let mut rep = self.report_before();
        let (k, cap) = (self.cfg.k, self.cfg.max_piece_len as u64);
        let (n1, n2) = self.layer_entries();

        // 0. Free the index shards and layer 3 (rebuilt from the pieces).
        let spill = self.cfg.spill_dir.clone();
        for s in self.l1.iter_mut().chain(self.l2s.iter_mut()) {
            *s.get_mut().unwrap_or_else(PoisonError::into_inner) =
                MinimizerIndex::with_buckets(1).with_spill_dir(spill.clone());
        }
        for s in self.l3.iter_mut() {
            s.get_mut().unwrap_or_else(PoisonError::into_inner).clear();
        }

        // 1-2. Chains, as KmerTable::defrag.
        let n = self.store.n_pids.load(Ordering::Relaxed) as usize;
        let chains = {
            let st = &self.store;
            plan_chains(
                n,
                k,
                cap,
                |p| st.view(p),
                |g, len| st.kmer_at_global(g, len),
            )
        };
        rep.joins = chains.joins;

        // 3. Write the chains into a new base array; piece i = chain i.
        let store = SharedStore::new();
        let mut stored = 0u64;
        {
            let old = &self.store;
            let mut w = SeqWriter::new(&store);
            for (i, chain) in chains.iter().enumerate() {
                let g = w.g;
                for (j, &(pid, rev)) in chain.iter().enumerate() {
                    let skip = if j == 0 { 0 } else { (k - 1) as u64 };
                    let (g0, len) = old.view(pid);
                    let len = len as u64;
                    if rev {
                        w.copy(old, g0, len - skip, true);
                    } else {
                        w.copy(old, g0 + skip, len - skip, false);
                    }
                }
                let len = w.g - g;
                stored += len - k as u64 + 1;
                store.set_view(i as PieceId, g, len as u32);
            }
            let end = w.finish();
            store
                .frontier
                .store(end.div_ceil(CHUNK_BASES) * CHUNK_BASES, Ordering::Relaxed);
            store.n_pids.store(chains.len() as u64, Ordering::Relaxed);
        }
        let n_new = chains.len();
        drop(chains);
        self.store = store;

        // 4. Re-index every piece into shards sized for their share at
        //    ≈ 80 % head load (as KmerTable::defrag); hot minimizers stay hot.
        let size = |e: u64, shards: usize| {
            ((e as f64 / shards as f64 / (0.8 * crate::index::BLOCK_SLOTS as f64)).ceil() as usize)
                .max(1)
        };
        let (b1, b2) = (size(n1, self.l1.len()), size(n2, self.l2s.len()));
        let hot1: Vec<u64> = rlock(&self.hot1).iter().copied().collect();
        let hot2: Vec<u64> = rlock(&self.hot2).iter().copied().collect();
        for (shards, b, hot) in [(&mut self.l1, b1, &hot1), (&mut self.l2s, b2, &hot2)] {
            for s in shards.iter_mut() {
                *s.get_mut().unwrap_or_else(PoisonError::into_inner) =
                    MinimizerIndex::with_buckets(b).with_spill_dir(spill.clone());
            }
            let ns = shards.len();
            for &m in hot {
                shards[shard_of(m, ns)]
                    .get_mut()
                    .unwrap_or_else(PoisonError::into_inner)
                    .mark_hot(m);
            }
        }
        self.reindex_pieces(n_new);
        for s in self.l1.iter_mut().chain(self.l2s.iter_mut()) {
            s.get_mut()
                .unwrap_or_else(PoisonError::into_inner)
                .mark_baseline();
        }

        // 5. Counters, as KmerTable::defrag.
        let c = self.carry.get_mut().unwrap_or_else(PoisonError::into_inner);
        c.stored_kmers = stored;
        c.pieces_created = 0;
        c.dirty = Vec::new();
        c.changed = Vec::new();
        c.pids_at_start = n_new as u64;
        c.pieces_at_defrag = n_new as u64;
        let e = &mut c.edits;
        e.defrags += 1;
        e.joins += rep.joins;
        e.dead_bases = 0;
        e.shared_bases = 0;
        e.cuts_since_defrag = 0;
        e.inserted_pieces_since_defrag = 0;
        self.report_after(rep)
    }

    /// File the index items of pieces `0..n` (every window: layer 1, or
    /// layer 2 / 3 if its minimizer is hot), from all cores.
    fn reindex_pieces(&self, n: usize) {
        let k = self.cfg.k;
        let hot1 = rlock(&self.hot1);
        let hot2 = rlock(&self.hot2);
        let (hot1, hot2): (&HashSet<u64>, &HashSet<u64>) = (&hot1, &hot2);
        let bases: u64 = (0..n).map(|p| self.store.view(p as PieceId).1 as u64).sum();
        let threads = conversion_threads(bases as usize / 8);
        par_ranges(n, 1024, threads, |a, b| {
            let mut items = Vec::new();
            let mut d = Del::default();
            for pid in a..b {
                let (g0, len) = self.store.view(pid as PieceId);
                let len = len as usize;
                if len < k {
                    continue;
                }
                items.clear();
                self.window_items(g0, 0, len - k, hot1, hot2, &mut items);
                items.sort_unstable();
                items.dedup();
                self.apply(&items, pid as PieceId, true, &mut d);
            }
        });
    }

    /// Compact every index shard, and re-size those whose heads would be at
    /// least [`INDEX_SHRINK_SLACK`] fewer at 80 % load.
    fn shrink_shards(&self) {
        let (l, l2) = (self.cfg.l, self.l2);
        let st = &self.store;
        for (shards, ll) in [(&self.l1, l), (&self.l2s, l2)] {
            let threads = conversion_threads(shards.len() * 256);
            par_ranges(shards.len(), 16, threads, |a, b| {
                for s in &shards[a..b] {
                    let mut ix = wlock(s);
                    if ix.orphans() > 0 {
                        ix.compact();
                    }
                    let slots = 0.8 * crate::index::BLOCK_SLOTS as f64;
                    let target = (ix.num_entries() as f64 / slots).ceil().max(1.0);
                    if target < ix.num_buckets() as f64 * (1.0 - INDEX_SHRINK_SLACK) {
                        ix.shrink_to_load(0.8, |p, q| {
                            order_value(codec::canonicalize(st.kmer_at(p, q, ll), ll))
                        });
                    }
                }
            });
        }
    }

    /// [`KmerTable::defrag_incremental`] (without merging) on the
    /// concurrent table, in place: rewrite the base array without dead
    /// bases (piece ids and index entries unchanged), drop unreachable
    /// overflow blocks, and re-size to 80 % head load the shards where that
    /// saves at least [`INDEX_SHRINK_SLACK`] of their heads. O(live bases +
    /// entries of the re-sized shards), all cores for the shards.
    pub fn defrag_incremental(&mut self) -> DefragReport {
        let rep = self.report_before();
        // Live views by start; overlapping views (cut parts share k − 1
        // bases) are copied once.
        let n = self.store.n_pids.load(Ordering::Relaxed);
        let mut order: Vec<(u64, u32, PieceId)> = (0..n as PieceId)
            .filter_map(|p| {
                let (s, l) = self.store.view(p);
                (l > 0).then_some((s, l, p))
            })
            .collect();
        order.sort_unstable();
        let store = SharedStore::new();
        let mut shared = 0u64;
        {
            let old = &self.store;
            let mut w = SeqWriter::new(&store);
            let (mut seg_lo, mut seg_hi, mut seg_new) = (0u64, 0u64, 0u64);
            let mut have = false;
            for &(s, len, pid) in &order {
                let e = s + len as u64;
                let ns = if have && s < seg_hi {
                    shared += seg_hi.min(e) - s;
                    if e > seg_hi {
                        w.copy(old, seg_hi, e - seg_hi, false);
                        seg_hi = e;
                    }
                    seg_new + (s - seg_lo)
                } else {
                    (seg_lo, seg_hi, seg_new, have) = (s, e, w.g, true);
                    w.copy(old, s, len as u64, false);
                    seg_new
                };
                store.set_view(pid, ns, len);
            }
            let end = w.finish();
            store
                .frontier
                .store(end.div_ceil(CHUNK_BASES) * CHUNK_BASES, Ordering::Relaxed);
        }
        drop(order);
        // Ids: all kept (emptied ones stay free for reuse).
        let old = std::mem::replace(&mut self.store, store);
        self.store.n_pids.store(n, Ordering::Relaxed);
        let lock_vec =
            |m: Mutex<Vec<PieceId>>| m.into_inner().unwrap_or_else(PoisonError::into_inner);
        let free = lock_vec(old.free);
        self.store
            .free_left
            .store(!free.is_empty(), Ordering::Relaxed);
        *lock(&self.store.free) = free;
        *lock(&self.store.reused) = lock_vec(old.reused);
        // Unused chunk tails pointed into the old array; id blocks stay valid.
        *lock(&self.store.pid_tails) = old
            .pid_tails
            .into_inner()
            .unwrap_or_else(PoisonError::into_inner);
        drop(old.data);

        // Index: drop orphans, and re-size (to ≈ 80 % head load) every
        // shard where that saves at least INDEX_SHRINK_SLACK of its blocks
        // (head slots freed by deletions, or a shard sized ahead by
        // with_capacity). All cores, one shard at a time each.
        self.shrink_shards();
        for s in self.l3.iter_mut() {
            s.get_mut()
                .unwrap_or_else(PoisonError::into_inner)
                .shrink_to_fit();
        }

        let c = self.carry.get_mut().unwrap_or_else(PoisonError::into_inner);
        let e = &mut c.edits;
        e.incremental_defrags += 1;
        e.dead_bases = 0;
        e.shared_bases = shared;
        self.report_after(rep)
    }
}
