//! A read-only view of a [`ConcurrentKmerTable`] whose lookups take no
//! locks.
//!
//! [`ConcurrentKmerTable::read_snapshot`] acquires the read lock of every
//! index shard once, in the table's lock order (layer 1, layer 2, layer 3),
//! and keeps them until the snapshot is dropped. Queries through the
//! snapshot ([`ReadSnapshot::contains`], [`ReadSnapshot::query_stream`]) then
//! read the shards without locking, from any number of threads (the snapshot
//! is `Sync`). Writers and deleters block until the snapshot is dropped, so
//! it suits read-only phases: the answers are those of the table at the
//! moment the snapshot was taken.
//!
//! The snapshot runs exactly the same lookup code on the same data as the
//! locked queries, without the per-query lock acquisitions. Comparing the two
//! measures what the shard read locks cost, an upper bound on what any
//! lock-free read path could gain on this structure.

use super::*;

/// Every index shard of a [`ConcurrentKmerTable`], read-locked once (see the
/// [module documentation](self)).
pub struct ReadSnapshot<'a> {
    pub(super) t: &'a ConcurrentKmerTable,
    pub(super) l1: Vec<RwLockReadGuard<'a, MinimizerIndex>>,
    l2s: Vec<RwLockReadGuard<'a, MinimizerIndex>>,
    l3: Vec<RwLockReadGuard<'a, Layer3>>,
}

impl ConcurrentKmerTable {
    /// Read-lock every index shard once and return a view whose queries
    /// take no locks. Writers and deleters block until it is dropped.
    pub fn read_snapshot(&self) -> ReadSnapshot<'_> {
        let l1 = self.l1.iter().map(|s| rlock(s)).collect();
        let l2s = self.l2s.iter().map(|s| rlock(s)).collect();
        let l3 = self.l3.iter().map(|s| rlock(s)).collect();
        ReadSnapshot {
            t: self,
            l1,
            l2s,
            l3,
        }
    }
}

impl<'a> ReadSnapshot<'a> {
    /// The table this snapshot reads.
    pub fn table(&self) -> &'a ConcurrentKmerTable {
        self.t
    }

    /// Exact membership of a packed k-mer (either orientation), as
    /// [`ConcurrentKmerTable::contains`], without locking.
    pub fn contains(&self, kmer: u64) -> bool {
        self.contains_traced(kmer, &mut QueryTrace::default())
    }

    /// [`contains`](Self::contains) that also reports the work done.
    pub fn contains_traced(&self, kmer: u64, trace: &mut QueryTrace) -> bool {
        let t = self.t;
        if kmer & !t.mask != 0 {
            return false;
        }
        let x = codec::canonicalize(kmer, t.cfg.k);
        let (m, ties) = kmer_minimizer_mask(x, t.cfg.k, t.cfg.l);
        let g = &self.l1[shard_of(m, self.l1.len())];
        if !(t.l2 > 0 && g.is_hot(m)) {
            return t
                .find(&t.store, g, x, m, ties, t.cfg.l, t.w, trace)
                .is_some();
        }
        self.locate_hot(x, trace).is_some()
    }

    /// [`ConcurrentKmerTable::locate_hot`] without locking.
    pub(super) fn locate_hot(&self, x: u64, trace: &mut QueryTrace) -> Option<Hit> {
        let t = self.t;
        let k = t.cfg.k;
        let (m2, ties2) = kmer_minimizer_mask(x, k, t.l2);
        trace.layer2 += 1;
        let g = &self.l2s[shard_of(m2, self.l2s.len())];
        if !g.is_hot(m2) {
            return t.find(&t.store, g, x, m2, ties2, t.l2, t.w2, trace);
        }
        trace.layer2 -= 1;
        trace.layer3 += 1;
        let loc = self.l3[shard_of(x, LAYER3_SHARDS)].get(x);
        loc.map(|(piece, pos)| Hit {
            piece,
            pos,
            kmer: t.store.kmer_at(piece, pos, k),
        })
    }
}
