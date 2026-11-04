//! In-order ("streaming") queries on a [`ConcurrentKmerTable`], safe while
//! other threads insert and delete.
//!
//! The same shortcuts as [`crate::StreamingQuery`], each kept exact under
//! concurrent changes:
//!
//! * **Extension of the previous hit** by one base: the piece's current
//!   view must still cover the next position (deletions only shrink views;
//!   bases never change), so a hit is a k-mer present when the view is read.
//! * **Candidate cache while the minimizer does not change**: the query
//!   keeps the read lock of the minimizer's layer-1 shard while successive
//!   k-mers share it, and the candidates collected under that lock are used
//!   only while it is held, so nobody can have changed the shard since.
//!
//! Every answer is therefore what [`ConcurrentKmerTable::contains`] would
//! return at some moment during that k-mer's lookup. A lock is held only
//! inside [`push`](ConcurrentStreamingQuery::push); the answers of a chunk
//! are handed to the callback after it is released, so the callback may
//! insert into or delete from the same table.
//!
//! Created from a [`ReadSnapshot`] ([`ReadSnapshot::query_stream`]), the
//! query runs the same steps but reads the snapshot's shards without
//! locking.

use super::snapshot::ReadSnapshot;
use super::*;
use crate::stream_query::StreamStats;

/// Answers the k-mers of one or more sequences in order, on a
/// [`ConcurrentKmerTable`] (see the [module documentation](self)). Create
/// one per thread with [`ConcurrentKmerTable::query_stream`].
pub struct ConcurrentStreamingQuery<'a> {
    t: &'a ConcurrentKmerTable,
    /// Read snapshot whose shards are read without locking, if any.
    snap: Option<&'a ReadSnapshot<'a>>,
    roll: RollingKmer,
    lroll: RollingKmer,
    wmin: WindowMin,
    prev: Option<Cursor>,
    /// Read lock of a layer-1 shard, and whether `cache` holds the
    /// candidates of `cache_m` collected under it.
    held: Option<(usize, RwLockReadGuard<'a, MinimizerIndex>)>,
    cache_m: u64,
    cache_valid: bool,
    cache: Vec<(PieceId, u32)>,
    /// Minimizers seen hot (hotness is permanent).
    known_hot: FastSet,
    /// Answers of the current chunk, handed out once no lock is held.
    out: Vec<(u64, bool)>,
    stats: StreamStats,
}

impl ConcurrentKmerTable {
    /// Answer the k-mers of sequences in order, reusing work between
    /// successive k-mers (see [`ConcurrentStreamingQuery`]).
    pub fn query_stream(&self) -> ConcurrentStreamingQuery<'_> {
        ConcurrentStreamingQuery::new(self, None)
    }
}

impl<'a> ReadSnapshot<'a> {
    /// Answer the k-mers of sequences in order, as
    /// [`ConcurrentKmerTable::query_stream`], without locking.
    pub fn query_stream(&self) -> ConcurrentStreamingQuery<'_> {
        ConcurrentStreamingQuery::new(self.t, Some(self))
    }
}

impl<'a> ConcurrentStreamingQuery<'a> {
    fn new(t: &'a ConcurrentKmerTable, snap: Option<&'a ReadSnapshot<'a>>) -> Self {
        let (k, l) = (t.cfg.k, t.cfg.l);
        ConcurrentStreamingQuery {
            t,
            snap,
            roll: RollingKmer::new(k),
            lroll: RollingKmer::new(l),
            wmin: WindowMin::new(k - l + 1),
            prev: None,
            held: None,
            cache_m: 0,
            cache_valid: false,
            cache: Vec::new(),
            known_hot: FastSet::default(),
            out: Vec::new(),
            stats: StreamStats::default(),
        }
    }
}

impl<'a> ConcurrentStreamingQuery<'a> {
    pub fn stats(&self) -> &StreamStats {
        &self.stats
    }

    /// End the current sequence: the next k-mer starts a new window.
    pub fn end(&mut self) {
        self.roll.reset();
        self.lroll.reset();
        self.wmin.clear();
        self.prev = None;
    }

    /// Feed the next chunk. `f(x, present)` is called for every complete
    /// k-mer in order (`x` canonical), after the whole chunk was looked up
    /// and every lock released. k-mers overlapping a non-ACGT character are
    /// skipped.
    pub fn push<F: FnMut(u64, bool)>(&mut self, chunk: &[u8], mut f: F) {
        self.out.clear();
        for &c in chunk {
            let b = encode(c);
            if b == INVALID {
                self.end();
                continue;
            }
            self.roll.push(b);
            if self.lroll.push(b) {
                self.wmin.push(order_value(self.lroll.canonical()));
            }
            if !self.roll.is_full() {
                continue;
            }
            let present = self.step(b);
            self.out.push((self.roll.canonical(), present));
        }
        self.release();
        for &(x, p) in &self.out {
            f(x, p);
        }
    }

    /// Drop the shard lock (and with it the candidate cache).
    fn release(&mut self) {
        self.held = None;
        self.cache_valid = false;
    }

    /// Answer the k-mer that ends with base `b` (already pushed).
    fn step(&mut self, b: u8) -> bool {
        let t = self.t;
        let fwd = self.roll.fwd();
        self.stats.queries += 1;
        if let Some(c) = self.prev {
            if let Some(c) = t.extend(&t.store, c, b) {
                self.prev = Some(c);
                self.stats.extensions += 1;
                self.stats.hits += 1;
                return true;
            }
        }
        let (m, ties) = self.wmin.current().expect("window complete");
        let x = self.roll.canonical();
        let mut trace = QueryTrace::default();
        let mut hot = t.l2 > 0 && self.known_hot.contains(&m);
        let found = if !hot {
            let sh = shard_of(m, t.l1.len());
            let g: &MinimizerIndex = match self.snap {
                Some(s) => &s.l1[sh],
                None => {
                    if !matches!(&self.held, Some((h, _)) if *h == sh) {
                        self.release();
                        self.held = Some((sh, rlock(&t.l1[sh])));
                    }
                    &self.held.as_ref().expect("lock held").1
                }
            };
            if t.l2 > 0 && g.is_hot(m) {
                hot = true;
                None
            } else {
                if self.cache_valid && self.cache_m == m {
                    self.stats.cache_reuses += 1;
                } else {
                    self.cache.clear();
                    let mut scan = ScanTrace::default();
                    g.collect_candidates(m, &mut scan, &mut self.cache);
                    self.stats.bucket_scans += 1;
                    self.stats.entries_scanned += scan.entries_scanned;
                    self.stats.tag_matches += scan.tag_matches;
                    self.cache_m = m;
                    self.cache_valid = true;
                }
                t.find_in(&t.store, &self.cache, x, ties, t.cfg.l, t.w, &mut trace)
            }
        } else {
            None
        };
        let found = if hot {
            self.known_hot.insert(m);
            self.release();
            self.stats.fallback_lookups += 1;
            let h = match self.snap {
                Some(s) => s.locate_hot(x, &mut trace),
                None => t.locate_hot(x, &mut trace),
            };
            self.stats.entries_scanned += trace.entries_scanned;
            self.stats.tag_matches += trace.tag_matches;
            h
        } else {
            found
        };
        self.stats.candidate_comparisons += trace.candidate_comparisons;
        self.stats.hits += found.is_some() as u64;
        self.prev = found.and_then(|h| t.cursor(&t.store, &h, fwd));
        found.is_some()
    }
}
