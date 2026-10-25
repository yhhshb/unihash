//! In-order ("streaming") queries over the successive k-mers of a sequence.
//!
//! Successive k-mers of a read overlap by `k - 1` bases, so most of the work
//! of one lookup can be reused for the next:
//!
//! 1. **Rolling minimizer.** The window minimizer `m` and its tie mask come
//!    from a monotone deque updated with one ℓ-mer per base (O(1)
//!    amortised) instead of being recomputed from the k ℓ-mers (O(k)).
//! 2. **Extension of the previous hit.** If the previous k-mer was found at
//!    `(piece, s)`, the next one is at `s + 1` (same orientation as the
//!    read) or `s - 1` (reverse complement) exactly when one stored base
//!    equals the incoming base (or its complement). That single comparison
//!    is exact, needs no index access, and does not depend on `m`.
//! 3. **Candidate cache while the minimizer does not change.** The
//!    tag-matching entries of `m`'s bucket are collected once and kept while
//!    successive k-mers share the same `m`, so the bucket is scanned once
//!    per minimizer run rather than once per k-mer. Only the exact
//!    candidate comparisons (with the current k-mer's tie offsets) remain.
//!
//! Answers are identical to [`KmerTable::contains`] for every k-mer.

use crate::codec::{canonical, encode, RollingKmer, INVALID};
use crate::index::ScanTrace;
use crate::minimizer::{order_value, WindowMin};
use crate::piece_store::PieceId;
use crate::table::{Hit, KmerTable, QueryTrace};
use std::fmt;

/// Work counters of a [`StreamingQuery`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct StreamStats {
    /// k-mers answered.
    pub queries: u64,
    pub hits: u64,
    /// Hits answered by extending the previous hit (one base comparison).
    pub extensions: u64,
    /// Index lookups whose minimizer equalled the cached one (no bucket scan).
    pub cache_reuses: u64,
    /// Index lookups that scanned a bucket (minimizer changed).
    pub bucket_scans: u64,
    /// Lookups routed to the fallback layers (hot ℓ-minimizer).
    pub fallback_lookups: u64,
    /// Entries visited / tag-matched during bucket scans.
    pub entries_scanned: u64,
    pub tag_matches: u64,
    /// Exact candidate k-mers compared (index path only).
    pub candidate_comparisons: u64,
}

impl StreamStats {
    fn per(&self, x: u64) -> f64 {
        x as f64 / self.queries.max(1) as f64
    }

    pub fn add(&mut self, o: &StreamStats) {
        self.queries += o.queries;
        self.hits += o.hits;
        self.extensions += o.extensions;
        self.cache_reuses += o.cache_reuses;
        self.bucket_scans += o.bucket_scans;
        self.fallback_lookups += o.fallback_lookups;
        self.entries_scanned += o.entries_scanned;
        self.tag_matches += o.tag_matches;
        self.candidate_comparisons += o.candidate_comparisons;
    }
}

impl fmt::Display for StreamStats {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "  queries                      {}", self.queries)?;
        writeln!(f, "  hits                         {}", self.hits)?;
        writeln!(
            f,
            "  answered by extension        {:.4}",
            self.per(self.extensions)
        )?;
        writeln!(
            f,
            "  index lookups, cached m      {:.4}",
            self.per(self.cache_reuses)
        )?;
        writeln!(
            f,
            "  index lookups, bucket scan   {:.4}",
            self.per(self.bucket_scans)
        )?;
        writeln!(
            f,
            "  index lookups, fallback      {:.4}",
            self.per(self.fallback_lookups)
        )?;
        writeln!(
            f,
            "  entries scanned / query      {:.3}",
            self.per(self.entries_scanned)
        )?;
        write!(
            f,
            "  candidate k-mers / query     {:.3}",
            self.per(self.candidate_comparisons)
        )
    }
}

/// Where the previous k-mer sits in the packed piece array.
#[derive(Clone, Copy, Debug)]
struct Cursor {
    /// Global base range of the piece.
    lo: u64,
    hi: u64,
    /// Global offset of the k-mer's first base.
    g: u64,
    /// The piece holds the k-mer in the read's orientation (`true`) or
    /// reverse-complemented (`false`).
    same: bool,
}

/// Answers the k-mers of one or more sequences, in order. Feed sequence with
/// [`push`](Self::push) (chunks may be cut anywhere) and call
/// [`end`](Self::end) between sequences. k-mers overlapping a non-ACGT
/// character are skipped, as in insertion.
#[derive(Debug)]
pub struct StreamingQuery<'a> {
    table: &'a KmerTable,
    k: usize,
    roll: RollingKmer,
    lroll: RollingKmer,
    wmin: WindowMin,
    /// Location of the previous k-mer (if it was found).
    prev: Option<Cursor>,
    /// Tag-matching candidates of `cache_m` (valid while `cache_valid`).
    cache_m: u64,
    cache_valid: bool,
    cache: Vec<(PieceId, u32)>,
    stats: StreamStats,
}

impl<'a> StreamingQuery<'a> {
    pub fn new(table: &'a KmerTable) -> Self {
        let k = table.config().k;
        let l = table.config().l;
        Self {
            table,
            k,
            roll: RollingKmer::new(k),
            lroll: RollingKmer::new(l),
            wmin: WindowMin::new(k - l + 1),
            prev: None,
            cache_m: 0,
            cache_valid: false,
            cache: Vec::new(),
            stats: StreamStats::default(),
        }
    }

    pub fn stats(&self) -> &StreamStats {
        &self.stats
    }

    /// End the current sequence: the next k-mer starts a new window.
    /// (The candidate cache is kept: it depends only on the minimizer.)
    pub fn end(&mut self) {
        self.roll.reset();
        self.lroll.reset();
        self.wmin.clear();
        self.prev = None;
    }

    /// Feed the next chunk. `f(x, present)` is called for every complete
    /// k-mer in order, with `x` canonical.
    pub fn push<F: FnMut(u64, bool)>(&mut self, chunk: &[u8], mut f: F) {
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
            f(canonical(self.roll.fwd(), self.roll.rc()), present);
        }
    }

    /// Answer the k-mer that ends with base `b` (already pushed).
    #[inline]
    fn step(&mut self, b: u8) -> bool {
        let k = self.k;
        let fwd = self.roll.fwd();
        let rc = self.roll.rc();
        self.stats.queries += 1;

        // 1. Extend the previous hit by one base (one 2-bit load, no index
        //    access and no piece-table lookup).
        if let Some(c) = self.prev.as_mut() {
            let store = self.table.pieces();
            let k = k as u64;
            let ok = if c.same {
                c.g + k < c.hi && store.code_at_global(c.g + k) == b
            } else {
                c.g > c.lo && store.code_at_global(c.g - 1) == 3 - b
            };
            if ok {
                if c.same {
                    c.g += 1;
                } else {
                    c.g -= 1;
                }
                self.stats.extensions += 1;
                self.stats.hits += 1;
                return true;
            }
        }

        // 2. Index lookup, reusing the candidates while m is unchanged.
        let (m, ties) = self.wmin.current().expect("window complete");
        let x = canonical(fwd, rc);
        if self.table.is_hot(m) {
            // Hot ℓ-minimizer: routed to layers 2/3 (no candidate cache).
            let mut tr = QueryTrace::default();
            let found = self.table.find_fallback(x, &mut tr);
            self.stats.fallback_lookups += 1;
            self.stats.entries_scanned += tr.entries_scanned;
            self.stats.tag_matches += tr.tag_matches;
            self.stats.candidate_comparisons += tr.candidate_comparisons;
            return self.finish(found, fwd);
        }
        if self.cache_valid && self.cache_m == m {
            self.stats.cache_reuses += 1;
        } else {
            self.cache.clear();
            let mut scan = ScanTrace::default();
            self.table
                .index()
                .collect_candidates(m, &mut scan, &mut self.cache);
            self.stats.bucket_scans += 1;
            self.stats.entries_scanned += scan.entries_scanned;
            self.stats.tag_matches += scan.tag_matches;
            self.cache_m = m;
            self.cache_valid = true;
        }
        let xr = fwd ^ rc ^ x; // the other orientation
        let ties = self.table.mirrored_ties(ties);
        let l = self.table.config().l;
        let mut comparisons = 0u64;
        let mut found = None;
        for &(pid, q) in &self.cache {
            found = self
                .table
                .check_candidate(x, xr, ties, pid, q, l, &mut comparisons);
            if found.is_some() {
                break;
            }
        }
        self.stats.candidate_comparisons += comparisons;
        self.finish(found, fwd)
    }

    /// Record the answer; a hit in a piece becomes the extension cursor.
    #[inline]
    fn finish(&mut self, found: Option<Hit>, fwd: u64) -> bool {
        let store = self.table.pieces();
        self.prev = match found {
            Some(h) => {
                let (lo, len) = store.view(h.piece);
                Some(Cursor {
                    lo,
                    hi: lo + len as u64,
                    g: lo + h.pos as u64,
                    same: h.kmer == fwd,
                })
            }
            None => None,
        };
        self.stats.hits += found.is_some() as u64;
        found.is_some()
    }
}
