//! The public k-mer presence table.

use crate::codec::{self, encode, kmer_mask, RollingKmer, INVALID, MAX_K};
use crate::index::{MinimizerIndex, ScanTrace, SpillStats};
use crate::minimizer::{
    kmer_minimizer_mask, lmer_values_into, order_value, piece_minimizer_positions_into, WindowMin,
};
use crate::piece_store::{PieceId, PieceStore, MAX_PIECE_LEN};
use crate::stream_query::StreamingQuery;
use std::collections::HashMap;
use std::fmt;
use std::path::PathBuf;

/// Construction parameters. See [`Config::validate`] for the allowed ranges.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Config {
    /// k-mer length, `2 ≤ k ≤ 32`.
    pub k: usize,
    /// Minimizer length ℓ, `1 ≤ ℓ ≤ k`.
    pub l: usize,
    /// The index starts with `2^bucket_bits` buckets (`0 ≤ bucket_bits ≤ 30`)
    /// and grows (×1.5) whenever its fixed-capacity bucket heads would be
    /// more than 85 % full, so it does not need to be sized in advance.
    pub bucket_bits: u32,
    /// Longest piece in bases (`k ≤ max_piece_len ≤ 2^24 - 1`, default
    /// 4096). Longer runs of absent k-mers are split into pieces overlapping
    /// by `k - 1` bases, and defragmentation does not merge pieces beyond
    /// this length. Deleting from a piece re-files the index entries of its
    /// right part, so this length bounds the cost of a deletion; a split
    /// costs about 15 B (piece record and overlap), i.e. < 0.1 % of the
    /// table at 4096.
    pub max_piece_len: u32,
    /// Which piece positions a tag-matching index entry is checked at.
    pub scan: CandidateScan,
    /// Minimizer length of the fallback index for hot minimizers
    /// (`ℓ < ℓ2 ≤ k`; `0` = automatic, `min(k, ℓ + 8)`).
    pub l2: usize,
    /// A minimizer with at least this many index entries is *hot*: its
    /// k-mers are re-indexed with ℓ2-minimizers (layer 2), and those whose
    /// ℓ2-minimizer is hot too go to a set of full k-mers (layer 3).
    /// `0` disables routing.
    pub hot_threshold: u32,
    /// Directory for temporary spill files written while the index grows
    /// (`None`: spill to a compact in-memory list).
    pub spill_dir: Option<PathBuf>,
    /// When the k-mers of a run of absent k-mers become visible to lookups.
    pub run_indexing: RunIndexing,
}

/// When insertion indexes the run of absent k-mers it is building.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum RunIndexing {
    /// The run is stored as a piece and indexed when it ends. A k-mer that
    /// occurs twice inside one run is stored twice (see
    /// [`KmerTable::redundant_kmers`]).
    Deferred,
    /// The run's piece is created at its first k-mer and every k-mer is
    /// indexed as soon as it is appended, so it is visible to the next
    /// lookup: no k-mer is ever stored twice, and an open read's k-mers are
    /// already queryable. Required by [`crate::ConcurrentKmerTable`].
    #[default]
    Eager,
}

impl std::str::FromStr for RunIndexing {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, String> {
        match s.to_ascii_lowercase().as_str() {
            "deferred" => Ok(Self::Deferred),
            "eager" => Ok(Self::Eager),
            _ => Err(format!(
                "unknown run indexing '{s}' (expected deferred|eager)"
            )),
        }
    }
}

/// Candidate positions checked for an index entry `(pid, q)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum CandidateScan {
    /// Every start `s` in `[q - (k-ℓ), q]` clipped to the piece, exactly as
    /// (up to `k - ℓ + 1` comparisons per entry).
    Range,
    /// Only `s = q - p` for the query's tied positions `p ∈ Q(x)` and their
    /// mirror images `(k-ℓ) - p` (the piece may hold the reverse
    /// complement). Exact as well: if `x` occurs at `s` in either
    /// orientation, `q - s` is one of those offsets. Usually one or two
    /// comparisons per entry instead of `k - ℓ + 1`.
    #[default]
    TieOffsets,
}

impl std::str::FromStr for CandidateScan {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, String> {
        match s.to_ascii_lowercase().as_str() {
            "range" => Ok(Self::Range),
            "ties" | "tie-offsets" => Ok(Self::TieOffsets),
            _ => Err(format!(
                "unknown candidate scan '{s}' (expected range|ties)"
            )),
        }
    }
}

impl Default for Config {
    fn default() -> Self {
        Self {
            k: 31,
            l: 17,
            bucket_bits: 10,
            max_piece_len: 4096,
            scan: CandidateScan::TieOffsets,
            l2: 0,
            hot_threshold: 64,
            spill_dir: None,
            run_indexing: RunIndexing::default(),
        }
    }
}

/// Construction-time misconfiguration.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ConfigError {
    KOutOfRange(usize),
    LOutOfRange { l: usize, k: usize },
    BucketBitsTooLarge(u32),
    MaxPieceLenTooSmall { max_piece_len: u32, k: usize },
    MaxPieceLenTooLarge(u32),
    L2OutOfRange { l2: usize, l: usize, k: usize },
}

impl fmt::Display for ConfigError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::KOutOfRange(k) => write!(
                f,
                "k = {k} is out of range: k must satisfy 2 ≤ k ≤ {MAX_K} so that a k-mer fits in 64 bits"
            ),
            Self::LOutOfRange { l, k } => write!(
                f,
                "minimizer length ℓ = {l} is out of range: ℓ must satisfy 1 ≤ ℓ ≤ k (k = {k})"
            ),
            Self::BucketBitsTooLarge(b) => write!(
                f,
                "bucket_bits = {b} is too large: at most 30 (2^30 buckets) is supported"
            ),
            Self::MaxPieceLenTooSmall { max_piece_len, k } => write!(
                f,
                "max_piece_len = {max_piece_len} is smaller than k = {k}: a piece must hold at least one k-mer"
            ),
            Self::MaxPieceLenTooLarge(m) => write!(
                f,
                "max_piece_len = {m} is too large: pieces are at most {MAX_PIECE_LEN} bases (24-bit length)"
            ),
            Self::L2OutOfRange { l2, l, k } => write!(
                f,
                "fallback minimizer length ℓ2 = {l2} is out of range: ℓ2 must satisfy ℓ < ℓ2 ≤ k (ℓ = {l}, k = {k}), or be 0 for automatic"
            ),
        }
    }
}

impl std::error::Error for ConfigError {}

impl Config {
    pub fn validate(&self) -> Result<(), ConfigError> {
        if !(2..=MAX_K).contains(&self.k) {
            return Err(ConfigError::KOutOfRange(self.k));
        }
        if !(1..=self.k).contains(&self.l) {
            return Err(ConfigError::LOutOfRange {
                l: self.l,
                k: self.k,
            });
        }
        if self.bucket_bits > 30 {
            return Err(ConfigError::BucketBitsTooLarge(self.bucket_bits));
        }
        if (self.max_piece_len as usize) < self.k {
            return Err(ConfigError::MaxPieceLenTooSmall {
                max_piece_len: self.max_piece_len,
                k: self.k,
            });
        }
        if self.max_piece_len > MAX_PIECE_LEN {
            return Err(ConfigError::MaxPieceLenTooLarge(self.max_piece_len));
        }
        if self.l2 != 0 && !(self.l + 1..=self.k).contains(&self.l2) {
            return Err(ConfigError::L2OutOfRange {
                l2: self.l2,
                l: self.l,
                k: self.k,
            });
        }
        Ok(())
    }

    /// The fallback minimizer length in use, or `None` if hot-minimizer
    /// routing is disabled (threshold 0, or no ℓ2 > ℓ fits in k).
    pub fn fallback_l(&self) -> Option<usize> {
        if self.hot_threshold == 0 {
            return None;
        }
        let l2 = if self.l2 == 0 {
            (self.l + 8).min(self.k)
        } else {
            self.l2
        };
        (l2 > self.l).then_some(l2)
    }
}

/// Error for string-based queries.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum KmerError {
    WrongLength { expected: usize, got: usize },
    InvalidBase { pos: usize, byte: u8 },
}

impl fmt::Display for KmerError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::WrongLength { expected, got } => {
                write!(f, "k-mer has length {got}, expected {expected}")
            }
            Self::InvalidBase { pos, byte } => {
                write!(
                    f,
                    "non-ACGT character {:?} at position {pos}",
                    *byte as char
                )
            }
        }
    }
}

impl std::error::Error for KmerError {}

/// Per-query work counters.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct QueryTrace {
    /// Index entries visited in the minimizer's bucket.
    pub entries_scanned: u64,
    /// Entries whose 8-bit tag matched.
    pub tag_matches: u64,
    /// Candidate k-mers read from pieces and compared exactly.
    pub candidate_comparisons: u64,
    /// Lookups routed to the fallback ℓ2 index (hot ℓ-minimizer).
    pub layer2: u64,
    /// Lookups routed to the full k-mer set (hot ℓ- and ℓ2-minimizers).
    pub layer3: u64,
}

impl QueryTrace {
    pub(crate) fn add_scan(&mut self, s: &ScanTrace) {
        self.entries_scanned += s.entries_scanned;
        self.tag_matches += s.tag_matches;
    }
}

/// Location of a k-mer found in the table.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Hit {
    pub piece: PieceId,
    /// Start of the k-mer in the piece.
    pub pos: u32,
    /// The k-mer as stored in the piece (forward orientation of the piece);
    /// equals the query or its reverse complement.
    pub kmer: u64,
}

/// Layer 3: full canonical k-mers (whose ℓ- and ℓ2-minimizers are hot)
/// with the piece location of each stored occurrence, so that deletion can
/// cut the piece. A k-mer repeated inside one piece has several locations;
/// the extra ones (rare) live in `more`.
#[derive(Clone, Debug, Default)]
pub(crate) struct Layer3 {
    map: HashMap<u64, (PieceId, u32)>,
    more: HashMap<u64, Vec<(PieceId, u32)>>,
}

impl Layer3 {
    #[inline]
    pub(crate) fn get(&self, x: u64) -> Option<(PieceId, u32)> {
        self.map.get(&x).copied()
    }

    pub(crate) fn len(&self) -> usize {
        self.map.len()
    }

    /// Record an occurrence (idempotent).
    pub(crate) fn add(&mut self, x: u64, loc: (PieceId, u32)) {
        match self.map.entry(x) {
            std::collections::hash_map::Entry::Vacant(e) => {
                e.insert(loc);
            }
            std::collections::hash_map::Entry::Occupied(e) => {
                if *e.get() != loc {
                    let v = self.more.entry(x).or_default();
                    if !v.contains(&loc) {
                        v.push(loc);
                    }
                }
            }
        }
    }

    /// Forget one occurrence; returns whether it was recorded.
    pub(crate) fn remove(&mut self, x: u64, loc: (PieceId, u32)) -> bool {
        let Some(&cur) = self.map.get(&x) else {
            return false;
        };
        if cur == loc {
            match self.more.get_mut(&x) {
                Some(v) => {
                    let next = v.pop().expect("non-empty");
                    if v.is_empty() {
                        self.more.remove(&x);
                    }
                    self.map.insert(x, next);
                }
                None => {
                    self.map.remove(&x);
                }
            }
            return true;
        }
        if let Some(v) = self.more.get_mut(&x) {
            if let Some(i) = v.iter().position(|&l| l == loc) {
                v.swap_remove(i);
                if v.is_empty() {
                    self.more.remove(&x);
                }
                return true;
            }
        }
        false
    }

    /// Every recorded `(k-mer, location)`.
    pub(crate) fn iter(&self) -> impl Iterator<Item = (u64, (PieceId, u32))> + '_ {
        self.map.iter().map(|(&x, &l)| (x, l)).chain(
            self.more
                .iter()
                .flat_map(|(&x, v)| v.iter().map(move |&l| (x, l))),
        )
    }

    pub(crate) fn clear(&mut self) {
        self.map = HashMap::new();
        self.more = HashMap::new();
    }

    pub(crate) fn heap_bytes(&self) -> usize {
        // SwissTable: (8 + 8)-byte slots + 1 control byte per bucket.
        let buckets = |cap: usize| {
            if cap == 0 {
                0
            } else if cap < 8 {
                if cap < 4 {
                    4
                } else {
                    8
                }
            } else {
                (cap * 8 / 7).next_power_of_two()
            }
        };
        let mut b = buckets(self.map.capacity()) * 17;
        b += buckets(self.more.capacity()) * 33;
        b += self.more.values().map(|v| v.capacity() * 8).sum::<usize>();
        b
    }

    pub(crate) fn shrink_to_fit(&mut self) {
        self.map.shrink_to_fit();
        self.more.shrink_to_fit();
    }
}

/// Aggregate of many [`QueryTrace`]s.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct QueryStats {
    pub queries: u64,
    pub hits: u64,
    pub entries_scanned: u64,
    pub tag_matches: u64,
    pub candidate_comparisons: u64,
    /// Number of queries that did at least one candidate comparison.
    pub queries_with_string_reads: u64,
    pub max_entries_scanned: u64,
    pub layer2: u64,
    pub layer3: u64,
}

impl QueryStats {
    /// Add the counters of `o` (e.g. another thread's).
    pub fn merge(&mut self, o: &QueryStats) {
        self.queries += o.queries;
        self.hits += o.hits;
        self.entries_scanned += o.entries_scanned;
        self.tag_matches += o.tag_matches;
        self.candidate_comparisons += o.candidate_comparisons;
        self.queries_with_string_reads += o.queries_with_string_reads;
        self.max_entries_scanned = self.max_entries_scanned.max(o.max_entries_scanned);
        self.layer2 += o.layer2;
        self.layer3 += o.layer3;
    }

    pub fn record(&mut self, hit: bool, t: &QueryTrace) {
        self.queries += 1;
        self.hits += hit as u64;
        self.entries_scanned += t.entries_scanned;
        self.tag_matches += t.tag_matches;
        self.candidate_comparisons += t.candidate_comparisons;
        self.queries_with_string_reads += (t.candidate_comparisons > 0) as u64;
        self.max_entries_scanned = self.max_entries_scanned.max(t.entries_scanned);
        self.layer2 += t.layer2;
        self.layer3 += t.layer3;
    }

    fn per_query(&self, x: u64) -> f64 {
        if self.queries == 0 {
            0.0
        } else {
            x as f64 / self.queries as f64
        }
    }
}

impl fmt::Display for QueryStats {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "  queries                      {}", self.queries)?;
        writeln!(f, "  hits                         {}", self.hits)?;
        writeln!(
            f,
            "  entries scanned / query      {:.3}  (max {})",
            self.per_query(self.entries_scanned),
            self.max_entries_scanned
        )?;
        writeln!(
            f,
            "  tag matches / query          {:.3}",
            self.per_query(self.tag_matches)
        )?;
        writeln!(
            f,
            "  candidate k-mers / query     {:.3}",
            self.per_query(self.candidate_comparisons)
        )?;
        writeln!(
            f,
            "  queries with string reads    {:.4}",
            self.per_query(self.queries_with_string_reads)
        )?;
        write!(
            f,
            "  routed to layer 2 / 3        {:.4} / {:.4}",
            self.per_query(self.layer2),
            self.per_query(self.layer3)
        )
    }
}

/// Counters accumulated by [`KmerTable::insert_read`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct InsertStats {
    pub reads: u64,
    /// Maximal ACGT segments of length ≥ k.
    pub segments: u64,
    /// ACGT segments shorter than k (produce nothing).
    pub short_segments: u64,
    /// k-mers looked up during insertion.
    pub kmers_processed: u64,
    /// Lookups that found the k-mer already present (skipped).
    pub present_hits: u64,
    /// Pieces emitted because a run reached `max_piece_len`.
    pub cap_splits: u64,
    /// Concurrent insertion: pieces closed because a writer's base chunk
    /// was full (the run continues in a new piece, overlapping by k - 1).
    pub chunk_splits: u64,
    /// ℓ-minimizers declared hot (entries moved to layer 2).
    pub hot1: u64,
    /// ℓ2-minimizers declared hot (k-mers moved to layer 3).
    pub hot2: u64,
    /// Lookup work done during insertion.
    pub lookup: QueryStats,
}

impl InsertStats {
    /// Add the counters of `o` (e.g. another thread's).
    pub fn merge(&mut self, o: &InsertStats) {
        self.reads += o.reads;
        self.segments += o.segments;
        self.short_segments += o.short_segments;
        self.kmers_processed += o.kmers_processed;
        self.present_hits += o.present_hits;
        self.cap_splits += o.cap_splits;
        self.chunk_splits += o.chunk_splits;
        self.hot1 += o.hot1;
        self.hot2 += o.hot2;
        self.lookup.merge(&o.lookup);
    }
}

/// Snapshot of table-level counters and memory.
#[derive(Clone, Debug)]
pub struct TableStats {
    pub k: usize,
    pub l: usize,
    pub pieces: u64,
    pub total_bases: u64,
    /// Σ over pieces of `len - k + 1`; exceeds the number of distinct
    /// k-mers only by intra-run repeats (see [`KmerTable::redundant_kmers`]).
    pub stored_kmers: u64,
    pub entries: u64,
    /// Current bucket count (the index grows automatically).
    pub buckets: u64,
    /// Overflow blocks in the flat index.
    pub overflow_blocks: u64,
    /// Number of times the index grew.
    pub index_grows: u32,
    /// Fallback minimizer length (0 if routing is disabled).
    pub l2: usize,
    /// Layer-2 (ℓ2) index entries and buckets.
    pub entries2: u64,
    pub buckets2: u64,
    /// Full k-mers in layer 3.
    pub layer3_kmers: u64,
    /// Heap bytes (by capacity) of the layer-2 index and layer-3 set.
    pub index2_bytes_reserved: usize,
    pub layer3_bytes_reserved: usize,
    /// Growth spill traffic of both indexes.
    pub spill: SpillStats,
    /// Σ piece lengths (`total_bases` minus dead bases, except that bases
    /// shared by two pieces after a deletion count twice).
    pub live_bases: u64,
    /// Pieces emptied by deletion (ids reclaimed by defragmentation).
    pub empty_pieces: u64,
    /// Deletion / defragmentation counters.
    pub edits: crate::edit::EditStats,
    /// Heap bytes (by capacity) of the bitset of pieces changed since the
    /// last defragmentation (1 bit per piece id).
    pub dirty_bytes_reserved: usize,
    pub insert: InsertStats,
    /// Heap bytes (by length) of piece storage / index.
    pub piece_bytes: usize,
    pub index_bytes: usize,
    /// Heap bytes (by capacity) of piece storage / index.
    pub piece_bytes_reserved: usize,
    pub index_bytes_reserved: usize,
}

impl TableStats {
    pub fn mean_piece_len(&self) -> f64 {
        if self.pieces == 0 {
            0.0
        } else {
            self.live_bases as f64 / (self.pieces - self.empty_pieces).max(1) as f64
        }
    }
}

impl fmt::Display for TableStats {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let i = &self.insert;
        writeln!(f, "table (k={}, ℓ={})", self.k, self.l)?;
        writeln!(f, "  reads                        {}", i.reads)?;
        writeln!(
            f,
            "  segments (≥k / <k)           {} / {}",
            i.segments, i.short_segments
        )?;
        writeln!(f, "  k-mers processed             {}", i.kmers_processed)?;
        writeln!(f, "  present on insert (skipped)  {}", i.present_hits)?;
        writeln!(
            f,
            "  pieces                       {}  ({} emptied by deletion)",
            self.pieces, self.empty_pieces
        )?;
        writeln!(
            f,
            "  pieces split at cap / chunk  {} / {}",
            i.cap_splits, i.chunk_splits
        )?;
        writeln!(
            f,
            "  total bases (array / live)   {} / {}",
            self.total_bases, self.live_bases
        )?;
        writeln!(
            f,
            "  mean piece length            {:.1}",
            self.mean_piece_len()
        )?;
        writeln!(f, "  stored k-mers                {}", self.stored_kmers)?;
        writeln!(f, "  index entries                {}", self.entries)?;
        writeln!(
            f,
            "  buckets (grew {}×)            {}",
            self.index_grows, self.buckets
        )?;
        writeln!(
            f,
            "  entries / bucket             {:.3}  (head load {:.3}, overflow blocks {})",
            self.entries as f64 / self.buckets.max(1) as f64,
            self.entries as f64 / (self.buckets.max(1) * crate::index::BLOCK_SLOTS as u64) as f64,
            self.overflow_blocks
        )?;
        writeln!(
            f,
            "  hot minimizers (ℓ / ℓ2)      {} / {}  (ℓ2 = {}; layer-2 entries {} in {} buckets; layer-3 k-mers {})",
            i.hot1, i.hot2, self.l2, self.entries2, self.buckets2, self.layer3_kmers
        )?;
        writeln!(
            f,
            "  index growth spill           {} rebuilds, {} B via disk, {} B via memory",
            self.spill.rebuilds, self.spill.bytes_to_disk, self.spill.bytes_in_memory
        )?;
        writeln!(
            f,
            "  piece bytes (used/reserved)  {} / {}",
            self.piece_bytes, self.piece_bytes_reserved
        )?;
        writeln!(
            f,
            "  index bytes (used/reserved)  {} / {}",
            self.index_bytes, self.index_bytes_reserved
        )?;
        if self.edits.delete_calls + self.edits.defrags + self.edits.incremental_defrags > 0 {
            writeln!(f, "{}", self.edits)?;
        }
        writeln!(f, "  insertion lookups:")?;
        write!(f, "{}", i.lookup)
    }
}

/// Streaming insertion state. Its memory is bounded by `max_piece_len`
/// (the pending run's bases plus the ℓ-mer values of the piece being
/// indexed), independently of read or chromosome length.
#[derive(Clone, Debug)]
pub(crate) struct Stream {
    /// Rolling k-mer and ℓ-mer of the current ACGT segment.
    roll: RollingKmer,
    lroll: RollingKmer,
    /// Window minimum over the ℓ-mers of the current k-mer.
    wmin: WindowMin,
    /// Bases seen in the current ACGT segment.
    seg_len: u64,
    /// Bases of the pending run of absent k-mers (uppercase ACGT).
    pending: Vec<u8>,
    /// Number of k-mers in the pending run (0 = no run).
    run_kmers: usize,
    /// A read has been started with `insert_chunk` and not yet ended.
    pub(crate) in_read: bool,
    /// Scratch for indexing an emitted piece.
    pub(crate) vals: Vec<u64>,
    marks: Vec<bool>,
    positions: Vec<u32>,
    /// Window minimum over a piece's ℓ-mers (routing path).
    pwin: WindowMin,
    /// Starts of the piece's windows whose ℓ-minimizer is hot.
    hot_windows: Vec<(PieceId, u32)>,
    /// Minimizers whose insertion overflowed (checked for hotness).
    check: Vec<u64>,
    /// Eager run indexing: the open piece holding the current run.
    open: Option<PieceId>,
    /// Eager run indexing: bit `i` of `marks1` / `marks2` is set if
    /// position `mark_base + i` of the open piece is already filed in
    /// layer 1 / layer 2 (positions below the current window are final).
    mark_base: u32,
    marks1: u64,
    marks2: u64,
}

impl Stream {
    pub(crate) fn new(k: usize, l: usize) -> Self {
        Self {
            roll: RollingKmer::new(k),
            lroll: RollingKmer::new(l),
            wmin: WindowMin::new(k - l + 1),
            seg_len: 0,
            pending: Vec::new(),
            run_kmers: 0,
            in_read: false,
            vals: Vec::new(),
            marks: Vec::new(),
            positions: Vec::new(),
            pwin: WindowMin::new(k - l + 1),
            hot_windows: Vec::new(),
            check: Vec::new(),
            open: None,
            mark_base: 0,
            marks1: 0,
            marks2: 0,
        }
    }

    fn heap_bytes(&self) -> usize {
        self.pending.capacity()
            + self.vals.capacity() * 8
            + self.marks.capacity()
            + self.positions.capacity() * 4
            + self.hot_windows.capacity() * 8
            + self.check.capacity() * 8
            + 2 * 33 * 16 // WindowMin deques: at most w + 1 ≤ 33 (u64, u64) pairs
    }

    pub(crate) fn shrink_to_fit(&mut self) {
        self.pending.shrink_to_fit();
        self.vals = Vec::new();
        self.marks = Vec::new();
        self.positions = Vec::new();
        self.hot_windows = Vec::new();
        self.check = Vec::new();
    }
}

/// Exact canonical k-mer set stored as pieces, with a minimizer index.
#[derive(Clone, Debug)]
pub struct KmerTable {
    pub(crate) cfg: Config,
    /// Window size in ℓ-mers: `k - ℓ + 1`.
    pub(crate) w: usize,
    pub(crate) mask: u64,
    pub(crate) store: PieceStore,
    pub(crate) index: MinimizerIndex,
    /// Fallback minimizer length (0: routing disabled) and its window size.
    pub(crate) l2: usize,
    pub(crate) w2: usize,
    /// Layer 2: ℓ2-minimizer index over the same pieces, for the k-mers
    /// whose ℓ-minimizer is hot.
    pub(crate) index2: MinimizerIndex,
    /// Layer 3: full canonical k-mers whose ℓ- and ℓ2-minimizers are hot,
    /// with their piece locations.
    pub(crate) layer3: Layer3,
    pub(crate) stats: InsertStats,
    pub(crate) stored_kmers: u64,
    /// Deletion and defragmentation counters.
    pub(crate) edits: crate::edit::EditStats,
    /// Non-empty pieces right after the last `defrag`.
    pub(crate) pieces_at_defrag: u64,
    /// Bit per piece id: created or cut since the last defragmentation
    /// (candidates for `defrag_incremental`).
    pub(crate) dirty: Vec<u64>,
    /// Always `Some` between calls; taken out while inserting.
    pub(crate) stream: Option<Stream>,
}

impl KmerTable {
    /// Build an empty table. Fails with a descriptive error on bad parameters.
    pub fn new(cfg: Config) -> Result<Self, ConfigError> {
        cfg.validate()?;
        let l2 = cfg.fallback_l().unwrap_or(0);
        Ok(Self {
            w: cfg.k - cfg.l + 1,
            mask: kmer_mask(cfg.k),
            store: PieceStore::new(),
            index: MinimizerIndex::new(cfg.bucket_bits).with_spill_dir(cfg.spill_dir.clone()),
            l2,
            w2: if l2 > 0 { cfg.k - l2 + 1 } else { 0 },
            index2: MinimizerIndex::new(0).with_spill_dir(cfg.spill_dir.clone()),
            layer3: Layer3::default(),
            stats: InsertStats::default(),
            stored_kmers: 0,
            edits: Default::default(),
            pieces_at_defrag: 0,
            dirty: Vec::new(),
            stream: Some(Stream::new(cfg.k, cfg.l)),
            cfg,
        })
    }

    pub fn config(&self) -> &Config {
        &self.cfg
    }
    pub fn k(&self) -> usize {
        self.cfg.k
    }
    pub fn pieces(&self) -> &PieceStore {
        &self.store
    }
    pub fn index(&self) -> &MinimizerIndex {
        &self.index
    }
    /// The layer-2 (ℓ2) fallback index.
    pub fn index2(&self) -> &MinimizerIndex {
        &self.index2
    }
    /// Number of full k-mers held in layer 3.
    pub fn layer3_len(&self) -> usize {
        self.layer3.len()
    }
    /// Window size in ℓ-mers (`k - ℓ + 1`).
    pub fn window(&self) -> usize {
        self.w
    }
    pub fn insert_stats(&self) -> &InsertStats {
        &self.stats
    }

    /// Heap bytes held by the insertion buffers (pending run, piece ℓ-mer
    /// values, deque). Bounded by about `14 × max_piece_len` bytes whatever
    /// the input length.
    pub fn insert_buffer_bytes(&self) -> usize {
        self.stream.as_ref().map_or(0, Stream::heap_bytes)
    }

    // ------------------------------------------------------------------
    // Insertion, streaming
    // ------------------------------------------------------------------

    /// Insert all k-mers of `read`. The read is split at non-ACGT characters;
    /// lowercase `acgt` is accepted and stored uppercase. Segments shorter
    /// than k contribute nothing. Equivalent to `insert_chunk(read)` followed
    /// by `end_read()`.
    pub fn insert_read(&mut self, read: &[u8]) {
        self.insert_chunk(read);
        self.end_read();
    }

    /// Feed the next chunk of the current read (or chromosome). Chunks may be
    /// cut anywhere, including inside a k-mer or a run of `N`s: the result is
    /// identical to inserting the concatenation with [`insert_read`]. Memory
    /// does not grow with the read length. Call [`end_read`] after the last
    /// chunk; until then, with [`RunIndexing::Deferred`], the last pending
    /// run of the read is not yet visible to queries.
    ///
    /// [`insert_read`]: Self::insert_read
    /// [`end_read`]: Self::end_read
    pub fn insert_chunk(&mut self, chunk: &[u8]) {
        let mut st = self.stream.take().expect("insertion state present");
        st.in_read = true;
        let k = self.cfg.k;
        let max_run = self.cfg.max_piece_len as usize - k + 1; // k-mers per piece
        let eager = self.cfg.run_indexing == RunIndexing::Eager;

        for &c in chunk {
            let b = encode(c);
            if b == INVALID {
                self.end_segment(&mut st);
                continue;
            }
            st.seg_len += 1;
            st.roll.push(b);
            if st.lroll.push(b) {
                st.wmin.push(order_value(st.lroll.canonical()));
            }
            if !st.roll.is_full() {
                continue;
            }
            // The k-mer is full exactly when its w ℓ-mers have been pushed.
            let (m, ties) = st.wmin.current().expect("window complete");
            let x = st.roll.canonical();
            self.stats.kmers_processed += 1;

            // A run that already holds max_run k-mers is closed *before* the
            // lookup, so the lookup sees it indexed. This keeps every
            // duplicate stored k-mer inside a single piece.
            let mut split = false;
            if st.run_kmers >= max_run {
                self.end_run(&mut st);
                split = true;
            }

            let mut trace = QueryTrace::default();
            let present = self.contains_with_minimizer(x, m, ties, &mut trace);
            self.stats.lookup.record(present, &trace);

            if present {
                self.stats.present_hits += 1;
                if st.run_kmers > 0 {
                    self.end_run(&mut st);
                }
            } else if eager {
                self.stats.cap_splits += split as u64;
                self.eager_push(&mut st, b, m, ties);
            } else {
                if st.run_kmers == 0 {
                    // Start a run with the k bases of the current k-mer.
                    let f = st.roll.fwd();
                    st.pending.clear();
                    st.pending
                        .extend((0..k).map(|i| codec::decode((f >> (2 * (k - 1 - i))) as u8)));
                    // Count a split only when the run actually continues.
                    self.stats.cap_splits += split as u64;
                } else {
                    st.pending.push(codec::decode(b));
                }
                st.run_kmers += 1;
                debug_assert_eq!(st.pending.len(), st.run_kmers + k - 1);
            }
        }
        self.stream = Some(st);
    }

    /// Finish the current read: emit its last pending run and reset the
    /// k-mer window. A no-op if no read is open.
    pub fn end_read(&mut self) {
        let mut st = self.stream.take().expect("insertion state present");
        self.end_segment(&mut st);
        if st.in_read {
            self.stats.reads += 1;
            st.in_read = false;
        }
        self.stream = Some(st);
    }

    /// Close the current ACGT segment (non-ACGT character or end of read).
    fn end_segment(&mut self, st: &mut Stream) {
        if st.run_kmers > 0 {
            self.end_run(st);
        }
        if st.seg_len >= self.cfg.k as u64 {
            self.stats.segments += 1;
        } else if st.seg_len > 0 {
            self.stats.short_segments += 1;
        }
        st.seg_len = 0;
        st.roll.reset();
        st.lroll.reset();
        st.wmin.clear();
    }

    /// Minimizer of index-1 entry `(p, q)`: the ℓ-mer at `q`.
    #[inline]
    fn remin1(store: &PieceStore, l: usize, p: PieceId, q: u32) -> u64 {
        order_value(codec::canonicalize(store.kmer_at(p, q, l), l))
    }

    /// Close the current run of absent k-mers (`run_kmers > 0`).
    fn end_run(&mut self, st: &mut Stream) {
        match st.open.take() {
            Some(pid) => {
                self.stored_kmers += st.run_kmers as u64;
                self.edits.inserted_pieces_since_defrag += 1;
                self.mark_dirty(pid);
                st.run_kmers = 0;
            }
            None => self.emit_pending(st),
        }
    }

    /// Eager run indexing: append the absent k-mer just read (its last base
    /// is `b`, its window minimizer `(m, ties)`) to the open piece, opening
    /// one if needed, and file its window in the index right away.
    fn eager_push(&mut self, st: &mut Stream, b: u8, m: u64, ties: u32) {
        let k = self.cfg.k;
        let f = st.roll.fwd();
        let (pid, s0) = match st.open {
            Some(pid) => {
                self.store.push_code(b);
                self.store.grow_piece(pid);
                (pid, st.run_kmers as u32)
            }
            None => {
                let start = self.store.end();
                for i in 0..k {
                    self.store.push_code((f >> (2 * (k - 1 - i))) as u8);
                }
                let pid = self.store.close_piece(start);
                st.open = Some(pid);
                st.mark_base = 0;
                st.marks1 = 0;
                st.marks2 = 0;
                (pid, 0)
            }
        };
        st.run_kmers += 1;
        let d = s0 - st.mark_base;
        st.marks1 = st.marks1.checked_shr(d).unwrap_or(0);
        st.marks2 = st.marks2.checked_shr(d).unwrap_or(0);
        st.mark_base = s0;
        if self.is_hot(m) {
            self.eager_layer2(st, pid, s0, f);
            return;
        }
        let new = ties as u64 & !st.marks1;
        if new == 0 {
            return;
        }
        st.marks1 |= new;
        let l = self.cfg.l;
        let store = &self.store;
        self.index
            .reserve(new.count_ones() as u64, |p, q| Self::remin1(store, l, p, q));
        let mut overflowed = false;
        let mut bits = new;
        while bits != 0 {
            overflowed |= self.index.insert(m, pid, s0 + bits.trailing_zeros());
            bits &= bits - 1;
        }
        if overflowed && self.l2 > 0 {
            self.maybe_make_hot1(m);
        }
    }

    /// Eager run indexing of window `(pid, s0)` (forward k-mer `f`) whose
    /// ℓ-minimizer is hot: layer 2, or layer 3 if its ℓ2-minimizer is hot.
    fn eager_layer2(&mut self, st: &mut Stream, pid: PieceId, s0: u32, f: u64) {
        let (k, l2) = (self.cfg.k, self.l2);
        let (m2, ties) = kmer_minimizer_mask(f, k, l2);
        if self.index2.is_hot(m2) {
            self.layer3.add(codec::canonicalize(f, k), (pid, s0));
            return;
        }
        let new = ties as u64 & !st.marks2;
        if new == 0 {
            return;
        }
        st.marks2 |= new;
        let store = &self.store;
        self.index2.reserve(new.count_ones() as u64, |p, q| {
            order_value(codec::canonicalize(store.kmer_at(p, q, l2), l2))
        });
        let mut overflowed = false;
        let mut bits = new;
        while bits != 0 {
            overflowed |= self.index2.insert(m2, pid, s0 + bits.trailing_zeros());
            bits &= bits - 1;
        }
        if overflowed {
            self.maybe_make_hot2(m2);
        }
    }

    /// Store the pending run as a piece and index its distinct minimizer
    /// occurrences `(m, q)`. With hot-minimizer routing, windows
    /// whose ℓ-minimizer is hot are indexed in layer 2 (or 3) instead.
    fn emit_pending(&mut self, st: &mut Stream) {
        let pid = self.store.add_piece(&st.pending);
        self.stored_kmers += st.run_kmers as u64;
        self.edits.inserted_pieces_since_defrag += 1;
        self.mark_dirty(pid);
        lmer_values_into(&st.pending, self.cfg.l, &mut st.vals);
        self.index_piece(pid, st);
        st.pending.clear();
        st.run_kmers = 0;
    }

    /// Index every window of piece `pid`, whose ℓ-mer order values are in
    /// `st.vals`: layer 1 for cold ℓ-minimizers, layer 2 (or 3) for hot
    /// ones, then promote minimizers that became hot.
    pub(crate) fn index_piece(&mut self, pid: PieceId, st: &mut Stream) {
        let l = self.cfg.l;
        st.hot_windows.clear();
        if self.l2 == 0 || self.index.num_hot() == 0 {
            piece_minimizer_positions_into(&st.vals, self.w, &mut st.marks, &mut st.positions);
        } else {
            // Per window: hot ℓ-minimizer → layer 2; otherwise mark its ties.
            st.marks.clear();
            st.marks.resize(st.vals.len(), false);
            st.pwin.clear();
            let mut last: Option<(u64, bool)> = None;
            for (i, &v) in st.vals.iter().enumerate() {
                st.pwin.push(v);
                let Some((m, ties)) = st.pwin.current() else {
                    continue;
                };
                let s0 = i + 1 - self.w;
                let hot = match last {
                    Some((lm, h)) if lm == m => h,
                    _ => {
                        let h = self.index.is_hot(m);
                        last = Some((m, h));
                        h
                    }
                };
                if hot {
                    st.hot_windows.push((pid, s0 as u32));
                } else {
                    let mut bits = ties;
                    while bits != 0 {
                        st.marks[s0 + bits.trailing_zeros() as usize] = true;
                        bits &= bits - 1;
                    }
                }
            }
            st.positions.clear();
            st.positions.extend(
                st.marks
                    .iter()
                    .enumerate()
                    .filter(|(_, &b)| b)
                    .map(|(p, _)| p as u32),
            );
        }
        // Grow the flat index first if needed; entries are re-placed using
        // minimizers recomputed from the (packed) pieces.
        let store = &self.store;
        self.index.reserve(st.positions.len() as u64, |p, q| {
            Self::remin1(store, l, p, q)
        });
        st.check.clear();
        for &q in &st.positions {
            let m = st.vals[q as usize];
            if self.index.insert(m, pid, q) && self.l2 > 0 {
                st.check.push(m);
            }
        }
        if !st.hot_windows.is_empty() {
            let wins = std::mem::take(&mut st.hot_windows);
            self.index_layer2(&wins);
            st.hot_windows = wins;
        }
        // Minimizers that overflowed both heads may have become hot.
        if !st.check.is_empty() {
            st.check.sort_unstable();
            st.check.dedup();
            let check = std::mem::take(&mut st.check);
            for &m in &check {
                self.maybe_make_hot1(m);
            }
            st.check = check;
        }
    }

    /// Index the given windows `(piece, start)` (k-mers whose ℓ-minimizer is
    /// hot) by their ℓ2-minimizers, or put them in layer 3 if that is hot too.
    fn index_layer2(&mut self, windows: &[(PieceId, u32)]) {
        let (k, l2) = (self.cfg.k, self.l2);
        let mut entries: Vec<(u64, PieceId, u32)> = Vec::new();
        for &(pid, s0) in windows {
            let f = self.store.kmer_at(pid, s0, k);
            let (m2, ties) = kmer_minimizer_mask(f, k, l2);
            if self.index2.is_hot(m2) {
                self.layer3.add(codec::canonicalize(f, k), (pid, s0));
                continue;
            }
            let mut bits = ties;
            while bits != 0 {
                entries.push((m2, pid, s0 + bits.trailing_zeros()));
                bits &= bits - 1;
            }
        }
        entries.sort_unstable_by_key(|&(_, p, q)| (p, q));
        entries.dedup();
        let store = &self.store;
        self.index2.reserve(entries.len() as u64, |p, q| {
            order_value(codec::canonicalize(store.kmer_at(p, q, l2), l2))
        });
        let mut check = Vec::new();
        for &(m2, pid, q) in &entries {
            if self.index2.insert(m2, pid, q) {
                check.push(m2);
            }
        }
        check.sort_unstable();
        check.dedup();
        for m2 in check {
            self.maybe_make_hot2(m2);
        }
    }

    /// Declare ℓ-minimizer `m` hot if it has reached the threshold: drain its
    /// entries from index 1 and re-index the k-mers they cover in layer 2.
    fn maybe_make_hot1(&mut self, m: u64) {
        let (k, l) = (self.cfg.k, self.cfg.l);
        if self.index.is_hot(m) {
            return;
        }
        let store = &self.store;
        let remin = |p, q| Self::remin1(store, l, p, q);
        if (self.index.count_minimizer(m, remin) as u32) < self.cfg.hot_threshold {
            return;
        }
        let drained = self.index.drain_minimizer(m, remin);
        self.stats.hot1 += 1;
        // Every window whose minimizer is m had all its tie positions
        // indexed, so the drained entries enumerate all of them.
        let mut wins: Vec<(PieceId, u32)> = Vec::new();
        for (pid, q) in drained {
            let len = self.store.length(pid) as usize;
            let q = q as usize;
            for s0 in q.saturating_sub(k - l)..=q.min(len - k) {
                let f = self.store.kmer_at(pid, s0 as u32, k);
                let (mm, ties) = kmer_minimizer_mask(f, k, l);
                if mm == m && ties >> (q - s0) & 1 == 1 {
                    wins.push((pid, s0 as u32));
                }
            }
        }
        wins.sort_unstable();
        wins.dedup();
        self.index_layer2(&wins);
    }

    /// Declare ℓ2-minimizer `m2` hot if it has reached the threshold: drain
    /// its entries from index 2 and put the k-mers they cover in layer 3.
    fn maybe_make_hot2(&mut self, m2: u64) {
        let (k, l, l2) = (self.cfg.k, self.cfg.l, self.l2);
        if self.index2.is_hot(m2) {
            return;
        }
        let store = &self.store;
        let remin = |p, q| order_value(codec::canonicalize(store.kmer_at(p, q, l2), l2));
        if (self.index2.count_minimizer(m2, remin) as u32) < self.cfg.hot_threshold {
            return;
        }
        let drained = self.index2.drain_minimizer(m2, remin);
        self.stats.hot2 += 1;
        for (pid, q) in drained {
            let len = self.store.length(pid) as usize;
            let q = q as usize;
            for s0 in q.saturating_sub(k - l2)..=q.min(len - k) {
                let f = self.store.kmer_at(pid, s0 as u32, k);
                let (mm2, ties) = kmer_minimizer_mask(f, k, l2);
                if mm2 == m2
                    && ties >> (q - s0) & 1 == 1
                    && self.index.is_hot(kmer_minimizer_mask(f, k, l).0)
                {
                    self.layer3.add(codec::canonicalize(f, k), (pid, s0 as u32));
                }
            }
        }
    }

    // ------------------------------------------------------------------
    // Queries
    // ------------------------------------------------------------------

    /// Exact membership of a packed k-mer (either orientation; it is
    /// canonicalised internally). Values with bits above `2k` return false.
    #[inline]
    pub fn contains(&self, kmer: u64) -> bool {
        self.contains_traced(kmer, &mut QueryTrace::default())
    }

    /// [`contains`](Self::contains) that also reports the work done.
    pub fn contains_traced(&self, kmer: u64, trace: &mut QueryTrace) -> bool {
        if kmer & !self.mask != 0 {
            return false;
        }
        let k = self.cfg.k;
        let x = codec::canonicalize(kmer, k);
        let (m, ties) = kmer_minimizer_mask(x, k, self.cfg.l);
        self.contains_with_minimizer(x, m, ties, trace)
    }

    /// Membership of an ASCII k-mer (case-insensitive, either strand).
    pub fn contains_seq(&self, kmer: &[u8]) -> Result<bool, KmerError> {
        if kmer.len() != self.cfg.k {
            return Err(KmerError::WrongLength {
                expected: self.cfg.k,
                got: kmer.len(),
            });
        }
        let mut x = 0u64;
        for (pos, &c) in kmer.iter().enumerate() {
            let b = encode(c);
            if b == INVALID {
                return Err(KmerError::InvalidBase { pos, byte: c });
            }
            x = (x << 2) | b as u64;
        }
        Ok(self.contains(x))
    }

    /// Core lookup: `x` is canonical, `m` is its window minimizer value and
    /// `ties` the bit mask of its tied minimizer positions `Q(x)`.
    fn contains_with_minimizer(&self, x: u64, m: u64, ties: u32, trace: &mut QueryTrace) -> bool {
        self.find_with_minimizer(x, m, ties, trace).is_some()
    }

    /// Lookup in layers 2 and 3, for a canonical `x` whose ℓ-minimizer is hot.
    pub(crate) fn find_fallback(&self, x: u64, trace: &mut QueryTrace) -> Option<Hit> {
        let (k, l2) = (self.cfg.k, self.l2);
        let (m2, ties) = kmer_minimizer_mask(x, k, l2);
        if self.index2.is_hot(m2) {
            trace.layer3 += 1;
            return self.layer3.get(x).map(|(piece, pos)| Hit {
                piece,
                pos,
                kmer: self.store.kmer_at(piece, pos, k),
            });
        }
        trace.layer2 += 1;
        let ties = Self::mirror(ties, self.w2);
        let xr = codec::revcomp_kmer(x, k);
        let mut comparisons = 0u64;
        let mut scan = ScanTrace::default();
        let mut found = None;
        self.index2.for_each_candidate(m2, &mut scan, |pid, q| {
            found = self.check_candidate(x, xr, ties, pid, q, l2, &mut comparisons);
            found.is_some()
        });
        trace.add_scan(&scan);
        trace.candidate_comparisons += comparisons;
        found
    }

    /// Like [`contains_with_minimizer`](Self::contains_with_minimizer) but
    /// returns where the k-mer was found.
    pub(crate) fn find_with_minimizer(
        &self,
        x: u64,
        m: u64,
        ties: u32,
        trace: &mut QueryTrace,
    ) -> Option<Hit> {
        if self.l2 > 0 && self.index.is_hot(m) {
            return self.find_fallback(x, trace);
        }
        let ties = self.mirrored_ties(ties);
        let xr = codec::revcomp_kmer(x, self.cfg.k);
        let mut comparisons = 0u64;
        let mut scan = ScanTrace::default();
        let mut found = None;
        let l = self.cfg.l;
        self.index.for_each_candidate(m, &mut scan, |pid, q| {
            found = self.check_candidate(x, xr, ties, pid, q, l, &mut comparisons);
            found.is_some()
        });
        trace.add_scan(&scan);
        trace.candidate_comparisons += comparisons;
        found
    }

    /// The piece may hold x in either orientation; on the reverse strand tie
    /// position p becomes (k - ℓ) - p, so take the union of Q(x) and its
    /// mirror image (w = k - ℓ + 1 ≤ 32 bits).
    #[inline]
    pub(crate) fn mirrored_ties(&self, ties: u32) -> u32 {
        Self::mirror(ties, self.w)
    }

    #[inline(always)]
    fn mirror(ties: u32, w: usize) -> u32 {
        ties | (ties.reverse_bits() >> (32 - w))
    }

    /// Whether routing sends ℓ-minimizer `m` to the fallback layers.
    #[inline(always)]
    pub(crate) fn is_hot(&self, m: u64) -> bool {
        self.l2 > 0 && self.index.is_hot(m)
    }

    /// Exact check of one index entry `(pid, q)` for canonical k-mer `x`
    /// (`xr` is its reverse complement). `ties` must already include the
    /// mirror image (see `mirrored_ties`). A stored k-mer `f` matches iff
    /// `f == x` or `f == xr`, so no per-candidate reverse complement is needed.
    #[inline]
    #[allow(clippy::too_many_arguments)] // hot path: plain scalars, no struct
    pub(crate) fn check_candidate(
        &self,
        x: u64,
        xr: u64,
        ties: u32,
        pid: PieceId,
        q: u32,
        l: usize,
        comparisons: &mut u64,
    ) -> Option<Hit> {
        let k = self.cfg.k;
        let q = q as usize;
        let (g0, len) = self.store.view(pid);
        let len = len as usize;
        let mut test = |s: usize| {
            *comparisons += 1;
            let f = self.store.kmer_at_global(g0 + s as u64, k);
            (f == x || f == xr).then_some(Hit {
                piece: pid,
                pos: s as u32,
                kmer: f,
            })
        };
        match self.cfg.scan {
            CandidateScan::Range => {
                let span = k - l; // = w - 1
                let s_lo = q.saturating_sub(span);
                if len < k {
                    return None;
                }
                let s_hi = q.min(len - k);
                (s_lo..=s_hi).find_map(test)
            }
            CandidateScan::TieOffsets => {
                let mut bits = ties;
                while bits != 0 {
                    let p = bits.trailing_zeros() as usize;
                    bits &= bits - 1;
                    if p > q {
                        break; // higher p only make s more negative
                    }
                    let s = q - p;
                    if s + k > len {
                        continue;
                    }
                    if let Some(h) = test(s) {
                        return Some(h);
                    }
                }
                None
            }
        }
    }

    /// Answer the k-mers of a sequence in order, reusing work between
    /// successive k-mers (see [`StreamingQuery`]).
    pub fn query_stream(&self) -> StreamingQuery<'_> {
        StreamingQuery::new(self)
    }

    // ------------------------------------------------------------------
    // Introspection
    // ------------------------------------------------------------------

    pub fn stats(&self) -> TableStats {
        let (pu, pr) = self.store.heap_bytes();
        let (iu, ir) = self.index.heap_bytes();
        TableStats {
            k: self.cfg.k,
            l: self.cfg.l,
            pieces: self.store.len() as u64,
            total_bases: self.store.total_bases(),
            stored_kmers: self.stored_kmers,
            entries: self.index.num_entries(),
            buckets: self.index.num_buckets() as u64,
            overflow_blocks: self.index.overflow_blocks() as u64,
            index_grows: self.index.grows(),
            l2: self.l2,
            entries2: self.index2.num_entries(),
            buckets2: self.index2.num_buckets() as u64,
            layer3_kmers: self.layer3.len() as u64,
            index2_bytes_reserved: self.index2.heap_bytes().1,
            layer3_bytes_reserved: self.layer3.heap_bytes(),
            live_bases: self.store.live_bases(),
            empty_pieces: self.store.empty_pieces(),
            edits: self.edits,
            dirty_bytes_reserved: self.dirty.capacity() * 8,
            spill: {
                let (a, b) = (self.index.spill_stats(), self.index2.spill_stats());
                SpillStats {
                    rebuilds: a.rebuilds + b.rebuilds,
                    bytes_to_disk: a.bytes_to_disk + b.bytes_to_disk,
                    bytes_in_memory: a.bytes_in_memory + b.bytes_in_memory,
                }
            },
            insert: self.stats,
            piece_bytes: pu,
            index_bytes: iu,
            piece_bytes_reserved: pr,
            index_bytes_reserved: ir,
        }
    }

    /// Visit every stored canonical k-mer occurrence (piece order).
    pub fn for_each_stored_kmer<F: FnMut(PieceId, u32, u64)>(&self, mut f: F) {
        let k = self.cfg.k;
        for pid in 0..self.store.len() as PieceId {
            let len = self.store.length(pid) as usize;
            if len < k {
                continue;
            }
            for s in 0..=(len - k) {
                f(
                    pid,
                    s as u32,
                    self.store.canonical_kmer_at(pid, s as u32, k),
                );
            }
        }
    }

    /// Number of stored k-mer occurrences that duplicate another stored
    /// occurrence (intra-run repeats). Diagnostic, O(n).
    pub fn redundant_kmers(&self) -> u64 {
        let mut seen = std::collections::HashSet::with_capacity(self.stored_kmers as usize);
        let mut dup = 0u64;
        self.for_each_stored_kmer(|_, _, x| {
            if !seen.insert(x) {
                dup += 1;
            }
        });
        dup
    }

    /// Number of distinct k-mers stored (diagnostic, O(n) memory).
    pub fn distinct_kmers(&self) -> u64 {
        self.stored_kmers - self.redundant_kmers()
    }

    /// Minimizer value of every index entry, recomputed from its piece.
    /// Diagnostic for skew histograms.
    pub fn entry_minimizers(&self) -> impl Iterator<Item = u64> + '_ {
        let l = self.cfg.l;
        self.index
            .iter_entries()
            .map(move |e| order_value(self.store.canonical_kmer_at(e.piece_id, e.q, l)))
    }

    /// Release unused capacity (useful before measuring memory).
    pub fn shrink_to_fit(&mut self) {
        self.store.shrink_to_fit();
        self.index.shrink_to_fit();
        self.index2.shrink_to_fit();
        self.layer3.shrink_to_fit();
        self.dirty.shrink_to_fit();
        if let Some(st) = self.stream.as_mut() {
            st.shrink_to_fit();
        }
    }
}
