//! Concurrent insertion, deletion and queries.
//!
//! [`ConcurrentKmerTable`] holds the same k-mer set representation as
//! [`KmerTable`] (packed pieces, minimizer index with hot-minimizer layers)
//! but can be shared between threads: any number of [`Writer`]s insert reads
//! while other threads call [`ConcurrentKmerTable::delete_read`] /
//! [`ConcurrentKmerTable::delete`] and [`ConcurrentKmerTable::contains`].
//! [`ConcurrentKmerTable::query_stream`] answers the k-mers of sequences in
//! order (see [`ConcurrentStreamingQuery`]). Defragmentation needs exclusive
//! access and is done in place ([`ConcurrentKmerTable::defrag`],
//! [`ConcurrentKmerTable::defrag_incremental`]);
//! [`ConcurrentKmerTable::into_table`] and
//! [`ConcurrentKmerTable::from_table`] convert to and from a [`KmerTable`].
//! Deletion, defragmentation and streaming queries are described in the
//! `concurrent_delete.rs`, `concurrent_defrag.rs` and `concurrent_stream.rs`
//! parts of this module.
//!
//! # How it works
//!
//! * **Pieces.** The base array and the piece records are segmented arrays of
//!   atomics whose segments never move, so reads need no lock. Each writer
//!   appends to its own chunk of 2^16 bases; a piece record (`start << 24 |
//!   len`) is one `u64`, published with release ordering after its bases.
//!   Insertion never modifies an existing piece, it only appends.
//! * **Index.** Every index layer is split into shards by a hash of the
//!   minimizer (layer 3 by a hash of the k-mer), each a [`MinimizerIndex`]
//!   behind its own `RwLock`, growing on its own.
//! * **Exactness.** A k-mer is looked up and, if absent, appended to the
//!   writer's open piece and filed in the index *while the writer holds the
//!   write lock of the k-mer's home shard*: the shard of its minimizer, or of
//!   its ℓ2-minimizer if the minimizer is hot, or its layer-3 shard. Every
//!   occurrence of a k-mer is filed in its home shard, so two writers cannot
//!   both store it (insertion is eager, see [`crate::RunIndexing::Eager`]).
//!   A writer keeps the lock while consecutive k-mers share the shard.
//! * **Hot minimizers.** A minimizer becomes hot under its layer-1 shard's
//!   write lock: its k-mers are first filed in layer 2 (or 3), then it is
//!   marked hot and its layer-1 entries are dropped, so a lookup never misses
//!   a k-mer in transit. Hotness is permanent, so writers cache it.
//! * **No false positives without locking pieces.** Every hit is checked
//!   against the bases of the piece's current view, so it is always a true
//!   member. Locks are always taken in the order layer 1 → layer 2 → layer 3,
//!   and a writer holds at most one shard lock between k-mers: no deadlock.

use crate::codec::{self, encode, RollingKmer, INVALID};
use crate::edit::EditStats;
use crate::index::{MinimizerIndex, ScanTrace};
use crate::minimizer::{kmer_minimizer_mask, order_value, WindowMin};
use crate::piece_store::{pack, PieceId, PieceStore, LEN_BITS, LEN_MASK};
use crate::seg::SegArray;
use crate::table::{
    CandidateScan, Config, ConfigError, Hit, InsertStats, KmerTable, Layer3, QueryTrace,
};
use std::collections::HashSet;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Mutex, PoisonError, RwLock, RwLockReadGuard, RwLockWriteGuard};
use std::time::Instant;

/// Bases per writer chunk (a multiple of 32: chunks never share a word).
const CHUNK_BASES: u64 = 1 << 16;
/// Base array: up to 2^35 words (2^40 bases); first segments of 2^12 words
/// (a multiple of a chunk, so a chunk never spans two segments), 16 per
/// size.
type BaseArray = SegArray<12, 4>;
const MAX_WORDS: u64 = 1 << 35;
/// Piece records: up to 2^32; first segments of 2^10, 16 per size.
type MetaArray = SegArray<10, 4>;
const MAX_PIECES: u64 = 1 << 32;
const LAYER3_SHARDS: usize = 64;
/// Piece ids handed to a writer at a time (writers update their pieces'
/// records constantly: separate blocks keep them on separate cache lines).
const PID_BLOCK: u64 = 16;
const SHARD_SEED: u64 = 0x9e37_79b9_7f4a_7c15;

/// Layer-1 shards used by [`ConcurrentKmerTable::new`].
pub const DEFAULT_SHARDS: usize = 1024;

#[inline(always)]
fn mix(mut x: u64) -> u64 {
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^ (x >> 31)
}

#[inline(always)]
fn shard_of(x: u64, n: usize) -> usize {
    (((mix(x ^ SHARD_SEED) >> 32) * n as u64) >> 32) as usize
}

/// Hasher for the writers' caches of hot minimizers (keys are already
/// well-mixed order values).
#[derive(Default)]
struct MixHasher(u64);

impl std::hash::Hasher for MixHasher {
    fn finish(&self) -> u64 {
        mix(self.0)
    }
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0 << 8) | b as u64;
        }
    }
    fn write_u64(&mut self, x: u64) {
        self.0 = x;
    }
}

type FastSet = HashSet<u64, std::hash::BuildHasherDefault<MixHasher>>;

fn wlock<T>(l: &RwLock<T>) -> RwLockWriteGuard<'_, T> {
    l.write().unwrap_or_else(PoisonError::into_inner)
}

fn rlock<T>(l: &RwLock<T>) -> RwLockReadGuard<'_, T> {
    l.read().unwrap_or_else(PoisonError::into_inner)
}

/// Where writers spend time holding or waiting for shard locks.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ContentionStats {
    /// Shard write locks taken by writers (layers 1–3).
    pub acquisitions: u64,
    /// Of those, how many were held by another thread at first try.
    pub contended: u64,
    /// Total time writers waited for a contended lock (ns, summed over
    /// threads).
    pub wait_ns: u64,
    /// Contended acquisitions that still had to sleep after spinning.
    pub parked: u64,
    /// Index shard rebuilds (growth), and their time (ns, under the lock).
    pub growths: u64,
    pub growth_ns: u64,
    /// Time spent checking and making minimizers hot (ns, under the lock).
    pub promotion_ns: u64,
    /// Deletions: piece locks found held by another cut, and the time spent
    /// waiting for them (ns).
    pub piece_contended: u64,
    pub piece_wait_ns: u64,
    /// Deletions: cuts that found their piece still open (a writer's run,
    /// or a right part being filed), and the time spent backing off (ns).
    pub open_waits: u64,
    pub open_wait_ns: u64,
    /// Deletions: cuts whose piece changed since the lookup (looked up
    /// again).
    pub stale_retries: u64,
}

impl ContentionStats {
    fn merge(&mut self, o: &ContentionStats) {
        self.acquisitions += o.acquisitions;
        self.contended += o.contended;
        self.wait_ns += o.wait_ns;
        self.parked += o.parked;
        self.growths += o.growths;
        self.growth_ns += o.growth_ns;
        self.promotion_ns += o.promotion_ns;
        self.piece_contended += o.piece_contended;
        self.piece_wait_ns += o.piece_wait_ns;
        self.open_waits += o.open_waits;
        self.open_wait_ns += o.open_wait_ns;
        self.stale_retries += o.stale_retries;
    }

    /// `self - before` (counters accumulated between two snapshots).
    pub fn since(&self, before: &ContentionStats) -> ContentionStats {
        ContentionStats {
            acquisitions: self.acquisitions - before.acquisitions,
            contended: self.contended - before.contended,
            wait_ns: self.wait_ns - before.wait_ns,
            parked: self.parked - before.parked,
            growths: self.growths - before.growths,
            growth_ns: self.growth_ns - before.growth_ns,
            promotion_ns: self.promotion_ns - before.promotion_ns,
            piece_contended: self.piece_contended - before.piece_contended,
            piece_wait_ns: self.piece_wait_ns - before.piece_wait_ns,
            open_waits: self.open_waits - before.open_waits,
            open_wait_ns: self.open_wait_ns - before.open_wait_ns,
            stale_retries: self.stale_retries - before.stale_retries,
        }
    }
}

impl std::fmt::Display for ContentionStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{} lock acquisitions, {} contended ({:.3} %, {} slept), {:.3} s waiting; {} shard growths, {:.3} s; promotions {:.3} s (thread-seconds)",
            self.acquisitions,
            self.contended,
            100.0 * self.contended as f64 / self.acquisitions.max(1) as f64,
            self.parked,
            self.wait_ns as f64 * 1e-9,
            self.growths,
            self.growth_ns as f64 * 1e-9,
            self.promotion_ns as f64 * 1e-9
        )?;
        if self.piece_contended + self.open_waits + self.stale_retries > 0 {
            write!(
                f,
                "; deletions: {} piece-lock waits ({:.3} s), {} waits on open pieces ({:.3} s), {} stale retries",
                self.piece_contended,
                self.piece_wait_ns as f64 * 1e-9,
                self.open_waits,
                self.open_wait_ns as f64 * 1e-9,
                self.stale_retries
            )?;
        }
        Ok(())
    }
}

/// A shard on its own cache lines (writers on neighbouring shards do not
/// invalidate each other's lock word and index header).
#[repr(align(128))]
pub(crate) struct Pad<T>(T);

impl<T> std::ops::Deref for Pad<T> {
    type Target = T;
    fn deref(&self) -> &T {
        &self.0
    }
}

impl<T> std::ops::DerefMut for Pad<T> {
    fn deref_mut(&mut self) -> &mut T {
        &mut self.0
    }
}

/// `try_write` attempts before a writer sleeps on a contended shard lock.
const SPIN: u32 = if cfg!(test) { 64 } else { 4000 };

/// Write-lock `l`, counting contention.
fn acquire<'b, T>(l: &'b RwLock<T>, c: &mut ContentionStats) -> RwLockWriteGuard<'b, T> {
    c.acquisitions += 1;
    match l.try_write() {
        Ok(g) => g,
        Err(std::sync::TryLockError::Poisoned(p)) => p.into_inner(),
        Err(std::sync::TryLockError::WouldBlock) => {
            c.contended += 1;
            let t0 = Instant::now();
            // Shard locks are held for about a microsecond: spin before
            // parking the thread (a futex sleep and wake-up costs tens of µs).
            for _ in 0..SPIN {
                std::hint::spin_loop();
                if let Ok(g) = l.try_write() {
                    c.wait_ns += t0.elapsed().as_nanos() as u64;
                    return g;
                }
            }
            let g = wlock(l);
            c.wait_ns += t0.elapsed().as_nanos() as u64;
            c.parked += 1;
            g
        }
    }
}

/// `g.reserve(n)`, timing it when it rebuilds.
fn reserve_timed<F: Fn(PieceId, u32) -> u64>(
    g: &mut MinimizerIndex,
    n: u64,
    remin: F,
    c: &mut ContentionStats,
) {
    if g.needs_growth(n) {
        let t0 = Instant::now();
        g.reserve(n, remin);
        c.growths += 1;
        c.growth_ns += t0.elapsed().as_nanos() as u64;
    }
}

fn lock<T>(l: &Mutex<T>) -> std::sync::MutexGuard<'_, T> {
    l.lock().unwrap_or_else(PoisonError::into_inner)
}

/// Piece storage shared by all threads (see the module documentation).
struct SharedStore {
    data: BaseArray,
    meta: MetaArray,
    /// First base never handed out in a chunk.
    frontier: AtomicU64,
    /// Piece ids handed out (`meta` length).
    n_pids: AtomicU64,
    /// Unused ends of the chunks and id blocks of finished writers.
    tails: Mutex<Vec<(u64, u64)>>,
    pid_tails: Mutex<Vec<(u64, u64)>>,
    /// Ids of emptied pieces (from a converted table), reused first.
    free: Mutex<Vec<PieceId>>,
    free_left: AtomicBool,
    /// Ids taken from `free` (marked changed on conversion back).
    reused: Mutex<Vec<PieceId>>,
    /// Bit per piece id: the piece is a writer's open run (it may still
    /// grow, so deletion waits for it to close).
    open: OpenBits,
}

type OpenBits = SegArray<8, 4>;

impl SharedStore {
    fn new() -> Self {
        Self {
            data: BaseArray::new(MAX_WORDS),
            meta: MetaArray::new(MAX_PIECES),
            frontier: AtomicU64::new(0),
            n_pids: AtomicU64::new(0),
            tails: Mutex::new(Vec::new()),
            pid_tails: Mutex::new(Vec::new()),
            free: Mutex::new(Vec::new()),
            free_left: AtomicBool::new(false),
            reused: Mutex::new(Vec::new()),
            open: OpenBits::new(MAX_PIECES / 64),
        }
    }

    /// Mark piece `pid` open (before anything of it is published).
    #[inline]
    fn set_open(&self, pid: PieceId) {
        self.open
            .slot(pid as u64 >> 6)
            .fetch_or(1 << (pid & 63), Ordering::Release);
    }

    /// Mark piece `pid` closed (after its last length update).
    #[inline]
    fn clear_open(&self, pid: PieceId) {
        self.open
            .slot(pid as u64 >> 6)
            .fetch_and(!(1 << (pid & 63)), Ordering::Release);
    }

    #[inline]
    fn is_open(&self, pid: PieceId) -> bool {
        self.open.load(pid as u64 >> 6, Ordering::Acquire) >> (pid & 63) & 1 == 1
    }

    /// 2-bit code of the base at global offset `g`.
    #[inline(always)]
    fn code_at_global(&self, g: u64) -> u8 {
        (self.data.load(g >> 5, Ordering::Relaxed) >> (62 - 2 * (g & 31)) & 3) as u8
    }

    #[inline(always)]
    fn set_view(&self, pid: PieceId, start: u64, len: u32) {
        self.meta
            .slot(pid as u64)
            .store(pack(start, len), Ordering::Release);
    }

    /// A piece id for a deletion's right part.
    fn alloc_pid(&self) -> PieceId {
        if let Some(id) = self.reuse_pid() {
            return id;
        }
        let id = self.n_pids.fetch_add(1, Ordering::Relaxed);
        assert!(
            id < u32::MAX as u64,
            "piece store full: piece ids are u32 (max {} pieces)",
            u32::MAX
        );
        id as PieceId
    }

    /// Give an emptied id back for reuse.
    fn free_pid(&self, pid: PieceId) {
        lock(&self.free).push(pid);
        self.free_left.store(true, Ordering::Relaxed);
    }

    #[inline(always)]
    fn view(&self, pid: PieceId) -> (u64, u32) {
        let m = self.meta.load(pid as u64, Ordering::Acquire);
        (m >> LEN_BITS, (m & LEN_MASK) as u32)
    }

    /// Forward packed k-mer at global base offset `g` (same layout as
    /// [`PieceStore::kmer_at_global`]).
    #[inline(always)]
    fn kmer_at_global(&self, g: u64, k: usize) -> u64 {
        let i = g >> 5;
        let r = (g & 31) as u32;
        let (hi, lo) = if 2 * r + 2 * k as u32 > 64 {
            self.data.load2(i)
        } else {
            (self.data.load(i, Ordering::Relaxed), 0)
        };
        let v = ((hi as u128) << 64) | lo as u128;
        ((v << (2 * r)) >> (128 - 2 * k as u32)) as u64
    }

    #[inline(always)]
    fn kmer_at(&self, pid: PieceId, s: u32, k: usize) -> u64 {
        self.kmer_at_global(self.view(pid).0 + s as u64, k)
    }

    /// An emptied id to reuse, if any.
    fn reuse_pid(&self) -> Option<PieceId> {
        if !self.free_left.load(Ordering::Relaxed) {
            return None;
        }
        let mut f = lock(&self.free);
        let id = f.pop()?;
        if f.is_empty() {
            self.free_left.store(false, Ordering::Relaxed);
        }
        lock(&self.reused).push(id);
        Some(id)
    }

    /// A block of fresh piece ids.
    fn alloc_pids(&self) -> (u64, u64) {
        if let Some(b) = lock(&self.pid_tails).pop() {
            return b;
        }
        let id = self.n_pids.fetch_add(PID_BLOCK, Ordering::Relaxed);
        assert!(
            id + PID_BLOCK <= u32::MAX as u64,
            "piece store full: piece ids are u32 (max {} pieces)",
            u32::MAX
        );
        (id, id + PID_BLOCK)
    }

    /// A chunk with room for at least `need` bases.
    fn alloc_chunk(&self, need: u64) -> (u64, u64) {
        {
            let mut t = lock(&self.tails);
            if let Some(i) = t.iter().position(|&(s, e)| e - s >= need) {
                return t.swap_remove(i);
            }
        }
        let s = self.frontier.fetch_add(CHUNK_BASES, Ordering::Relaxed);
        assert!(
            s + CHUNK_BASES <= self.data.capacity() * 32,
            "piece store full: at most 2^40 bases"
        );
        (s, s + CHUNK_BASES)
    }

    fn heap_bytes(&self) -> usize {
        self.data.heap_bytes()
            + self.meta.heap_bytes()
            + (lock(&self.tails).capacity() + lock(&self.pid_tails).capacity()) * 16
            + (lock(&self.free).capacity() + lock(&self.reused).capacity()) * 4
            + self.open.heap_bytes()
    }
}

/// Read access to piece views and bases, for candidate checks.
trait Bases {
    fn view(&self, pid: PieceId) -> (u64, u32);
    fn kmer_at_global(&self, g: u64, k: usize) -> u64;
}

impl Bases for SharedStore {
    #[inline(always)]
    fn view(&self, pid: PieceId) -> (u64, u32) {
        SharedStore::view(self, pid)
    }
    #[inline(always)]
    fn kmer_at_global(&self, g: u64, k: usize) -> u64 {
        SharedStore::kmer_at_global(self, g, k)
    }
}

/// A writer's view of the store: its open piece and its chunk's segment
/// are read directly (in eager insertion most candidates are in them).
struct Local<'s> {
    store: &'s SharedStore,
    open: Option<(PieceId, u64, u32)>,
    seg_first: u64,
    seg: &'s [AtomicU64],
}

impl Bases for Local<'_> {
    #[inline(always)]
    fn view(&self, pid: PieceId) -> (u64, u32) {
        match self.open {
            Some((p, s, l)) if p == pid => (s, l),
            _ => self.store.view(pid),
        }
    }
    #[inline(always)]
    fn kmer_at_global(&self, g: u64, k: usize) -> u64 {
        let i = (g >> 5).wrapping_sub(self.seg_first);
        if i < (self.seg.len() as u64).saturating_sub(1) {
            let r = (g & 31) as u32;
            let hi = self.seg[i as usize].load(Ordering::Relaxed);
            let lo = self.seg[i as usize + 1].load(Ordering::Relaxed);
            let v = ((hi as u128) << 64) | lo as u128;
            ((v << (2 * r)) >> (128 - 2 * k as u32)) as u64
        } else {
            self.store.kmer_at_global(g, k)
        }
    }
}

/// Counters carried over from / back to a [`KmerTable`], and merged from
/// finished writers.
#[derive(Default)]
struct Carry {
    stats: InsertStats,
    stored_kmers: u64,
    pieces_created: u64,
    edits: EditStats,
    pieces_at_defrag: u64,
    dirty: Vec<u64>,
    /// Piece ids below this existed before (the rest were created here).
    pids_at_start: u64,
    contention: ContentionStats,
    /// Pieces cut by deletions (marked changed on conversion back).
    changed: Vec<PieceId>,
}

/// Snapshot of a [`ConcurrentKmerTable`] (counters of finished writers).
#[derive(Clone, Debug, Default)]
pub struct ConcurrentStats {
    pub shards: usize,
    pub pieces: u64,
    /// Bases handed out to writers' chunks (live, dead and not yet used).
    pub reserved_bases: u64,
    pub live_bases: u64,
    pub stored_kmers: u64,
    pub entries: u64,
    pub entries2: u64,
    pub layer3_kmers: u64,
    pub hot1: u64,
    pub hot2: u64,
    /// Heap bytes (by capacity) of the whole table.
    pub heap_bytes: usize,
    pub insert: InsertStats,
    pub contention: ContentionStats,
    /// Deletion counters (deletions merge them when each call returns).
    pub edits: EditStats,
}

/// Exact canonical k-mer set that several threads can insert into and query
/// at the same time. See the [module documentation](self).
pub struct ConcurrentKmerTable {
    cfg: Config,
    w: usize,
    l2: usize,
    w2: usize,
    mask: u64,
    store: SharedStore,
    l1: Box<[Pad<RwLock<MinimizerIndex>>]>,
    l2s: Box<[Pad<RwLock<MinimizerIndex>>]>,
    l3: Box<[Pad<RwLock<Layer3>>]>,
    /// Union of the layer-1 hot minimizers (read when a layer-2 minimizer
    /// becomes hot; never locked while another lock is taken).
    hot1: RwLock<HashSet<u64>>,
    /// Union of the layer-2 hot minimizers.
    hot2: RwLock<HashSet<u64>>,
    /// Held (shared) by each deletion cut; minimizers become hot only when
    /// it can be taken exclusively (a promotion is skipped otherwise and
    /// retried later), so hotness is frozen while a cut computes its items.
    gate: RwLock<()>,
    /// Striped locks serializing cuts of the same piece.
    piece_locks: Box<[Pad<Mutex<()>>]>,
    /// The table may hold k-mers stored twice (converted from a table built
    /// with deferred run indexing): deletions re-check their input.
    may_have_duplicates: bool,
    carry: Mutex<Carry>,
}

/// Piece lock stripes.
const PIECE_LOCKS: usize = 1024;

impl ConcurrentKmerTable {
    /// An empty table with [`DEFAULT_SHARDS`] layer-1 shards.
    pub fn new(cfg: Config) -> Result<Self, ConfigError> {
        Self::with_shards(cfg, DEFAULT_SHARDS)
    }

    /// An empty table whose layer-1 index has `shards` shards (layer 2 has
    /// `shards / 8`, layer 3 has 64). More shards mean less waiting between
    /// writers; each costs about 250 B.
    pub fn with_shards(cfg: Config, shards: usize) -> Result<Self, ConfigError> {
        let shards = shards.max(1);
        let buckets = 1usize
            << cfg
                .bucket_bits
                .saturating_sub(shards.next_power_of_two().trailing_zeros());
        Self::build(cfg, shards, buckets)
    }

    /// An empty table sized for about `kmers` distinct k-mers: the index
    /// shards start large enough not to grow (a shard grows under its write
    /// lock, which stalls every writer that needs it meanwhile), and there
    /// are more shards for larger tables (fewer writers meeting on one
    /// lock): `kmers / 8192` rounded up to a power of two, between 256 and
    /// 65536. An over-estimate costs about 1.6 B per k-mer of excess index
    /// space; an under-estimate only means that shards grow as usual.
    pub fn with_capacity(cfg: Config, kmers: u64) -> Result<Self, ConfigError> {
        let shards = ((kmers / 8192) as usize)
            .next_power_of_two()
            .clamp(256, 65536);
        // Index entries per k-mer ≈ 2 / (w + 1) (density of a random
        // minimizer scheme) plus about one per piece: 0.128 on chr21, 0.151
        // on a read set (k = 31, ℓ = 17). Size for 2.5 / (w + 1) at a head
        // load of 80 %: read sets fill the heads to about 77 % (growth starts
        // at 85 %, shards vary by a few %), genomes to about 66 %. Under-sizing
        // costs more than it saves: a shard that grows ends at ×1.5 its size.
        let w = (cfg.k + 1).saturating_sub(cfg.l).max(1) as f64;
        let entries = kmers as f64 * 2.5 / (w + 1.0);
        let buckets = (entries / (crate::index::BLOCK_SLOTS as f64 * 0.8) / shards as f64).ceil();
        Self::build(cfg, shards, (buckets as usize).max(1))
    }

    fn build(cfg: Config, shards: usize, buckets: usize) -> Result<Self, ConfigError> {
        cfg.validate()?;
        let l2 = cfg.fallback_l().unwrap_or(0);
        let mk = |n: usize, buckets: usize| -> Box<[Pad<RwLock<MinimizerIndex>>]> {
            (0..n)
                .map(|_| {
                    Pad(RwLock::new(
                        MinimizerIndex::with_buckets(buckets).with_spill_dir(cfg.spill_dir.clone()),
                    ))
                })
                .collect()
        };
        Ok(Self {
            w: cfg.k - cfg.l + 1,
            l2,
            w2: if l2 > 0 { cfg.k - l2 + 1 } else { 0 },
            mask: codec::kmer_mask(cfg.k),
            store: SharedStore::new(),
            l1: mk(shards, buckets),
            l2s: mk((shards / 8).max(1), 1),
            l3: (0..LAYER3_SHARDS)
                .map(|_| Pad(RwLock::new(Layer3::default())))
                .collect(),
            hot1: RwLock::new(HashSet::new()),
            hot2: RwLock::new(HashSet::new()),
            gate: RwLock::new(()),
            piece_locks: (0..PIECE_LOCKS).map(|_| Pad(Mutex::new(()))).collect(),
            may_have_duplicates: false,
            carry: Mutex::new(Carry::default()),
            cfg,
        })
    }

    pub fn config(&self) -> &Config {
        &self.cfg
    }

    pub fn k(&self) -> usize {
        self.cfg.k
    }

    /// A handle that inserts reads from the current thread. Create one per
    /// thread; dropping it ends its open read and records its counters.
    pub fn writer(&self) -> Writer<'_> {
        Writer::new(self)
    }

    // ------------------------------------------------------------------
    // Queries
    // ------------------------------------------------------------------

    /// Exact membership of a packed k-mer (either orientation), safe to call
    /// while writers insert. A k-mer inserted by a writer is visible as soon
    /// as that writer has processed it.
    pub fn contains(&self, kmer: u64) -> bool {
        self.contains_traced(kmer, &mut QueryTrace::default())
    }

    /// [`contains`](Self::contains) that also reports the work done.
    pub fn contains_traced(&self, kmer: u64, trace: &mut QueryTrace) -> bool {
        kmer & !self.mask == 0
            && self
                .locate_canonical(codec::canonicalize(kmer, self.cfg.k), trace)
                .is_some()
    }

    /// Where a k-mer (either orientation) is stored, if present.
    pub fn locate(&self, kmer: u64) -> Option<Hit> {
        if kmer & !self.mask != 0 {
            return None;
        }
        self.locate_canonical(
            codec::canonicalize(kmer, self.cfg.k),
            &mut QueryTrace::default(),
        )
    }

    pub(crate) fn locate_canonical(&self, x: u64, trace: &mut QueryTrace) -> Option<Hit> {
        let (k, l) = (self.cfg.k, self.cfg.l);
        let (m, ties) = kmer_minimizer_mask(x, k, l);
        {
            let g = rlock(&self.l1[shard_of(m, self.l1.len())]);
            if !(self.l2 > 0 && g.is_hot(m)) {
                return self.find(&self.store, &g, x, m, ties, l, self.w, trace);
            }
        }
        self.locate_hot(x, trace)
    }

    /// Where canonical `x`, whose ℓ-minimizer is hot, is stored (layers 2
    /// and 3), if present.
    pub(crate) fn locate_hot(&self, x: u64, trace: &mut QueryTrace) -> Option<Hit> {
        let k = self.cfg.k;
        let (m2, ties2) = kmer_minimizer_mask(x, k, self.l2);
        trace.layer2 += 1;
        {
            let g = rlock(&self.l2s[shard_of(m2, self.l2s.len())]);
            if !g.is_hot(m2) {
                return self.find(&self.store, &g, x, m2, ties2, self.l2, self.w2, trace);
            }
        }
        trace.layer2 -= 1;
        trace.layer3 += 1;
        let loc = rlock(&self.l3[shard_of(x, LAYER3_SHARDS)]).get(x);
        loc.map(|(piece, pos)| Hit {
            piece,
            pos,
            kmer: self.store.kmer_at(piece, pos, k),
        })
    }

    /// Look `x` (canonical) up in one locked index layer.
    #[allow(clippy::too_many_arguments)]
    fn find<B: Bases>(
        &self,
        bases: &B,
        idx: &MinimizerIndex,
        x: u64,
        m: u64,
        ties: u32,
        l: usize,
        w: usize,
        trace: &mut QueryTrace,
    ) -> Option<Hit> {
        let k = self.cfg.k;
        let ties = ties | (ties.reverse_bits() >> (32 - w));
        let xr = codec::revcomp_kmer(x, k);
        let mut comparisons = 0u64;
        let mut scan = ScanTrace::default();
        let mut found = None;
        idx.for_each_candidate(m, &mut scan, |pid, q| {
            found = self.check_candidate(bases, x, xr, ties, pid, q, l, &mut comparisons);
            found.is_some()
        });
        trace.add_scan(&scan);
        trace.candidate_comparisons += comparisons;
        found
    }

    /// [`find`](Self::find) over a list of candidates (collected from a
    /// locked shard).
    #[allow(clippy::too_many_arguments)]
    fn find_in<B: Bases>(
        &self,
        bases: &B,
        cands: &[(PieceId, u32)],
        x: u64,
        ties: u32,
        l: usize,
        w: usize,
        trace: &mut QueryTrace,
    ) -> Option<Hit> {
        let ties = ties | (ties.reverse_bits() >> (32 - w));
        let xr = codec::revcomp_kmer(x, self.cfg.k);
        let mut comparisons = 0u64;
        let found = cands.iter().find_map(|&(pid, q)| {
            self.check_candidate(bases, x, xr, ties, pid, q, l, &mut comparisons)
        });
        trace.candidate_comparisons += comparisons;
        found
    }

    /// The extension cursor of hit `h` (forward k-mer `fwd`), if the bases
    /// at the position computed from the piece's current view are the hit's
    /// (bases never change once written, so the cursor stays valid as long
    /// as some view covers it).
    ///
    /// Deletion never moves the start of a piece it keeps, so the check
    /// fails only when the piece was emptied and its id reused for a new
    /// piece between the lookup and this call (possible for a hit found in
    /// layers 2/3, without the layer-1 shard lock). That window is too
    /// narrow for the tests to hit; the check keeps the cursor exact there.
    #[inline]
    fn cursor<B: Bases>(&self, bases: &B, h: &Hit, fwd: u64) -> Option<Cursor> {
        let g = bases.view(h.piece).0 + h.pos as u64;
        (bases.kmer_at_global(g, self.cfg.k) == h.kmer).then_some(Cursor {
            pid: h.piece,
            g,
            same: h.kmer == fwd,
        })
    }

    /// The next k-mer (last base `b`) is right after the cursor's (or,
    /// reverse-complemented, right before) in the cursor's piece, whose
    /// current view still covers it: then it is present (at the moment the
    /// view is read). Returns the moved cursor.
    #[inline]
    fn extend<B: Bases>(&self, bases: &B, c: Cursor, b: u8) -> Option<Cursor> {
        let k = self.cfg.k as u64;
        let (lo, len) = bases.view(c.pid);
        let hi = lo + len as u64;
        let (g, at, want) = if c.same {
            (c.g + 1, c.g + k, b)
        } else {
            let g = c.g.checked_sub(1)?;
            (g, g, 3 - b)
        };
        (g >= lo && g + k <= hi && self.store.code_at_global(at) == want)
            .then_some(Cursor { g, ..c })
    }

    /// Same check as [`KmerTable`]'s: every hit is read from the piece's
    /// current view, so it is a true member.
    #[allow(clippy::too_many_arguments)]
    #[inline]
    fn check_candidate<B: Bases>(
        &self,
        bases: &B,
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
        let (g0, len) = bases.view(pid);
        let len = len as usize;
        if len < k {
            return None;
        }
        let mut test = |s: usize| {
            *comparisons += 1;
            let f = bases.kmer_at_global(g0 + s as u64, k);
            (f == x || f == xr).then_some(Hit {
                piece: pid,
                pos: s as u32,
                kmer: f,
            })
        };
        match self.cfg.scan {
            CandidateScan::Range => (q.saturating_sub(k - l)..=q.min(len - k)).find_map(test),
            CandidateScan::TieOffsets => {
                let mut bits = ties;
                while bits != 0 {
                    let p = bits.trailing_zeros() as usize;
                    bits &= bits - 1;
                    if p > q {
                        break;
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

    // ------------------------------------------------------------------
    // Hot minimizers (called with the shard's write lock held)
    // ------------------------------------------------------------------

    fn remin(&self, l: usize) -> impl Fn(PieceId, u32) -> u64 + '_ {
        move |p, q| order_value(codec::canonicalize(self.store.kmer_at(p, q, l), l))
    }

    /// Layer-1 shard `g` holds minimizer `m`: make it hot if it reached the
    /// threshold (file its k-mers in layer 2/3 first, then drop its entries).
    fn promote1(&self, g: &mut MinimizerIndex, m: u64, stats: &mut InsertStats) {
        let (k, l) = (self.cfg.k, self.cfg.l);
        if g.is_hot(m) || (g.count_minimizer(m, self.remin(l)) as u32) < self.cfg.hot_threshold {
            return;
        }
        // Not while a deletion is cutting a piece (retried at the next
        // overflow of m).
        let Ok(_gate) = self.gate.try_write() else {
            return;
        };
        // Collect the windows before draining (the shard lock is held, so
        // nobody observes the intermediate state).
        let mut wins: Vec<(PieceId, u32)> = Vec::new();
        g.for_each_candidate(m, &mut ScanTrace::default(), |pid, q| {
            if self.remin(l)(pid, q) != m {
                return false;
            }
            let len = self.store.view(pid).1 as usize;
            let q = q as usize;
            for s0 in q.saturating_sub(k - l)..=q.min(len - k) {
                let f = self.store.kmer_at(pid, s0 as u32, k);
                let (mm, ties) = kmer_minimizer_mask(f, k, l);
                if mm == m && ties >> (q - s0) & 1 == 1 {
                    wins.push((pid, s0 as u32));
                }
            }
            false
        });
        wins.sort_unstable();
        wins.dedup();
        stats.hot1 += 1;
        wlock(&self.hot1).insert(m);
        self.index_layer2(&wins, stats);
        g.drain_minimizer(m, self.remin(l));
    }

    /// File windows `(piece, start)` (k-mers whose ℓ-minimizer is hot) in
    /// layer 2 by their ℓ2-minimizers, or in layer 3 if that is hot too.
    fn index_layer2(&self, windows: &[(PieceId, u32)], stats: &mut InsertStats) {
        let (k, l2) = (self.cfg.k, self.l2);
        let n = self.l2s.len();
        let mut items: Vec<(usize, u64, PieceId, u32, u32)> = windows
            .iter()
            .map(|&(pid, s0)| {
                let f = self.store.kmer_at(pid, s0, k);
                let (m2, ties) = kmer_minimizer_mask(f, k, l2);
                (shard_of(m2, n), m2, pid, s0, ties)
            })
            .collect();
        items.sort_unstable_by_key(|it| it.0);
        let mut entries: Vec<(u64, PieceId, u32)> = Vec::new();
        let mut i = 0;
        while i < items.len() {
            let j = i + items[i..]
                .iter()
                .take_while(|it| it.0 == items[i].0)
                .count();
            let group = &items[i..j];
            i = j;
            let mut g = wlock(&self.l2s[group[0].0]);
            entries.clear();
            for &(_, m2, pid, s0, ties) in group {
                if g.is_hot(m2) {
                    let x = codec::canonicalize(self.store.kmer_at(pid, s0, k), k);
                    wlock(&self.l3[shard_of(x, LAYER3_SHARDS)]).add(x, (pid, s0));
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
            g.reserve(entries.len() as u64, self.remin(l2));
            let mut check = Vec::new();
            for &(m2, pid, q) in &entries {
                if g.insert(m2, pid, q) {
                    check.push(m2);
                }
            }
            check.sort_unstable();
            check.dedup();
            for m2 in check {
                self.promote2_gated(&mut g, m2, stats);
            }
        }
    }

    /// Layer-2 shard `g` holds ℓ2-minimizer `m2`: make it hot if it reached
    /// the threshold (k-mers to layer 3 first, then drop its entries).
    fn promote2(&self, g: &mut MinimizerIndex, m2: u64, stats: &mut InsertStats) {
        if let Ok(_gate) = self.gate.try_write() {
            self.promote2_gated(g, m2, stats);
        }
    }

    /// [`promote2`](Self::promote2) with the gate already held.
    fn promote2_gated(&self, g: &mut MinimizerIndex, m2: u64, stats: &mut InsertStats) {
        let (k, l, l2) = (self.cfg.k, self.cfg.l, self.l2);
        if g.is_hot(m2) || (g.count_minimizer(m2, self.remin(l2)) as u32) < self.cfg.hot_threshold {
            return;
        }
        let mut add: Vec<(u64, PieceId, u32)> = Vec::new();
        {
            let hot1 = rlock(&self.hot1);
            g.for_each_candidate(m2, &mut ScanTrace::default(), |pid, q| {
                if self.remin(l2)(pid, q) != m2 {
                    return false;
                }
                let len = self.store.view(pid).1 as usize;
                let q = q as usize;
                for s0 in q.saturating_sub(k - l2)..=q.min(len - k) {
                    let f = self.store.kmer_at(pid, s0 as u32, k);
                    let (mm2, ties) = kmer_minimizer_mask(f, k, l2);
                    if mm2 == m2
                        && ties >> (q - s0) & 1 == 1
                        && hot1.contains(&kmer_minimizer_mask(f, k, l).0)
                    {
                        add.push((codec::canonicalize(f, k), pid, s0 as u32));
                    }
                }
                false
            });
        }
        for (x, pid, s0) in add {
            wlock(&self.l3[shard_of(x, LAYER3_SHARDS)]).add(x, (pid, s0));
        }
        stats.hot2 += 1;
        wlock(&self.hot2).insert(m2);
        g.drain_minimizer(m2, self.remin(l2));
    }

    // ------------------------------------------------------------------
    // Introspection and conversion
    // ------------------------------------------------------------------

    /// Counters and memory. Locks every shard in turn (for reading); the
    /// insertion counters include finished writers only.
    pub fn stats(&self) -> ConcurrentStats {
        let n_pids = self.store.n_pids.load(Ordering::Relaxed);
        let mut st = ConcurrentStats {
            shards: self.l1.len(),
            reserved_bases: self.store.frontier.load(Ordering::Relaxed),
            ..Default::default()
        };
        for pid in 0..n_pids {
            let len = self.store.view(pid as PieceId).1 as u64;
            st.pieces += (len > 0) as u64;
            st.live_bases += len;
        }
        let mut heap = self.store.heap_bytes() + hashset_bytes(rlock(&self.hot1).capacity());
        for s in self.l1.iter() {
            let g = rlock(s);
            st.entries += g.num_entries();
            st.hot1 += g.num_hot() as u64;
            heap += g.heap_bytes().1;
        }
        for s in self.l2s.iter() {
            let g = rlock(s);
            st.entries2 += g.num_entries();
            st.hot2 += g.num_hot() as u64;
            heap += g.heap_bytes().1;
        }
        for s in self.l3.iter() {
            let g = rlock(s);
            st.layer3_kmers += g.len() as u64;
            heap += g.heap_bytes();
        }
        heap += (self.l1.len() + self.l2s.len())
            * std::mem::size_of::<Pad<RwLock<MinimizerIndex>>>()
            + self.l3.len() * std::mem::size_of::<Pad<RwLock<Layer3>>>();
        st.heap_bytes = heap;
        let c = lock(&self.carry);
        st.stored_kmers = c.stored_kmers;
        st.insert = c.stats;
        st.contention = c.contention;
        st.edits = c.edits;
        st
    }

    /// Hand the content over to a single-threaded [`KmerTable`] (for
    /// deletion, defragmentation, streaming queries, or lock-free queries
    /// from many threads). Unused chunk space is compacted away and the
    /// index shards are merged; pieces created here count as changed for
    /// [`KmerTable::defrag_incremental`].
    pub fn into_table(self) -> KmerTable {
        let Self {
            cfg,
            store,
            l1,
            l2s,
            l3,
            carry,
            l2,
            ..
        } = self;
        let carry = carry.into_inner().unwrap_or_else(PoisonError::into_inner);
        let n_pids = store.n_pids.load(Ordering::Relaxed);
        let n_bases = store.frontier.load(Ordering::Relaxed);
        let data = store
            .data
            .to_vec((n_bases as usize).div_ceil(32) as u64 + 1);
        let meta = store.meta.to_vec(n_pids);
        let free = store
            .free
            .into_inner()
            .unwrap_or_else(PoisonError::into_inner);
        let reused = store
            .reused
            .into_inner()
            .unwrap_or_else(PoisonError::into_inner);
        drop(store.data);
        drop(store.meta);
        let mut ps = PieceStore::from_raw(data, n_bases, meta, free);
        ps.compact_bases();
        ps.shrink_to_fit();

        let mut t = KmerTable::new(cfg).expect("validated configuration");
        t.store = ps;
        let (k, l) = (t.cfg.k, t.cfg.l);
        t.index = merge_shards(l1, &t.store, l, &t.cfg);
        if l2 > 0 {
            t.index2 = merge_shards(l2s, &t.store, l2, &t.cfg);
        }
        for s in l3.into_vec() {
            let s = s.0.into_inner().unwrap_or_else(PoisonError::into_inner);
            for (x, loc) in s.iter() {
                debug_assert_eq!(codec::canonicalize(t.store.kmer_at(loc.0, loc.1, k), k), x);
                t.layer3.add(x, loc);
            }
        }
        t.stats = carry.stats;
        t.stored_kmers = carry.stored_kmers;
        t.edits = carry.edits;
        t.edits.inserted_pieces_since_defrag += carry.pieces_created;
        t.pieces_at_defrag = carry.pieces_at_defrag;
        t.dirty = carry.dirty;
        for pid in carry.pids_at_start..n_pids {
            t.mark_dirty(pid as PieceId);
        }
        for pid in reused.into_iter().chain(carry.changed) {
            if (pid as usize) < t.store.len() && t.store.length(pid) > 0 {
                t.mark_dirty(pid);
            }
        }
        t.shrink_to_fit();
        t
    }

    /// Take over a [`KmerTable`] (no read may be open) to insert into it
    /// from several threads, with [`DEFAULT_SHARDS`] shards.
    pub fn from_table(t: KmerTable) -> Self {
        Self::from_table_with_shards(t, DEFAULT_SHARDS)
    }

    /// [`from_table`](Self::from_table) with `shards` layer-1 shards.
    pub fn from_table_with_shards(t: KmerTable, shards: usize) -> Self {
        assert!(
            !t.stream.as_ref().is_some_and(|s| s.in_read),
            "from_table: a read is open (call end_read first)"
        );
        let KmerTable {
            cfg,
            store,
            index,
            index2,
            layer3,
            stats,
            stored_kmers,
            edits,
            pieces_at_defrag,
            dirty,
            ..
        } = t;
        let mut c = Self::with_shards(cfg, shards).expect("validated configuration");
        c.may_have_duplicates = c.cfg.run_indexing == crate::table::RunIndexing::Deferred;
        let (data, n_bases, meta, free) = store.into_raw();
        c.store
            .data
            .fill_from(&data[..(n_bases as usize).div_ceil(32)]);
        drop(data);
        c.store.meta.fill_from(&meta);
        let n_pids = meta.len() as u64;
        drop(meta);
        // Chunks start at multiples of their size (never across segments).
        c.store.frontier.store(
            n_bases.div_ceil(CHUNK_BASES) * CHUNK_BASES,
            Ordering::Relaxed,
        );
        c.store.n_pids.store(n_pids, Ordering::Relaxed);
        c.store.free_left.store(!free.is_empty(), Ordering::Relaxed);
        *lock(&c.store.free) = free;

        let (l, l2) = (c.cfg.l, c.l2);
        let hot: Vec<u64> = index.hot_minimizers().copied().collect();
        split_into_shards(index, &c.store, l, &mut c.l1);
        for &m in &hot {
            c.l1[shard_of(m, c.l1.len())]
                .get_mut()
                .unwrap_or_else(PoisonError::into_inner)
                .mark_hot(m);
        }
        *c.hot1.get_mut().unwrap_or_else(PoisonError::into_inner) = hot.into_iter().collect();
        if l2 > 0 {
            let hot2: Vec<u64> = index2.hot_minimizers().copied().collect();
            split_into_shards(index2, &c.store, l2, &mut c.l2s);
            for &m2 in &hot2 {
                c.l2s[shard_of(m2, c.l2s.len())]
                    .get_mut()
                    .unwrap_or_else(PoisonError::into_inner)
                    .mark_hot(m2);
            }
            *c.hot2.get_mut().unwrap_or_else(PoisonError::into_inner) = hot2.into_iter().collect();
        }
        for (x, loc) in layer3.iter() {
            c.l3[shard_of(x, LAYER3_SHARDS)]
                .get_mut()
                .unwrap_or_else(PoisonError::into_inner)
                .add(x, loc);
        }
        *c.carry.get_mut().unwrap_or_else(PoisonError::into_inner) = Carry {
            stats,
            stored_kmers,
            pieces_created: 0,
            edits,
            pieces_at_defrag,
            dirty,
            pids_at_start: n_pids,
            contention: ContentionStats::default(),
            changed: Vec::new(),
        };
        c
    }
}

fn hashset_bytes(cap: usize) -> usize {
    crate::index::hashset_u64_bytes(cap)
}

/// Position `j` of a scattered order of `0..n` (`p` = `n` rounded up to a
/// power of two): a fixed odd stride modulo `p`, skipping values `≥ n`.
/// Entries listed bucket by bucket and re-inserted in that order would defeat
/// two-choice placement and overflow far more.
#[inline(always)]
fn scattered_at(j: usize, p: usize) -> usize {
    const STRIDE: usize = 0x9e37_79b9; // odd: a permutation modulo p
    j.wrapping_mul(STRIDE) & (p - 1)
}

/// Threads used by conversions: all available cores, or one for small jobs.
fn conversion_threads(work: usize) -> usize {
    if work < 1 << 16 {
        1
    } else {
        std::thread::available_parallelism().map_or(1, |n| n.get())
    }
}

/// Work unit of the parallel conversion loops (scattered positions).
const CONVERT_CHUNK: usize = 1 << 18;

/// Run `f(start, end)` over `0..n` split into chunks, on `threads` threads.
fn par_ranges<F: Fn(usize, usize) + Sync>(n: usize, chunk: usize, threads: usize, f: F) {
    if threads <= 1 || n <= chunk {
        let mut s = 0;
        while s < n {
            f(s, (s + chunk).min(n));
            s += chunk;
        }
        return;
    }
    let next = AtomicU64::new(0);
    std::thread::scope(|sc| {
        for _ in 0..threads {
            let (next, f) = (&next, &f);
            sc.spawn(move || loop {
                let s = next.fetch_add(chunk as u64, Ordering::Relaxed) as usize;
                if s >= n {
                    break;
                }
                f(s, (s + chunk).min(n));
            });
        }
    });
}

/// Move every entry of `index` into the shard of its minimizer. The shards
/// are re-sized first for their share of the entries (no growth), then
/// filled by all cores in a scattered order.
fn split_into_shards(
    index: MinimizerIndex,
    store: &SharedStore,
    l: usize,
    shards: &mut [Pad<RwLock<MinimizerIndex>>],
) {
    let list: Vec<(PieceId, u32)> = index.iter_entries().map(|e| (e.piece_id, e.q)).collect();
    drop(index);
    let n = shards.len();
    // Shard sizes are close to the mean (thousands of entries each); size
    // for a 75 % head load.
    let per = list.len() as f64 / n as f64;
    let buckets = ((per / (crate::index::BLOCK_SLOTS as f64 * 0.75)).ceil() as usize).max(1);
    for s in shards.iter_mut() {
        let g = s.get_mut().unwrap_or_else(PoisonError::into_inner);
        let dir = g.spill_dir_owned();
        *g = MinimizerIndex::with_buckets(buckets).with_spill_dir(dir);
    }
    let remin = |p: PieceId, q: u32| order_value(codec::canonicalize(store.kmer_at(p, q, l), l));
    let p = list.len().next_power_of_two().max(1);
    let shards: &[Pad<RwLock<MinimizerIndex>>] = shards;
    par_ranges(p, CONVERT_CHUNK, conversion_threads(list.len()), |a, b| {
        let mut c = ContentionStats::default();
        for j in a..b {
            let i = scattered_at(j, p);
            if i >= list.len() {
                continue;
            }
            let (pid, q) = list[i];
            let m = remin(pid, q);
            let mut g = acquire(&shards[shard_of(m, n)], &mut c);
            g.reserve(1, remin);
            g.insert(m, pid, q);
        }
    });
}

/// Merge index shards into one index over `store`. Entries are listed (8 B
/// each) as each shard is dropped, so shards and result do not coexist in
/// full; their minimizers are recomputed from the pieces by all cores, chunk
/// by chunk, and each chunk is inserted in a scattered order, its head
/// blocks fetched in batches first.
fn merge_shards(
    shards: Box<[Pad<RwLock<MinimizerIndex>>]>,
    store: &PieceStore,
    l: usize,
    cfg: &Config,
) -> MinimizerIndex {
    let mut total = 0u64;
    for s in shards.iter() {
        total += rlock(s).num_entries();
    }
    let mut list: Vec<(PieceId, u32)> = Vec::with_capacity(total as usize);
    let mut hot: Vec<u64> = Vec::new();
    for s in shards.into_vec() {
        let s = s.0.into_inner().unwrap_or_else(PoisonError::into_inner);
        list.extend(s.iter_entries().map(|e| (e.piece_id, e.q)));
        hot.extend(s.hot_minimizers().copied());
    }
    let remin = |p: PieceId, q: u32| order_value(codec::canonicalize(store.kmer_at(p, q, l), l));
    // Sized for the exact entry count at ≈ 80 % head load (as defrag), not
    // grown from 2^bucket_bits by ×1.5 steps (which land anywhere between
    // 57 and 85 % load).
    let buckets = (list.len() as f64 / (0.8 * crate::index::BLOCK_SLOTS as f64)).ceil() as usize;
    let mut out =
        MinimizerIndex::with_buckets(buckets.max(1)).with_spill_dir(cfg.spill_dir.clone());
    out.reserve(list.len() as u64, remin);
    let n = list.len();
    let p = n.next_power_of_two().max(1);
    let threads = conversion_threads(n);
    // Minimizers of one chunk of scattered positions (u64::MAX: skipped).
    let chunk = (CONVERT_CHUNK * threads).min(p);
    let mut ms = vec![0u64; chunk];
    let mut j0 = 0;
    while j0 < p {
        let len = chunk.min(p - j0);
        let sub = len.div_ceil(threads).max(1);
        std::thread::scope(|sc| {
            for (t, part) in ms[..len].chunks_mut(sub).enumerate() {
                let (list, remin) = (&list, &remin);
                let start = j0 + t * sub;
                let mut job = move || {
                    for (o, slot) in part.iter_mut().enumerate() {
                        let i = scattered_at(start + o, p);
                        *slot = if i < n {
                            let (pid, q) = list[i];
                            remin(pid, q)
                        } else {
                            u64::MAX
                        };
                    }
                };
                if threads > 1 {
                    sc.spawn(job);
                } else {
                    job();
                }
            }
        });
        const BATCH: usize = 32;
        let mut sink = 0usize;
        for (b, ms_b) in ms[..len].chunks(BATCH).enumerate() {
            for &m in ms_b {
                if m != u64::MAX {
                    sink = sink.wrapping_add(out.touch(m));
                }
            }
            for (o, &m) in ms_b.iter().enumerate() {
                if m != u64::MAX {
                    let (pid, q) = list[scattered_at(j0 + b * BATCH + o, p)];
                    out.insert(m, pid, q);
                }
            }
        }
        std::hint::black_box(sink);
        j0 += len;
    }
    for m in hot {
        out.mark_hot(m);
    }
    out.mark_baseline();
    out
}

/// Where a writer's previous k-mer was found: piece, global offset of its
/// first base, and whether the piece holds it in the read's orientation.
#[derive(Clone, Copy)]
struct Cursor {
    pid: PieceId,
    g: u64,
    same: bool,
}

/// The lock a writer keeps between consecutive k-mers of the same shard.
enum Held<'a> {
    None,
    L1(usize, RwLockWriteGuard<'a, MinimizerIndex>),
    L2(usize, RwLockWriteGuard<'a, MinimizerIndex>),
}

/// Inserts reads into a [`ConcurrentKmerTable`] from one thread. Same
/// semantics as [`KmerTable::insert_read`] / `insert_chunk` / `end_read`.
pub struct Writer<'a> {
    t: &'a ConcurrentKmerTable,
    held: Held<'a>,
    s: WriterState<'a>,
}

/// Everything a writer owns except its lock.
struct WriterState<'a> {
    store: &'a SharedStore,
    roll: RollingKmer,
    lroll: RollingKmer,
    wmin: WindowMin,
    seg_len: u64,
    in_read: bool,
    run_kmers: usize,
    /// Open piece: id, start, length, and its record.
    open: Option<(PieceId, u64, u32, &'a AtomicU64)>,
    mark_base: u32,
    marks1: u64,
    marks2: u64,
    /// Current chunk: next free base and end, and the segment of the base
    /// array holding it (index of its first word).
    cur: u64,
    end: u64,
    seg: (u64, &'a [AtomicU64]),
    /// Current block of piece ids: next and end.
    pid_next: u64,
    pid_end: u64,
    /// Minimizers seen hot (hotness is permanent), and the last answer.
    known_hot1: FastSet,
    known_hot2: FastSet,
    last1: Option<(u64, bool)>,
    /// Where the previous k-mer was found (see [`Writer::kmer`]).
    prev: Option<Cursor>,
    /// Shard locks taken so far (`Held` changes); the candidates of
    /// `cache_m` were collected under acquisition `cache_acq`.
    acquisitions: u64,
    cache_m: u64,
    cache_acq: u64,
    cache_valid: bool,
    cache: Vec<(PieceId, u32)>,
    stats: InsertStats,
    cont: ContentionStats,
    stored: u64,
    pieces: u64,
}

impl<'a> Writer<'a> {
    fn new(t: &'a ConcurrentKmerTable) -> Self {
        let (k, l) = (t.cfg.k, t.cfg.l);
        Self {
            t,
            held: Held::None,
            s: WriterState {
                store: &t.store,
                roll: RollingKmer::new(k),
                lroll: RollingKmer::new(l),
                wmin: WindowMin::new(k - l + 1),
                seg_len: 0,
                in_read: false,
                run_kmers: 0,
                open: None,
                mark_base: 0,
                marks1: 0,
                marks2: 0,
                cur: 0,
                end: 0,
                seg: (0, &[]),
                pid_next: 0,
                pid_end: 0,
                known_hot1: FastSet::default(),
                known_hot2: FastSet::default(),
                last1: None,
                prev: None,
                acquisitions: 0,
                cache_m: 0,
                cache_acq: 0,
                cache_valid: false,
                cache: Vec::new(),
                stats: InsertStats::default(),
                cont: ContentionStats::default(),
                stored: 0,
                pieces: 0,
            },
        }
    }

    /// Insert all k-mers of `read` (see [`KmerTable::insert_read`]).
    pub fn insert_read(&mut self, read: &[u8]) {
        self.insert_chunk(read);
        self.end_read();
    }

    /// Feed the next chunk of the current read. Every k-mer completed by the
    /// chunk is visible to all threads when this returns.
    pub fn insert_chunk(&mut self, chunk: &[u8]) {
        let s = &mut self.s;
        s.in_read = true;
        let k = self.t.cfg.k;
        let max_run = self.t.cfg.max_piece_len as usize - k + 1;
        for &c in chunk {
            let b = encode(c);
            let s = &mut self.s;
            if b == INVALID {
                s.end_segment(k);
                continue;
            }
            s.seg_len += 1;
            s.roll.push(b);
            if s.lroll.push(b) {
                s.wmin.push(order_value(s.lroll.canonical()));
            }
            if !s.roll.is_full() {
                continue;
            }
            let (m, ties) = s.wmin.current().expect("window complete");
            let x = s.roll.canonical();
            s.stats.kmers_processed += 1;
            let mut split = false;
            if s.run_kmers >= max_run {
                s.end_run();
                split = true;
            }
            self.kmer(x, m, ties, b, split);
        }
        // Do not keep a shard locked while the caller does something else.
        self.held = Held::None;
    }

    /// Finish the current read (see [`KmerTable::end_read`]).
    pub fn end_read(&mut self) {
        let s = &mut self.s;
        s.end_segment(self.t.cfg.k);
        if s.in_read {
            s.stats.reads += 1;
            s.in_read = false;
        }
        self.held = Held::None;
    }

    /// Look up canonical `x` (window minimizer `(m, ties)`, last base `b`)
    /// in its home shard and, if absent, append it and file it there.
    /// Shortcuts (as in [`KmerTable`]'s insertion): if the previous k-mer
    /// was found, this one is present when the piece's current view still
    /// covers the next position and one base matches (no lock, no index
    /// access); and the layer-1 candidates of `m` are reused while this
    /// writer has held `m`'s shard lock since collecting them (nobody else
    /// changed the shard; this writer adds its own entries to them).
    fn kmer(&mut self, x: u64, m: u64, ties: u32, b: u8, split: bool) {
        let t: &'a ConcurrentKmerTable = self.t;
        let (k, l) = (t.cfg.k, t.cfg.l);
        let s = &mut self.s;
        let f = s.roll.fwd();
        let mut trace = QueryTrace::default();

        if let Some(c) = s.prev {
            if let Some(c) = t.extend(&s.local(t), c, b) {
                s.prev = Some(c);
                s.stats.extensions += 1;
                s.present();
                return;
            }
        }

        // Layer 1, unless m is known to be hot.
        if t.l2 == 0 || !s.known_hot1(m) {
            let sh = shard_of(m, t.l1.len());
            if !matches!(&self.held, Held::L1(h, _) if *h == sh) {
                self.held = Held::None;
                self.held = Held::L1(sh, acquire(&t.l1[sh], &mut s.cont));
                s.acquisitions += 1;
            }
            let Held::L1(_, g) = &mut self.held else {
                unreachable!()
            };
            if !(t.l2 > 0 && g.is_hot(m)) {
                if s.cache_valid && s.cache_m == m && s.cache_acq == s.acquisitions {
                    s.stats.cache_reuses += 1;
                } else {
                    s.cache.clear();
                    let mut scan = ScanTrace::default();
                    g.collect_candidates(m, &mut scan, &mut s.cache);
                    trace.add_scan(&scan);
                    s.cache_m = m;
                    s.cache_acq = s.acquisitions;
                    s.cache_valid = true;
                }
                let local = s.local(t);
                let found = t.find_in(&local, &s.cache, x, ties, l, t.w, &mut trace);
                s.stats.lookup.record(found.is_some(), &trace);
                if let Some(h) = found {
                    s.prev = t.cursor(&local, &h, f);
                    s.present();
                } else {
                    s.prev = None;
                    let (pid, s0) = s.append(t, f, b, split);
                    let new = ties as u64 & !s.marks1;
                    if new != 0 {
                        s.marks1 |= new;
                        reserve_timed(g, new.count_ones() as u64, t.remin(l), &mut s.cont);
                        let mut overflowed = false;
                        let mut bits = new;
                        while bits != 0 {
                            let q = s0 + bits.trailing_zeros();
                            overflowed |= g.insert(m, pid, q);
                            if s.cache_m == m {
                                s.cache.push((pid, q));
                            }
                            bits &= bits - 1;
                        }
                        if overflowed && t.l2 > 0 {
                            let t0 = Instant::now();
                            t.promote1(g, m, &mut s.stats);
                            s.cont.promotion_ns += t0.elapsed().as_nanos() as u64;
                            s.cache_valid = false;
                        }
                    }
                }
                return;
            }
            self.held = Held::None;
            s.known_hot1.insert(m);
            s.last1 = Some((m, true));
        }

        // Layer 2, unless its ℓ2-minimizer is known to be hot.
        let l2 = t.l2;
        let (m2, ties2) = kmer_minimizer_mask(f, k, l2);
        trace.layer2 += 1;
        if !s.known_hot2.contains(&m2) {
            let sh = shard_of(m2, t.l2s.len());
            if !matches!(&self.held, Held::L2(h, _) if *h == sh) {
                self.held = Held::None;
                self.held = Held::L2(sh, acquire(&t.l2s[sh], &mut s.cont));
                s.acquisitions += 1;
            }
            let Held::L2(_, g) = &mut self.held else {
                unreachable!()
            };
            if !g.is_hot(m2) {
                let local = s.local(t);
                let found = t.find(&local, g, x, m2, ties2, l2, t.w2, &mut trace);
                s.stats.lookup.record(found.is_some(), &trace);
                if let Some(h) = found {
                    s.prev = t.cursor(&local, &h, f);
                    s.present();
                } else {
                    s.prev = None;
                    let (pid, s0) = s.append(t, f, b, split);
                    let new = ties2 as u64 & !s.marks2;
                    if new != 0 {
                        s.marks2 |= new;
                        reserve_timed(g, new.count_ones() as u64, t.remin(l2), &mut s.cont);
                        let mut overflowed = false;
                        let mut bits = new;
                        while bits != 0 {
                            overflowed |= g.insert(m2, pid, s0 + bits.trailing_zeros());
                            bits &= bits - 1;
                        }
                        if overflowed {
                            let t0 = Instant::now();
                            t.promote2(g, m2, &mut s.stats);
                            s.cont.promotion_ns += t0.elapsed().as_nanos() as u64;
                        }
                    }
                }
                return;
            }
            self.held = Held::None;
            s.known_hot2.insert(m2);
        }

        // Layer 3.
        self.held = Held::None;
        trace.layer2 -= 1;
        trace.layer3 += 1;
        let mut g = acquire(&t.l3[shard_of(x, LAYER3_SHARDS)], &mut s.cont);
        let loc = g.get(x);
        s.stats.lookup.record(loc.is_some(), &trace);
        if let Some((piece, pos)) = loc {
            drop(g);
            let h = Hit {
                piece,
                pos,
                kmer: t.store.kmer_at(piece, pos, k),
            };
            s.prev = t.cursor(&s.local(t), &h, f);
            s.present();
        } else {
            s.prev = None;
            let loc = s.append(t, f, b, split);
            g.add(x, loc);
        }
    }
}

impl<'a> WriterState<'a> {
    fn end_segment(&mut self, k: usize) {
        self.end_run();
        if self.seg_len >= k as u64 {
            self.stats.segments += 1;
        } else if self.seg_len > 0 {
            self.stats.short_segments += 1;
        }
        self.seg_len = 0;
        self.prev = None;
        self.roll.reset();
        self.lroll.reset();
        self.wmin.clear();
    }

    fn end_run(&mut self) {
        if self.run_kmers > 0 {
            if let Some((pid, ..)) = self.open {
                self.store.clear_open(pid);
            }
            self.stored += self.run_kmers as u64;
            self.pieces += 1;
            self.open = None;
            self.run_kmers = 0;
        }
    }

    fn present(&mut self) {
        self.stats.present_hits += 1;
        self.end_run();
    }

    #[inline(always)]
    fn local(&self, t: &'a ConcurrentKmerTable) -> Local<'a> {
        Local {
            store: &t.store,
            open: self.open.map(|(p, s, l, _)| (p, s, l)),
            seg_first: self.seg.0,
            seg: self.seg.1,
        }
    }

    #[inline]
    fn known_hot1(&mut self, m: u64) -> bool {
        match self.last1 {
            Some((lm, h)) if lm == m => h,
            _ => {
                let h = self.known_hot1.contains(&m);
                self.last1 = Some((m, h));
                h
            }
        }
    }

    /// Write base `b` at global offset `g` of the current chunk (only the
    /// chunk's owner writes its words).
    #[inline(always)]
    fn write_code(&self, g: u64, b: u8) {
        let a = &self.seg.1[((g >> 5) - self.seg.0) as usize];
        let w = a.load(Ordering::Relaxed);
        a.store(
            w | ((b as u64 & 3) << (62 - 2 * (g & 31))),
            Ordering::Relaxed,
        );
    }

    /// Append the absent k-mer just read (forward `f`, last base `b`) to the
    /// open piece, or open a new one; returns its `(piece, start)`.
    fn append(&mut self, t: &'a ConcurrentKmerTable, f: u64, b: u8, split: bool) -> (PieceId, u32) {
        let k = t.cfg.k as u64;
        if let Some((pid, start, len, rec)) = self.open {
            if self.cur < self.end {
                self.write_code(self.cur, b);
                self.cur += 1;
                // Release: the base is visible to whoever sees the length.
                rec.store(pack(start, len + 1), Ordering::Release);
                self.open = Some((pid, start, len + 1, rec));
                let s0 = self.run_kmers as u32;
                self.run_kmers += 1;
                let d = s0 - self.mark_base;
                self.marks1 = self.marks1.checked_shr(d).unwrap_or(0);
                self.marks2 = self.marks2.checked_shr(d).unwrap_or(0);
                self.mark_base = s0;
                return (pid, s0);
            }
            // Chunk full: continue the run in a new piece.
            self.end_run();
            self.stats.chunk_splits += 1;
        } else {
            self.stats.cap_splits += split as u64;
        }
        if self.end - self.cur < k {
            let (s, e) = t.store.alloc_chunk(k);
            self.cur = s;
            self.end = e;
            self.seg = t.store.data.segment(s >> 5);
            debug_assert!((e - 1) >> 5 < self.seg.0 + self.seg.1.len() as u64);
        }
        let start = self.cur;
        for i in 0..k {
            self.write_code(start + i, (f >> (2 * (k - 1 - i))) as u8);
        }
        self.cur += k;
        let pid = match t.store.reuse_pid() {
            Some(id) => id,
            None => {
                if self.pid_next == self.pid_end {
                    (self.pid_next, self.pid_end) = t.store.alloc_pids();
                }
                self.pid_next += 1;
                (self.pid_next - 1) as PieceId
            }
        };
        t.store.set_open(pid);
        let rec = t.store.meta.slot(pid as u64);
        rec.store(pack(start, k as u32), Ordering::Release);
        self.open = Some((pid, start, k as u32, rec));
        self.run_kmers = 1;
        self.mark_base = 0;
        self.marks1 = 0;
        self.marks2 = 0;
        (pid, 0)
    }
}

impl Drop for Writer<'_> {
    fn drop(&mut self) {
        self.end_read();
        let (t, s) = (self.t, &self.s);
        if s.end - s.cur >= t.cfg.k as u64 {
            lock(&t.store.tails).push((s.cur, s.end));
        }
        if s.pid_next < s.pid_end {
            lock(&t.store.pid_tails).push((s.pid_next, s.pid_end));
        }
        let mut c = lock(&t.carry);
        c.stats.merge(&s.stats);
        c.contention.merge(&s.cont);
        c.stored_kmers = c.stored_kmers.wrapping_add(s.stored);
        c.pieces_created += s.pieces;
    }
}

#[path = "concurrent_delete.rs"]
mod delete;

#[path = "concurrent_defrag.rs"]
mod defrag;

#[path = "concurrent_stream.rs"]
mod stream;
pub use stream::ConcurrentStreamingQuery;

#[path = "concurrent_snapshot.rs"]
mod snapshot;
pub use snapshot::ReadSnapshot;
