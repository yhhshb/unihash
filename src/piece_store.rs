//! M2 + M4 — packed piece storage.
//!
//! All bases live back to back in one array of 2-bit bases (32 bases per
//! `u64`, first base in the most significant bits). A piece is a *view*
//! `(start, len)` into that array, packed into one `u64` (40-bit start,
//! 24-bit length). Pieces are not word-aligned, so there is no padding
//! between them. One zero word is kept past the end so that any k-mer
//! (k ≤ 32) can be read with two word loads and a shift, without branching
//! on word boundaries.
//!
//! Views make deletion cheap: splitting a piece only shortens it and adds a
//! second view on the same bases, without copying. Bases no longer covered
//! by any piece are *dead* until [`defragmentation`](crate::KmerTable::defrag)
//! rewrites the array.
//!
//! Cost: 2 bits per base plus 8 bytes per piece.

use crate::codec::{self, encode, INVALID};

pub type PieceId = u32;

/// Longest piece (24-bit length field).
pub const MAX_PIECE_LEN: u32 = (1 << 24) - 1;
/// Largest base array (40-bit start field).
pub const MAX_STORE_BASES: u64 = 1 << 40;

const LEN_BITS: u32 = 24;
const LEN_MASK: u64 = (1 << LEN_BITS) - 1;

#[inline(always)]
fn pack(start: u64, len: u32) -> u64 {
    debug_assert!(start < MAX_STORE_BASES && len <= MAX_PIECE_LEN);
    (start << LEN_BITS) | len as u64
}

#[derive(Clone, Debug)]
pub struct PieceStore {
    /// Packed bases; always `ceil(n_bases / 32) + 1` words (one pad word).
    data: Vec<u64>,
    /// Bases written to `data` (live and dead).
    n_bases: u64,
    /// Piece `i` is `meta[i] = start << 24 | len`.
    meta: Vec<u64>,
    /// Σ piece lengths (bases shared by overlapping views count twice).
    live_bases: u64,
    /// Pieces of length 0 (fully deleted or merged away).
    empty: u64,
    /// Ids of emptied pieces, reused by the next pieces created.
    free: Vec<PieceId>,
}

impl Default for PieceStore {
    fn default() -> Self {
        Self::new()
    }
}

impl PieceStore {
    pub fn new() -> Self {
        Self {
            data: vec![0],
            n_bases: 0,
            meta: Vec::new(),
            live_bases: 0,
            empty: 0,
            free: Vec::new(),
        }
    }

    /// An empty store whose base array can hold `bases` without growing.
    pub fn with_capacity(bases: u64, pieces: usize) -> Self {
        let mut s = Self::new();
        s.data.reserve_exact((bases as usize).div_ceil(32) + 1);
        s.meta.reserve_exact(pieces);
        s
    }

    fn grow_data(&mut self, new_len: u64) {
        assert!(
            new_len <= MAX_STORE_BASES,
            "piece store full: at most 2^40 bases"
        );
        let words = (new_len as usize).div_ceil(32) + 1;
        // Grow by ~12.5 % instead of doubling, to bound slack (peak memory).
        if words > self.data.capacity() {
            let extra = (words - self.data.len()).max(self.data.len() / 8);
            self.data.reserve_exact(extra);
        }
        self.data.resize(words, 0);
    }

    fn push_meta(&mut self, start: u64, len: u32) -> PieceId {
        if len > 0 {
            if let Some(id) = self.free.pop() {
                debug_assert_eq!(self.length(id), 0);
                self.meta[id as usize] = pack(start, len);
                self.live_bases += len as u64;
                self.empty -= 1;
                return id;
            }
        }
        let id = self.meta.len();
        assert!(
            id < u32::MAX as usize,
            "piece store full: piece ids are u32 (max {} pieces)",
            u32::MAX
        );
        if self.meta.len() == self.meta.capacity() {
            self.meta.reserve_exact((self.meta.len() / 8).max(16));
        }
        self.meta.push(pack(start, len));
        self.live_bases += len as u64;
        self.empty += (len == 0) as u64;
        id as PieceId
    }

    /// Append a piece (ASCII `ACGTacgt`) and return its id.
    ///
    /// # Panics
    /// If `bases` contains a non-ACGT character, is longer than
    /// [`MAX_PIECE_LEN`], or if more than `u32::MAX` pieces would exist.
    pub fn add_piece(&mut self, bases: &[u8]) -> PieceId {
        assert!(
            bases.len() <= MAX_PIECE_LEN as usize,
            "piece of {} bases exceeds the limit of {MAX_PIECE_LEN}",
            bases.len()
        );
        let start = self.n_bases;
        self.grow_data(start + bases.len() as u64);
        for (g, &c) in (start as usize..).zip(bases) {
            let b = encode(c);
            assert!(
                b != INVALID,
                "add_piece: non-ACGT character {:?}",
                c as char
            );
            self.data[g >> 5] |= (b as u64) << (62 - 2 * (g & 31));
        }
        self.n_bases = start + bases.len() as u64;
        self.push_meta(start, bases.len() as u32)
    }

    /// Append one base (2-bit code) to the array without creating a piece;
    /// see [`close_piece`](Self::close_piece).
    #[inline]
    pub fn push_code(&mut self, b: u8) {
        let g = self.n_bases;
        if (g as usize) >> 5 >= self.data.len() - 1 {
            self.grow_data(g + 1);
        }
        self.data[(g >> 5) as usize] |= (b as u64 & 3) << (62 - 2 * (g & 31));
        self.n_bases = g + 1;
    }

    /// Make the bases pushed since global offset `start` a piece.
    pub fn close_piece(&mut self, start: u64) -> PieceId {
        let len = self.n_bases - start;
        assert!(len <= MAX_PIECE_LEN as u64, "piece too long");
        self.push_meta(start, len as u32)
    }

    /// Bases written so far (the next base goes to this global offset).
    #[inline]
    pub fn end(&self) -> u64 {
        self.n_bases
    }

    /// Add a piece viewing existing bases `start..start + len`.
    pub fn add_view(&mut self, start: u64, len: u32) -> PieceId {
        assert!(start + len as u64 <= self.n_bases);
        self.push_meta(start, len)
    }

    /// Re-point piece `id` to `start..start + len` (bases must exist).
    /// Setting the length to 0 frees the id for reuse; a freed id cannot be
    /// re-pointed (it is handed out again by the next piece creation).
    pub fn set_view(&mut self, id: PieceId, start: u64, len: u32) {
        assert!(start + len as u64 <= self.n_bases);
        let old = self.length(id);
        assert!(old > 0 || len == 0, "set_view on an emptied piece");
        if old > 0 && len == 0 {
            self.free.push(id);
        }
        self.live_bases = self.live_bases - old as u64 + len as u64;
        self.empty = self.empty - (old == 0) as u64 + (len == 0) as u64;
        self.meta[id as usize] = pack(start, len);
    }

    /// Number of pieces (including emptied ones).
    #[inline]
    pub fn len(&self) -> usize {
        self.meta.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Pieces of length 0 (deleted entirely; ids are reclaimed by
    /// defragmentation).
    pub fn empty_pieces(&self) -> u64 {
        self.empty
    }

    /// Bases in the array, live or dead.
    #[inline]
    pub fn total_bases(&self) -> u64 {
        self.n_bases
    }

    /// Σ piece lengths.
    #[inline]
    pub fn live_bases(&self) -> u64 {
        self.live_bases
    }

    /// Global offset of the first base of piece `id`.
    #[inline(always)]
    pub fn start(&self, id: PieceId) -> u64 {
        self.meta[id as usize] >> LEN_BITS
    }

    /// Length (in bases) of piece `id`.
    #[inline(always)]
    pub fn length(&self, id: PieceId) -> u32 {
        (self.meta[id as usize] & LEN_MASK) as u32
    }

    /// `(start, len)` of piece `id` with one load.
    #[inline(always)]
    pub fn view(&self, id: PieceId) -> (u64, u32) {
        let m = self.meta[id as usize];
        (m >> LEN_BITS, (m & LEN_MASK) as u32)
    }

    /// 2-bit code of the base at global offset `g`.
    #[inline(always)]
    pub fn code_at_global(&self, g: u64) -> u8 {
        ((self.data[(g >> 5) as usize] >> (62 - 2 * (g & 31))) & 3) as u8
    }

    /// Forward packed k-mer (`k ≤ 32`) starting at global offset `g`.
    #[inline(always)]
    pub fn kmer_at_global(&self, g: u64, k: usize) -> u64 {
        debug_assert!((1..=32).contains(&k));
        let i = (g >> 5) as usize;
        let r = (g & 31) as u32;
        let v = ((self.data[i] as u128) << 64) | self.data[i + 1] as u128;
        // 2r ≤ 62 and 128 - 2k ∈ [64, 126]: no shift reaches 128.
        ((v << (2 * r)) >> (128 - 2 * k as u32)) as u64
    }

    /// 2-bit code of base `i` of piece `id`.
    #[inline]
    pub fn base_code(&self, id: PieceId, i: u32) -> u8 {
        self.code_at_global(self.start(id) + i as u64)
    }

    /// Base character (uppercase) at offset `i` of piece `id`.
    #[inline]
    pub fn get_base(&self, id: PieceId, i: u32) -> u8 {
        codec::decode(self.base_code(id, i))
    }

    /// Forward packed k-mer starting at offset `s` of piece `id`.
    #[inline]
    pub fn kmer_at(&self, id: PieceId, s: u32, k: usize) -> u64 {
        debug_assert!(s as usize + k <= self.length(id) as usize);
        self.kmer_at_global(self.start(id) + s as u64, k)
    }

    /// Canonical packed k-mer starting at offset `s` of piece `id`.
    #[inline]
    pub fn canonical_kmer_at(&self, id: PieceId, s: u32, k: usize) -> u64 {
        codec::canonicalize(self.kmer_at(id, s, k), k)
    }

    /// Decoded bases of piece `id` (uppercase ACGT). Allocates; meant for
    /// tests and diagnostics.
    pub fn piece_seq(&self, id: PieceId) -> Vec<u8> {
        let s = self.start(id);
        (0..self.length(id) as u64)
            .map(|i| codec::decode(self.code_at_global(s + i)))
            .collect()
    }

    /// Iterate over all pieces, decoded (emptied pieces yield empty vectors).
    pub fn iter(&self) -> impl Iterator<Item = Vec<u8>> + '_ {
        (0..self.len()).map(move |i| self.piece_seq(i as PieceId))
    }

    /// Heap bytes in use (by length) and reserved (by capacity).
    pub fn heap_bytes(&self) -> (usize, usize) {
        let used = (self.data.len() + self.meta.len()) * 8 + self.free.len() * 4;
        let reserved = (self.data.capacity() + self.meta.capacity()) * 8 + self.free.capacity() * 4;
        (used, reserved)
    }

    pub fn shrink_to_fit(&mut self) {
        self.data.shrink_to_fit();
        self.meta.shrink_to_fit();
        self.free.shrink_to_fit();
    }

    /// Set the base at global offset `g` (other bases untouched).
    #[inline(always)]
    fn set_code_at_global(&mut self, g: u64, b: u8) {
        let sh = 62 - 2 * (g & 31);
        let w = &mut self.data[(g >> 5) as usize];
        *w = (*w & !(3u64 << sh)) | ((b as u64 & 3) << sh);
    }

    /// Rewrite the base array in place without the bases no piece covers
    /// (pieces keep their ids and lengths; only their starts move). Bases
    /// shared by overlapping pieces stay shared. O(bases + pieces log
    /// pieces), no index work. Returns `(bases removed, shared bases)`.
    pub fn compact_bases(&mut self) -> (u64, u64) {
        let mut order: Vec<(u64, u32, PieceId)> = (0..self.meta.len() as PieceId)
            .filter_map(|id| {
                let (s, l) = self.view(id);
                (l > 0).then_some((s, l, id))
            })
            .collect();
        order.sort_unstable();
        // Merge overlapping views into segments, sliding each to the write
        // cursor. The cursor never passes the read position, and bases are
        // written one at a time, so the copy is safe in place.
        let mut w = 0u64;
        let (mut seg_lo, mut seg_hi, mut seg_new) = (0u64, 0u64, 0u64);
        let mut have = false;
        let mut shared = 0u64;
        for (s, len, id) in order {
            let e = s + len as u64;
            let ns = if have && s < seg_hi {
                shared += seg_hi.min(e) - s;
                if e > seg_hi {
                    self.move_bases(seg_hi, w, e - seg_hi);
                    w += e - seg_hi;
                    seg_hi = e;
                }
                seg_new + (s - seg_lo)
            } else {
                (seg_lo, seg_hi, seg_new, have) = (s, e, w, true);
                self.move_bases(s, w, len as u64);
                w += len as u64;
                seg_new
            };
            self.meta[id as usize] = pack(ns, len);
        }
        for m in &mut self.meta {
            if *m & LEN_MASK == 0 {
                *m = 0;
            }
        }
        let removed = self.n_bases - w;
        // Zero everything past the new end (appends OR bases into words).
        let words = (w as usize).div_ceil(32) + 1;
        self.data.truncate(words);
        let last = (w >> 5) as usize;
        let r = (w & 31) as u32;
        self.data[last] = if r == 0 {
            0
        } else {
            self.data[last] & !(u64::MAX >> (2 * r))
        };
        for x in &mut self.data[last + 1..] {
            *x = 0;
        }
        self.data.shrink_to_fit();
        self.n_bases = w;
        (removed, shared)
    }

    /// Move `n` bases from global offset `from` down to `to` (< `from`),
    /// 32 bases per step: read with two word loads, written with masked
    /// two-word updates that touch only the destination bases. Safe in place
    /// because every write lands at or before bases already read.
    fn move_bases(&mut self, from: u64, to: u64, n: u64) {
        if from == to || n == 0 {
            return;
        }
        debug_assert!(to < from);
        let mut i = 0u64;
        while i + 32 <= n {
            let x = self.kmer_at_global(from + i, 32);
            self.write32(to + i, x);
            i += 32;
        }
        for j in i..n {
            let b = self.code_at_global(from + j);
            self.set_code_at_global(to + j, b);
        }
    }

    /// Write 32 bases (`x`, first base in the top bits) at global offset `g`.
    #[inline(always)]
    fn write32(&mut self, g: u64, x: u64) {
        let w = (g >> 5) as usize;
        let r = 2 * (g & 31) as u32;
        if r == 0 {
            self.data[w] = x;
        } else {
            let hi_mask = u64::MAX >> r; // low 64 - r bits of word w
            self.data[w] = (self.data[w] & !hi_mask) | (x >> r);
            let lo_mask = u64::MAX << (64 - r); // high r bits of word w + 1
            self.data[w + 1] = (self.data[w + 1] & !lo_mask) | (x << (64 - r));
        }
    }
}
