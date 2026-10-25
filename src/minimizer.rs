//! M3 — canonical minimizers with tied positions.
//!
//! For a window (a k-mer starting at `s`), its ℓ-mers start at `s..=s+k-ℓ`.
//! The window minimizer value `m` is the minimum *order value* over the
//! canonical ℓ-mers of the window. Every position whose order value equals
//! `m` is a minimizer position (ties are kept, not broken). Because the order
//! value is a function of the canonical ℓ-mer only, `m(x) == m(revcomp(x))`
//! and the tie set of `revcomp(x)` is the mirror image of that of `x`.
//!
//! The order value of a canonical ℓ-mer is a fixed bijective 64-bit mix of
//! it ([`order_value`]). A bijection keeps exactly the tie structure of the
//! plain value (two positions tie iff their canonical ℓ-mers are equal), so
//! every correctness argument of the plan carries over; it only removes the
//! skew of the lexicographic order (poly-A and other low-complexity ℓ-mers
//! being the minimum of every window that contains them). The plain
//! lexicographic order was removed after the chr21 benchmark (4.8× slower
//! insertion, 2.5× slower queries).

use crate::codec::{encode, RollingKmer};
use std::collections::VecDeque;

/// SplitMix64 finalizer: a bijection on `u64`.
#[inline(always)]
pub fn mix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// Inverse of [`mix64`].
#[inline]
pub fn unmix64(mut z: u64) -> u64 {
    z = z ^ (z >> 31) ^ (z >> 62);
    z = z.wrapping_mul(0x3196_42b2_d24d_8ec3);
    z = z ^ (z >> 27) ^ (z >> 54);
    z = z.wrapping_mul(0x96de_1b17_3f11_9089);
    z ^ (z >> 30) ^ (z >> 60)
}

const ORDER_SEED: u64 = 0x9e37_79b9_7f4a_7c15;

/// Order value of a canonical ℓ-mer (injective). The XOR keeps the all-A
/// ℓ-mer (0) from mapping to 0.
#[inline(always)]
pub fn order_value(canonical_lmer: u64) -> u64 {
    mix64(canonical_lmer ^ ORDER_SEED)
}

/// Canonical ℓ-mer of an order value (inverse of [`order_value`]); used to
/// print hot minimizers.
#[inline]
pub fn lmer_of_order_value(m: u64) -> u64 {
    unmix64(m) ^ ORDER_SEED
}

/// Order values of the canonical ℓ-mer at every position `p` of `seq`
/// (`seq.len() - l + 1` values; empty if `seq` is shorter than `l`).
/// `seq` must be ACGT only (case-insensitive).
pub fn lmer_values(seq: &[u8], l: usize) -> Vec<u64> {
    let mut out = Vec::with_capacity(seq.len().saturating_sub(l - 1));
    lmer_values_into(seq, l, &mut out);
    out
}

/// Like [`lmer_values`] but reuses `out` (cleared first).
pub fn lmer_values_into(seq: &[u8], l: usize, out: &mut Vec<u64>) {
    out.clear();
    let mut roll = RollingKmer::new(l);
    for &c in seq {
        let b = encode(c);
        assert!(b < 4, "lmer_values: non-ACGT character {:?}", c as char);
        if roll.push(b) {
            out.push(order_value(roll.canonical()));
        }
    }
}

/// Window minima over windows of `w` consecutive values, via a monotone
/// deque. Output has `vals.len() - w + 1` entries (empty if too short).
pub fn window_minima(vals: &[u64], w: usize) -> Vec<u64> {
    let mut out = Vec::with_capacity(vals.len().saturating_sub(w - 1));
    window_minima_into(vals, w, &mut out);
    out
}

/// Like [`window_minima`] but reuses `out` (cleared first).
pub fn window_minima_into(vals: &[u64], w: usize, out: &mut Vec<u64>) {
    assert!(w >= 1);
    out.clear();
    if vals.len() < w {
        return;
    }
    let mut dq: VecDeque<usize> = VecDeque::with_capacity(w + 1);
    for (i, &v) in vals.iter().enumerate() {
        // Keep equal values (pop strictly greater only) so ties survive.
        while let Some(&b) = dq.back() {
            if vals[b] > v {
                dq.pop_back();
            } else {
                break;
            }
        }
        dq.push_back(i);
        if i + 1 >= w {
            let s = i + 1 - w;
            while dq[0] < s {
                dq.pop_front();
            }
            out.push(vals[dq[0]]);
        }
    }
}

/// Brute-force reference for [`window_minima`].
pub fn window_minima_brute(vals: &[u64], w: usize) -> Vec<u64> {
    if vals.len() < w {
        return Vec::new();
    }
    vals.windows(w)
        .map(|win| *win.iter().min().unwrap())
        .collect()
}

/// All positions `p` in `[s, s + w - 1]` whose value equals the window
/// minimum (brute force; `w = k - ℓ + 1`).
pub fn minimizer_positions_in_window(vals: &[u64], s: usize, w: usize) -> (u64, Vec<usize>) {
    let win = &vals[s..s + w];
    let m = *win.iter().min().expect("empty window");
    let pos = (s..s + w).filter(|&p| vals[p] == m).collect();
    (m, pos)
}

/// Every position of `vals` that is a (tied) minimizer position of at least
/// one window of `w` values. Sorted and deduplicated. The minimizer value of
/// position `q` is `vals[q]`, so these positions are exactly the distinct
/// `(m, q)` pairs to index for a piece (spec §5.3).
pub fn piece_minimizer_positions(vals: &[u64], w: usize) -> Vec<u32> {
    let mut out = Vec::new();
    piece_minimizer_positions_into(vals, w, &mut Vec::new(), &mut out);
    out
}

/// Like [`piece_minimizer_positions`] with caller-provided scratch buffers.
pub fn piece_minimizer_positions_into(
    vals: &[u64],
    w: usize,
    marks: &mut Vec<bool>,
    out: &mut Vec<u32>,
) {
    assert!(w >= 1);
    out.clear();
    if vals.len() < w {
        return;
    }
    marks.clear();
    marks.resize(vals.len(), false);
    let mut dq: VecDeque<usize> = VecDeque::with_capacity(w + 1);
    for (i, &v) in vals.iter().enumerate() {
        while let Some(&b) = dq.back() {
            if vals[b] > v {
                dq.pop_back();
            } else {
                break;
            }
        }
        dq.push_back(i);
        if i + 1 >= w {
            let s = i + 1 - w;
            while dq[0] < s {
                dq.pop_front();
            }
            // The deque is non-decreasing and retains every position holding
            // the window minimum, so the ties are exactly its equal prefix.
            let m = vals[dq[0]];
            for &p in dq.iter() {
                if vals[p] != m {
                    break;
                }
                marks[p] = true;
            }
        }
    }
    out.extend(
        marks
            .iter()
            .enumerate()
            .filter(|(_, &b)| b)
            .map(|(p, _)| p as u32),
    );
}

/// Minimizer value and tied positions (relative to the k-mer start) of a
/// packed k-mer `x`. Used by queries (only `m` is needed) and by tests.
pub fn kmer_minimizer_positions(x: u64, k: usize, l: usize) -> (u64, Vec<usize>) {
    let mut vals = Vec::with_capacity(k - l + 1);
    let mut roll = RollingKmer::new(l);
    for i in 0..k {
        let b = ((x >> (2 * (k - 1 - i))) & 3) as u8;
        if roll.push(b) {
            vals.push(order_value(roll.canonical()));
        }
    }
    minimizer_positions_in_window(&vals, 0, k - l + 1)
}

/// Window minimizer value `m(x)` of a packed k-mer together with its tied
/// positions `Q(x)` as a bit mask (bit `p` set ⇔ position `p` ties;
/// `k - ℓ + 1 ≤ 32` so the mask fits in a `u32`). Allocation-free, O(k).
#[inline]
pub fn kmer_minimizer_mask(x: u64, k: usize, l: usize) -> (u64, u32) {
    let mut roll = RollingKmer::new(l);
    let mut m = u64::MAX;
    let mut mask = 0u32;
    for i in 0..k {
        let b = ((x >> (2 * (k - 1 - i))) & 3) as u8;
        if roll.push(b) {
            let p = i + 1 - l;
            let v = order_value(roll.canonical());
            if v < m {
                m = v;
                mask = 1 << p;
            } else if v == m {
                mask |= 1 << p;
            }
        }
    }
    (m, mask)
}

/// Tie mask of window `s` given its minimum `m` (bit `p` ⇔ `vals[s+p] == m`).
#[inline]
pub fn window_tie_mask(vals: &[u64], s: usize, w: usize, m: u64) -> u32 {
    let mut mask = 0u32;
    for (p, &v) in vals[s..s + w].iter().enumerate() {
        mask |= ((v == m) as u32) << p;
    }
    mask
}

/// Window minimizer value `m(x)` of a packed k-mer. Allocation-free, O(k).
#[inline]
pub fn kmer_minimizer(x: u64, k: usize, l: usize) -> u64 {
    let mut roll = RollingKmer::new(l);
    let mut m = u64::MAX;
    for i in 0..k {
        let b = ((x >> (2 * (k - 1 - i))) & 3) as u8;
        if roll.push(b) {
            m = m.min(order_value(roll.canonical()));
        }
    }
    m
}

/// Streaming window minimum with ties, for ℓ-mer order values fed one at a
/// time. Memory is O(w) regardless of the sequence length.
#[derive(Clone, Debug)]
pub struct WindowMin {
    w: usize,
    /// (ℓ-mer position, value); values non-decreasing, ties retained.
    dq: VecDeque<(u64, u64)>,
    next: u64,
}

impl WindowMin {
    /// `w = k - ℓ + 1` (1 ≤ w ≤ 32).
    pub fn new(w: usize) -> Self {
        assert!((1..=32).contains(&w), "window of {w} ℓ-mers not supported");
        Self {
            w,
            dq: VecDeque::with_capacity(w + 1),
            next: 0,
        }
    }

    /// Start a new sequence.
    #[inline]
    pub fn clear(&mut self) {
        self.dq.clear();
        self.next = 0;
    }

    /// Feed the order value of the next ℓ-mer.
    #[inline]
    pub fn push(&mut self, v: u64) {
        while let Some(&(_, b)) = self.dq.back() {
            if b > v {
                self.dq.pop_back();
            } else {
                break;
            }
        }
        self.dq.push_back((self.next, v));
        self.next += 1;
        if self.next >= self.w as u64 {
            let s = self.next - self.w as u64;
            while self.dq[0].0 < s {
                self.dq.pop_front();
            }
        }
    }

    /// `(m, ties)` of the window made of the last `w` ℓ-mers, where bit `p`
    /// of `ties` is set iff the ℓ-mer at window offset `p` has value `m`.
    /// `None` until `w` values have been pushed.
    #[inline]
    pub fn current(&self) -> Option<(u64, u32)> {
        if self.next < self.w as u64 {
            return None;
        }
        let s = self.next - self.w as u64;
        let m = self.dq[0].1;
        let mut mask = 0u32;
        for &(p, v) in &self.dq {
            if v != m {
                break;
            }
            mask |= 1 << (p - s);
        }
        Some((m, mask))
    }
}
