//! k-mer codec.
//!
//! Base encoding is A=0, C=1, G=2, T=3, so the complement of `b` is `3 - b`.
//! A k-mer (k ≤ 32) is packed MSB-first into the low `2k` bits of a `u64`:
//! the first base occupies bits `2(k-1)..2k`.
//!
//! Every shift in this module is bounded by 62 (`2 * (32 - 1)`), so no
//! shift by 64 is ever executed, including for `k = 32`.

/// Sentinel returned by [`encode`] for characters outside `ACGTacgt`.
pub const INVALID: u8 = 4;

/// Largest supported k-mer length (a canonical k-mer must fit in a `u64`).
pub const MAX_K: usize = 32;

static ENCODE_TABLE: [u8; 256] = {
    let mut t = [INVALID; 256];
    t[b'A' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b'a' as usize] = 0;
    t[b'c' as usize] = 1;
    t[b'g' as usize] = 2;
    t[b't' as usize] = 3;
    t
};

/// Uppercase ASCII letter for each 2-bit code.
pub const DECODE: [u8; 4] = *b"ACGT";

/// Encode one character. Lowercase (soft-masked) bases are accepted.
/// Returns [`INVALID`] for anything that is not `ACGT` in either case.
#[inline(always)]
pub fn encode(c: u8) -> u8 {
    ENCODE_TABLE[c as usize]
}

/// Decode a 2-bit code (only the low two bits are used).
#[inline(always)]
pub fn decode(b: u8) -> u8 {
    DECODE[(b & 3) as usize]
}

/// Mask with the low `2k` bits set. Valid for `1 ≤ k ≤ 32`.
#[inline(always)]
pub fn kmer_mask(k: usize) -> u64 {
    debug_assert!((1..=MAX_K).contains(&k));
    if k == MAX_K {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    }
}

/// Canonical form: the numerically smaller of a k-mer and its reverse complement.
#[inline(always)]
pub fn canonical(fwd: u64, rc: u64) -> u64 {
    fwd.min(rc)
}

/// Reverse complement of a packed k-mer.
#[inline]
pub fn revcomp_kmer(x: u64, k: usize) -> u64 {
    debug_assert!((1..=MAX_K).contains(&k));
    // Complement every base (3 - b == b ^ 3), then reverse the 2-bit groups
    // of the whole word and shift the k-mer back down to the low bits.
    let mut y = !x;
    y = ((y >> 2) & 0x3333_3333_3333_3333) | ((y & 0x3333_3333_3333_3333) << 2);
    y = ((y >> 4) & 0x0F0F_0F0F_0F0F_0F0F) | ((y & 0x0F0F_0F0F_0F0F_0F0F) << 4);
    y = y.swap_bytes();
    // 64 - 2k is in [0, 62] for 1 ≤ k ≤ 32.
    y >> (64 - 2 * k)
}

/// Canonical form of a packed k-mer.
#[inline]
pub fn canonicalize(x: u64, k: usize) -> u64 {
    canonical(x, revcomp_kmer(x, k))
}

/// Encode a k-mer given as ASCII directly (no rolling). Returns `None` if
/// the slice is longer than 32 or contains a non-ACGT character.
#[inline]
pub fn encode_kmer(seq: &[u8]) -> Option<u64> {
    if seq.is_empty() || seq.len() > MAX_K {
        return None;
    }
    let mut x = 0u64;
    for &c in seq {
        let b = encode(c);
        if b == INVALID {
            return None;
        }
        x = (x << 2) | b as u64;
    }
    Some(x)
}

/// Encode the forward and reverse-complement k-mers of a slice known to be
/// valid ACGT (as piece storage guarantees). Allocation-free.
#[inline(always)]
pub fn encode_kmer_both(seq: &[u8]) -> (u64, u64) {
    let k = seq.len();
    debug_assert!((1..=MAX_K).contains(&k));
    let shift = 2 * (k - 1);
    let mut fwd = 0u64;
    let mut rc = 0u64;
    for &c in seq {
        let b = encode(c) as u64;
        debug_assert!(b < 4, "invalid base in stored sequence");
        fwd = (fwd << 2) | b;
        rc = (rc >> 2) | ((3 - b) << shift);
    }
    (fwd, rc)
}

/// Decode a packed k-mer to an uppercase ASCII string.
pub fn decode_kmer(x: u64, k: usize) -> String {
    (0..k)
        .map(|i| decode((x >> (2 * (k - 1 - i))) as u8) as char)
        .collect()
}

/// Reverse complement of an ASCII sequence (output is uppercase ACGT).
/// Non-ACGT characters are mapped to `N`.
pub fn revcomp_seq(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&c| match encode(c) {
            INVALID => b'N',
            b => decode(3 - b),
        })
        .collect()
}

/// Rolling forward/reverse-complement k-mer over a stream of bases.
///
/// ```text
/// fwd = ((fwd << 2) | b) & mask_k
/// rc  = (rc >> 2) | ((3 - b) << 2(k-1))
/// ```
#[derive(Clone, Debug)]
pub struct RollingKmer {
    k: usize,
    mask: u64,
    shift: u32,
    fwd: u64,
    rc: u64,
    filled: usize,
}

impl RollingKmer {
    pub fn new(k: usize) -> Self {
        assert!((1..=MAX_K).contains(&k), "k must be in 1..=32, got {k}");
        Self {
            k,
            mask: kmer_mask(k),
            shift: 2 * (k as u32 - 1),
            fwd: 0,
            rc: 0,
            filled: 0,
        }
    }

    /// Forget all bases (e.g. after a non-ACGT character).
    #[inline]
    pub fn reset(&mut self) {
        self.fwd = 0;
        self.rc = 0;
        self.filled = 0;
    }

    /// Push one 2-bit code (must be < 4). Returns true once k bases are loaded.
    #[inline(always)]
    pub fn push(&mut self, b: u8) -> bool {
        debug_assert!(b < 4);
        let b = b as u64;
        self.fwd = ((self.fwd << 2) | b) & self.mask;
        self.rc = (self.rc >> 2) | ((3 - b) << self.shift);
        if self.filled < self.k {
            self.filled += 1;
        }
        self.filled == self.k
    }

    #[inline(always)]
    pub fn is_full(&self) -> bool {
        self.filled == self.k
    }
    #[inline(always)]
    pub fn fwd(&self) -> u64 {
        self.fwd
    }
    #[inline(always)]
    pub fn rc(&self) -> u64 {
        self.rc
    }
    #[inline(always)]
    pub fn canonical(&self) -> u64 {
        canonical(self.fwd, self.rc)
    }
}

/// Canonical k-mers of a read delivered in arbitrary chunks. k-mers spanning
/// chunk boundaries are reported; non-ACGT characters reset the window.
#[derive(Clone, Debug)]
pub struct KmerStream {
    roll: RollingKmer,
}

impl KmerStream {
    pub fn new(k: usize) -> Self {
        Self {
            roll: RollingKmer::new(k),
        }
    }

    /// Feed the next chunk of the current read.
    #[inline]
    pub fn push<F: FnMut(u64)>(&mut self, chunk: &[u8], mut f: F) {
        for &c in chunk {
            let b = encode(c);
            if b == INVALID {
                self.roll.reset();
            } else if self.roll.push(b) {
                f(self.roll.canonical());
            }
        }
    }

    /// End the current read.
    #[inline]
    pub fn end(&mut self) {
        self.roll.reset();
    }
}

/// Iterate over the canonical k-mers of an arbitrary read, skipping every
/// k-mer that overlaps a non-ACGT character. Used by the oracle and tests.
pub fn canonical_kmers(read: &[u8], k: usize) -> impl Iterator<Item = u64> + '_ {
    let mut roll = RollingKmer::new(k);
    read.iter().filter_map(move |&c| {
        let b = encode(c);
        if b == INVALID {
            roll.reset();
            None
        } else if roll.push(b) {
            Some(roll.canonical())
        } else {
            None
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mask_bounds() {
        assert_eq!(kmer_mask(1), 0b11);
        assert_eq!(kmer_mask(31), (1u64 << 62) - 1);
        assert_eq!(kmer_mask(32), u64::MAX);
    }

    #[test]
    fn revcomp_small() {
        let x = encode_kmer(b"AACGT").unwrap();
        assert_eq!(decode_kmer(revcomp_kmer(x, 5), 5), "ACGTT");
        let x = encode_kmer(b"A").unwrap();
        assert_eq!(decode_kmer(revcomp_kmer(x, 1), 1), "T");
    }
}
