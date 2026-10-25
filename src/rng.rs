//! Tiny deterministic PRNG and synthetic sequence helpers for tests and
//! benchmarks (no external dependencies, reproducible from a logged seed).

use crate::codec::DECODE;

/// SplitMix64.
#[derive(Clone, Debug)]
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Self(seed)
    }

    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        crate::minimizer::mix64(self.0)
    }

    /// Uniform in `0..n` (`n > 0`), via 128-bit multiply.
    #[inline]
    pub fn below(&mut self, n: u64) -> u64 {
        assert!(n > 0);
        ((self.next_u64() as u128 * n as u128) >> 64) as u64
    }

    /// Uniform in `lo..=hi`.
    pub fn range(&mut self, lo: usize, hi: usize) -> usize {
        lo + self.below((hi - lo + 1) as u64) as usize
    }

    pub fn chance(&mut self, p: f64) -> bool {
        ((self.next_u64() >> 11) as f64 / (1u64 << 53) as f64) < p
    }

    pub fn base(&mut self) -> u8 {
        DECODE[self.below(4) as usize]
    }

    /// Uniform random DNA of length `n`.
    pub fn dna(&mut self, n: usize) -> Vec<u8> {
        (0..n).map(|_| self.base()).collect()
    }

    /// DNA over a reduced alphabet (e.g. `b"AC"`) to force low complexity.
    pub fn dna_from(&mut self, n: usize, alphabet: &[u8]) -> Vec<u8> {
        (0..n)
            .map(|_| alphabet[self.below(alphabet.len() as u64) as usize])
            .collect()
    }

    /// Shuffle in place (Fisher–Yates).
    pub fn shuffle<T>(&mut self, v: &mut [T]) {
        for i in (1..v.len()).rev() {
            let j = self.below(i as u64 + 1) as usize;
            v.swap(i, j);
        }
    }
}

/// A random genome with planted repeats: copies of earlier segments
/// (forward and reverse complement), tandem repeats and homopolymers.
pub fn repeat_rich_genome(rng: &mut Rng, len: usize) -> Vec<u8> {
    let mut g = Vec::with_capacity(len + 1024);
    while g.len() < len {
        match rng.below(10) {
            0..=4 => {
                let n = rng.range(20, 400);
                g.extend(rng.dna(n));
            }
            5 | 6 if g.len() > 100 => {
                // copy an earlier segment, possibly reverse-complemented
                let n = rng.range(20, 300).min(g.len());
                let s = rng.below((g.len() - n + 1) as u64) as usize;
                let seg: Vec<u8> = g[s..s + n].to_vec();
                if rng.chance(0.5) {
                    g.extend(crate::codec::revcomp_seq(&seg));
                } else {
                    g.extend(seg);
                }
            }
            7 => {
                // tandem repeat
                let ulen = rng.range(1, 6);
                let unit = rng.dna(ulen);
                let copies = rng.range(5, 40);
                for _ in 0..copies {
                    g.extend_from_slice(&unit);
                }
            }
            8 => {
                let b = rng.base();
                let n = rng.range(5, 60);
                g.extend(std::iter::repeat(b).take(n));
            }
            _ => {
                let n = rng.range(10, 100);
                g.extend(rng.dna_from(n, b"AC"));
            }
        }
    }
    g.truncate(len);
    g
}

/// Sample `n` reads of length in `[min_len, max_len]` from `genome`, each
/// reverse-complemented with probability 1/2 and with a small substitution
/// rate to create novel k-mers.
pub fn sample_reads(
    rng: &mut Rng,
    genome: &[u8],
    n: usize,
    min_len: usize,
    max_len: usize,
    error_rate: f64,
) -> Vec<Vec<u8>> {
    (0..n)
        .map(|_| {
            let len = rng.range(min_len, max_len).min(genome.len());
            let s = rng.below((genome.len() - len + 1) as u64) as usize;
            let mut r = genome[s..s + len].to_vec();
            for c in r.iter_mut() {
                if rng.chance(error_rate) {
                    *c = rng.base();
                }
            }
            if rng.chance(0.5) {
                r = crate::codec::revcomp_seq(&r);
            }
            r
        })
        .collect()
}
