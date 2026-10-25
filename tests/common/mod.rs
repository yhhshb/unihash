#![allow(dead_code)]
//! Shared helpers for integration tests.
//!
//! Every randomized test iterates over a fixed list of seeds and prints the
//! seed in every assertion message. To replay one failing case:
//!
//! ```text
//! UNIHASH_SEED=1234 cargo test <test_name>
//! ```

use unihash::codec::{canonical_kmers, encode_kmer};
use unihash::rng::Rng;
use unihash::{KmerTable, Oracle};

pub fn seeds(defaults: &[u64]) -> Vec<u64> {
    match std::env::var("UNIHASH_SEED") {
        Ok(s) => vec![s.trim().parse().expect("UNIHASH_SEED must be a u64")],
        Err(_) => defaults.to_vec(),
    }
}

/// Compare table and oracle on: every k-mer of every read, every oracle
/// k-mer, single-substitution neighbours of read k-mers (near misses that
/// share minimizers), and random k-mers. Returns the number of queries.
pub fn check_against_oracle(
    t: &KmerTable,
    o: &Oracle,
    reads: &[Vec<u8>],
    rng: &mut Rng,
    n_random: usize,
    ctx: &str,
) -> usize {
    let k = t.k();
    let mut n = 0usize;
    for (ri, r) in reads.iter().enumerate() {
        for x in canonical_kmers(r, k) {
            assert!(
                t.contains(x),
                "{ctx}: false negative for k-mer {} of read {ri}",
                unihash::codec::decode_kmer(x, k)
            );
            n += 1;
            // near miss: substitute one base
            let pos = rng.below(k as u64) as usize;
            let nb = rng.below(4);
            let y = (x & !(3u64 << (2 * pos))) | (nb << (2 * pos));
            assert_eq!(
                t.contains(y),
                o.contains(y),
                "{ctx}: mismatch on near-miss k-mer {}",
                unihash::codec::decode_kmer(y, k)
            );
            n += 1;
        }
    }
    for &x in o.iter() {
        assert!(t.contains(x), "{ctx}: oracle k-mer missing from table");
        n += 1;
    }
    for _ in 0..n_random {
        let x = rng.next_u64() & unihash::codec::kmer_mask(k);
        assert_eq!(
            t.contains(x),
            o.contains(x),
            "{ctx}: mismatch on random k-mer {}",
            unihash::codec::decode_kmer(x, k)
        );
        n += 1;
    }
    n
}

/// Structural invariants of a table built from the same reads as `o`.
pub fn check_invariants(t: &KmerTable, o: &Oracle, ctx: &str) {
    let k = t.k();
    let max = t.config().max_piece_len;
    for (pid, p) in t.pieces().iter().enumerate() {
        assert!(p.len() >= k, "{ctx}: piece {pid} shorter than k");
        assert!(
            p.len() <= max as usize,
            "{ctx}: piece {pid} exceeds max_piece_len"
        );
    }
    // Every absent k-mer goes to exactly one piece; duplicates can only be
    // intra-run repeats, i.e. within the same piece.
    let mut owner = std::collections::HashMap::new();
    t.for_each_stored_kmer(|pid, _, x| {
        let prev = *owner.entry(x).or_insert(pid);
        assert_eq!(
            prev,
            pid,
            "{ctx}: k-mer {} stored in two different pieces",
            unihash::codec::decode_kmer(x, k)
        );
        assert!(o.contains(x), "{ctx}: stored k-mer not in oracle");
    });
    assert_eq!(
        owner.len(),
        o.len(),
        "{ctx}: distinct stored k-mers != oracle size"
    );
    let st = t.stats();
    assert_eq!(
        st.stored_kmers - t.redundant_kmers(),
        o.len() as u64,
        "{ctx}: stored - redundant != distinct"
    );
}

pub fn kmer(s: &str) -> u64 {
    encode_kmer(s.as_bytes()).unwrap()
}
