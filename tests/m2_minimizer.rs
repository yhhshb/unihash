//! Milestone 2 — minimizers with tied positions (hashed order).

mod common;
use common::seeds;
use std::collections::BTreeSet;
use unihash::codec::*;
use unihash::minimizer::*;
use unihash::rng::Rng;

/// Brute-force union of tied positions over all windows.
fn brute_piece_positions(vals: &[u64], w: usize) -> Vec<u32> {
    let mut set = BTreeSet::new();
    if vals.len() >= w {
        for s in 0..=vals.len() - w {
            let (_, pos) = minimizer_positions_in_window(vals, s, w);
            set.extend(pos.into_iter().map(|p| p as u32));
        }
    }
    set.into_iter().collect()
}

fn low_complexity_sequences(rng: &mut Rng) -> Vec<Vec<u8>> {
    let mut v = vec![
        vec![b'A'; 100],
        vec![b'T'; 100],
        b"AC".repeat(60),
        b"ACG".repeat(40),
        b"ACGT".repeat(30),
        b"AACCT".repeat(25),
        b"AT".repeat(50), // palindromic dinucleotide
        b"ACGTTGCA".repeat(15),
    ];
    for _ in 0..10 {
        let unit_len = rng.range(1, 7);
        let unit = rng.dna(unit_len);
        v.push(unit.repeat(120 / unit_len + 1));
        v.push(rng.dna_from(150, b"AT"));
    }
    v
}

#[test]
fn deque_equals_brute_force() {
    for seed in seeds(&[21, 22, 23]) {
        let mut rng = Rng::new(seed);
        // Random values from a small range → many ties.
        for w in 1..=20 {
            let n = rng.range(0, 120);
            let vals: Vec<u64> = (0..n).map(|_| rng.below(6)).collect();
            assert_eq!(
                window_minima(&vals, w),
                window_minima_brute(&vals, w),
                "seed={seed} w={w}"
            );
            assert_eq!(
                piece_minimizer_positions(&vals, w),
                brute_piece_positions(&vals, w),
                "seed={seed} w={w} positions"
            );
        }
        // Real ℓ-mer values from random and low-complexity DNA.
        let mut seqs = low_complexity_sequences(&mut rng);
        seqs.push(rng.dna(500));
        for seq in &seqs {
            for &(k, l) in &[(31, 17), (31, 15), (32, 19), (5, 3), (6, 6), (8, 1)] {
                let vals = lmer_values(seq, l);
                let w = k - l + 1;
                assert_eq!(window_minima(&vals, w), window_minima_brute(&vals, w));
                assert_eq!(
                    piece_minimizer_positions(&vals, w),
                    brute_piece_positions(&vals, w),
                    "seed={seed} k={k} l={l} seq={}",
                    String::from_utf8_lossy(seq)
                );
            }
        }
    }
}

#[test]
fn streaming_window_min_equals_brute_force() {
    for seed in seeds(&[24, 25]) {
        let mut rng = Rng::new(seed);
        for w in 1..=32 {
            let n = rng.range(0, 200);
            let vals: Vec<u64> = (0..n).map(|_| rng.below(5)).collect();
            let mut wm = WindowMin::new(w);
            for (i, &v) in vals.iter().enumerate() {
                wm.push(v);
                match wm.current() {
                    None => assert!(i + 1 < w),
                    Some((m, ties)) => {
                        let s = i + 1 - w;
                        let (bm, bpos) = minimizer_positions_in_window(&vals, s, w);
                        assert_eq!(m, bm, "seed={seed} w={w} s={s}");
                        let bmask = bpos.iter().fold(0u32, |a, &p| a | 1 << (p - s));
                        assert_eq!(ties, bmask, "seed={seed} w={w} s={s} ties");
                        assert_eq!(ties, window_tie_mask(&vals, s, w, m));
                    }
                }
            }
            wm.clear();
            assert!(wm.current().is_none());
        }
    }
}

#[test]
fn strand_invariance_of_window_minimizer() {
    for seed in seeds(&[31, 32]) {
        let mut rng = Rng::new(seed);
        for &(k, l) in &[
            (31, 17),
            (31, 15),
            (31, 19),
            (32, 19),
            (32, 32),
            (21, 11),
            (5, 3),
            (2, 1),
        ] {
            for _ in 0..1000 {
                let x = rng.next_u64() & kmer_mask(k);
                let r = revcomp_kmer(x, k);
                let (mx, px) = kmer_minimizer_positions(x, k, l);
                let (mr, pr) = kmer_minimizer_positions(r, k, l);
                assert_eq!(mx, mr, "seed={seed} k={k} l={l}");
                assert_eq!(kmer_minimizer(x, k, l), mx);
                let (mm, mask) = kmer_minimizer_mask(x, k, l);
                assert_eq!(mm, mx);
                assert_eq!(mask, px.iter().fold(0u32, |a, &p| a | 1 << p));
                // Tied positions mirror: p in x ↔ (k-ℓ-p) in rc(x).
                let mut mirrored: Vec<usize> = px.iter().map(|&p| k - l - p).collect();
                mirrored.sort_unstable();
                assert_eq!(mirrored, pr, "seed={seed} k={k} l={l} mirror");
            }
        }
    }
}

#[test]
fn ties_are_all_returned_on_low_complexity() {
    let k = 31;
    let l = 15;
    let w = k - l + 1;
    // poly-A: every position ties in every window.
    let s = vec![b'A'; 60];
    let vals = lmer_values(&s, l);
    for st in 0..=vals.len() - w {
        let (_, pos) = minimizer_positions_in_window(&vals, st, w);
        assert_eq!(pos, (st..st + w).collect::<Vec<_>>());
    }
    assert_eq!(
        piece_minimizer_positions(&vals, w),
        (0..vals.len() as u32).collect::<Vec<_>>()
    );
    // ACAC...: ℓ odd ⇒ ACA..A and CAC..C alternate; one class ties at
    // every other position.
    let s = b"AC".repeat(40);
    let vals = lmer_values(&s, l);
    let (m, pos) = minimizer_positions_in_window(&vals, 0, w);
    assert!(pos.len() >= w / 2, "expected ~half of positions tied");
    for p in &pos {
        assert_eq!(vals[*p], m);
    }
    let gaps: BTreeSet<usize> = pos.windows(2).map(|p| p[1] - p[0]).collect();
    assert_eq!(gaps, BTreeSet::from([2]));
}

#[test]
fn kmer_minimizer_matches_sequence_path() {
    for seed in seeds(&[41]) {
        let mut rng = Rng::new(seed);
        for &(k, l) in &[(31, 17), (32, 19), (5, 3), (4, 4)] {
            let seq = rng.dna(300);
            let vals = lmer_values(&seq, l);
            let mins = window_minima(&vals, k - l + 1);
            let mut roll = RollingKmer::new(k);
            for (j, &c) in seq.iter().enumerate() {
                if roll.push(encode(c)) {
                    let s = j + 1 - k;
                    assert_eq!(kmer_minimizer(roll.fwd(), k, l), mins[s]);
                    assert_eq!(kmer_minimizer(roll.canonical(), k, l), mins[s]);
                }
            }
        }
    }
}

#[test]
fn order_value_is_a_bijection_with_inverse() {
    let mut rng = Rng::new(99);
    for _ in 0..100_000 {
        let z = rng.next_u64();
        assert_eq!(unmix64(mix64(z)), z);
        assert_eq!(mix64(unmix64(z)), z);
        let x = z & kmer_mask(17);
        assert_eq!(lmer_of_order_value(order_value(x)), x);
    }
    for z in [0u64, 1, u64::MAX, 1 << 63] {
        assert_eq!(unmix64(mix64(z)), z);
    }
    // The all-A ℓ-mer must not be the global minimum by construction.
    assert_ne!(order_value(0), 0);
}
