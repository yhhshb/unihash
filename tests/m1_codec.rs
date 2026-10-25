//! Codec and piece sequence access.

mod common;
use common::seeds;
use unihash::codec::*;
use unihash::piece_store::PieceStore;
use unihash::rng::Rng;

#[test]
fn encode_round_trip() {
    for c in 0u8..=255 {
        let b = encode(c);
        match c {
            b'A' | b'C' | b'G' | b'T' => assert_eq!(decode(b), c),
            b'a' | b'c' | b'g' | b't' => assert_eq!(decode(b), c.to_ascii_uppercase()),
            _ => assert_eq!(b, INVALID, "char {c} should be invalid"),
        }
    }
    assert_eq!(encode(b'A'), 0);
    assert_eq!(encode(b'C'), 1);
    assert_eq!(encode(b'G'), 2);
    assert_eq!(encode(b'T'), 3);
}

#[test]
fn rolling_equals_direct_all_k() {
    for seed in seeds(&[1, 2, 3]) {
        let mut rng = Rng::new(seed);
        for k in 1..=32 {
            let seq = rng.dna(200);
            let mut roll = RollingKmer::new(k);
            for (j, &c) in seq.iter().enumerate() {
                if roll.push(encode(c)) {
                    let s = j + 1 - k;
                    let win = &seq[s..s + k];
                    let f = encode_kmer(win).unwrap();
                    let r = encode_kmer(&revcomp_seq(win)).unwrap();
                    assert_eq!(roll.fwd(), f, "seed={seed} k={k} s={s} fwd");
                    assert_eq!(roll.rc(), r, "seed={seed} k={k} s={s} rc");
                    assert_eq!(roll.canonical(), f.min(r), "seed={seed} k={k} canonical");
                    assert_eq!(revcomp_kmer(f, k), r, "seed={seed} k={k} revcomp_kmer");
                    assert_eq!(encode_kmer_both(win), (f, r));
                }
            }
        }
    }
}

#[test]
fn kmer_at_equals_direct_around_word_boundaries() {
    for seed in seeds(&[7, 8]) {
        let mut rng = Rng::new(seed);
        let mut store = PieceStore::new();
        let mut originals = Vec::new();
        for len in (30..=34).chain(62..=66).chain(126..=130) {
            let s = rng.dna(len);
            originals.push(s.clone());
            store.add_piece(&s);
        }
        for k in [2usize, 5, 16, 17, 31, 32] {
            for (pid, s) in originals.iter().enumerate() {
                let pid = pid as u32;
                assert_eq!(store.length(pid) as usize, s.len());
                if s.len() < k {
                    continue;
                }
                for st in 0..=s.len() - k {
                    let direct = encode_kmer(&s[st..st + k]).unwrap();
                    assert_eq!(
                        store.kmer_at(pid, st as u32, k),
                        direct,
                        "seed={seed} k={k}"
                    );
                    assert_eq!(
                        store.canonical_kmer_at(pid, st as u32, k),
                        canonicalize(direct, k)
                    );
                }
            }
        }
    }
}

#[test]
fn canonical_is_min_of_fwd_and_rc() {
    for seed in seeds(&[11]) {
        let mut rng = Rng::new(seed);
        for k in 1..=32 {
            for _ in 0..200 {
                let x = rng.next_u64() & kmer_mask(k);
                let r = revcomp_kmer(x, k);
                assert_eq!(canonicalize(x, k), x.min(r));
                assert_eq!(
                    canonicalize(x, k),
                    canonicalize(r, k),
                    "strand independence"
                );
            }
        }
    }
}

#[test]
fn revcomp_involution() {
    for seed in seeds(&[13]) {
        let mut rng = Rng::new(seed);
        for k in 1..=32 {
            for _ in 0..200 {
                let x = rng.next_u64() & kmer_mask(k);
                assert_eq!(revcomp_kmer(revcomp_kmer(x, k), k), x, "seed={seed} k={k}");
                assert!(revcomp_kmer(x, k) <= kmer_mask(k));
                let s = decode_kmer(x, k);
                assert_eq!(encode_kmer(s.as_bytes()), Some(x));
                assert_eq!(
                    encode_kmer(&revcomp_seq(s.as_bytes())),
                    Some(revcomp_kmer(x, k))
                );
            }
        }
        let s = rng.dna(1000);
        assert_eq!(revcomp_seq(&revcomp_seq(&s)), s);
    }
}

#[test]
fn k32_extremes() {
    let all_t = encode_kmer(&[b'T'; 32]).unwrap();
    assert_eq!(all_t, u64::MAX);
    assert_eq!(revcomp_kmer(all_t, 32), 0);
    assert_eq!(canonicalize(all_t, 32), 0);
    let mut r = RollingKmer::new(32);
    for _ in 0..40 {
        r.push(3);
    }
    assert_eq!((r.fwd(), r.rc()), (u64::MAX, 0));
}

#[test]
fn canonical_kmers_skip_non_acgt() {
    let v: Vec<u64> = canonical_kmers(b"ACGNACGTxAC", 3).collect();
    let expect: Vec<u64> = ["ACG", "ACG", "CGT"]
        .iter()
        .map(|s| canonicalize(encode_kmer(s.as_bytes()).unwrap(), 3))
        .collect();
    assert_eq!(v, expect);
    // lowercase is accepted
    let a: Vec<u64> = canonical_kmers(b"acgtac", 4).collect();
    let b: Vec<u64> = canonical_kmers(b"ACGTAC", 4).collect();
    assert_eq!(a, b);
}
