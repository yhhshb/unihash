//! Packed piece store and flat minimizer index.

mod common;
use common::seeds;
use std::collections::{BTreeSet, HashMap};
use unihash::codec::{canonicalize, encode_kmer};
use unihash::index::{MinimizerIndex, ScanTrace, BLOCK_SLOTS, MAX_LOAD};
use unihash::piece_store::PieceStore;
use unihash::rng::Rng;

#[test]
fn piece_round_trip() {
    for seed in seeds(&[51, 52]) {
        let mut rng = Rng::new(seed);
        let mut store = PieceStore::new();
        let mut v = Vec::new();
        for _ in 0..300 {
            // lengths around word boundaries and arbitrary offsets
            let n = match rng.below(3) {
                0 => rng.range(1, 3),
                1 => rng.range(30, 34),
                _ => rng.range(1, 300),
            };
            let s = rng.dna(n);
            let id = store.add_piece(&s);
            assert_eq!(id as usize, v.len());
            v.push(s);
        }
        assert_eq!(store.len(), v.len());
        let total: usize = v.iter().map(|s| s.len()).sum();
        assert_eq!(store.total_bases(), total as u64);
        for (i, s) in v.iter().enumerate() {
            let id = i as u32;
            assert_eq!(store.piece_seq(id), &s[..], "seed={seed} piece {i}");
            assert_eq!(store.length(id) as usize, s.len());
            for (j, &c) in s.iter().enumerate() {
                assert_eq!(store.get_base(id, j as u32), c);
            }
            // every k-mer of every piece, k = 1..=32 (spans word boundaries)
            for k in [1usize, 2, 7, 16, 17, 31, 32] {
                if s.len() < k {
                    continue;
                }
                for st in 0..=s.len() - k {
                    let direct = encode_kmer(&s[st..st + k]).unwrap();
                    assert_eq!(store.kmer_at(id, st as u32, k), direct, "seed={seed} k={k}");
                    assert_eq!(
                        store.canonical_kmer_at(id, st as u32, k),
                        canonicalize(direct, k)
                    );
                }
            }
        }
        assert!(store.iter().zip(v.iter()).all(|(a, b)| a == *b));
        // 2 bits per base plus one pad word, plus 8 bytes per piece.
        let (used, _) = store.heap_bytes();
        assert_eq!(used, (total.div_ceil(32) + 1) * 8 + v.len() * 8);
    }
}

#[test]
fn piece_lowercase_is_normalised() {
    let mut store = PieceStore::new();
    let id = store.add_piece(b"acgtACGT");
    assert_eq!(store.piece_seq(id), b"ACGTACGT");
}

#[test]
#[should_panic(expected = "non-ACGT")]
fn piece_rejects_non_acgt() {
    PieceStore::new().add_piece(b"ACGNT");
}

/// Visited set of `for_each_candidate(m)` must be exactly the entries whose
/// tag equals `m`'s tag and that sit in one of `m`'s two bucket chains.
fn check_candidates(idx: &MinimizerIndex, triples: &[(u64, u32, u32)], ms: &[u64], ctx: &str) {
    let placed: Vec<(usize, u8, u32, u32)> = idx
        .iter_entries_with_bucket()
        .map(|(b, e)| (b, e.tag, e.piece_id, e.q))
        .collect();
    assert_eq!(placed.len(), triples.len(), "{ctx}: entry count");
    for &m in ms {
        let s = idx.slot(m);
        let expected: BTreeSet<(u32, u32)> = placed
            .iter()
            .filter(|&&(b, t, _, _)| t == s.tag && (b == s.b1 || b == s.b2))
            .map(|&(_, _, p, q)| (p, q))
            .collect();
        // all true entries of m are included
        for &(m2, p, q) in triples {
            if m2 == m {
                assert!(
                    expected.contains(&(p, q)),
                    "{ctx}: true entry of m={m} not in its chains"
                );
            }
        }
        let mut got = BTreeSet::new();
        let mut tr = ScanTrace::default();
        let stopped = idx.for_each_candidate(m, &mut tr, |p, q| {
            assert!(got.insert((p, q)), "{ctx}: duplicate visit");
            false
        });
        assert!(!stopped);
        assert_eq!(got, expected, "{ctx}: m={m}");
        assert_eq!(tr.tag_matches, expected.len() as u64);
        let mut cached = Vec::new();
        idx.collect_candidates(m, &mut ScanTrace::default(), &mut cached);
        assert_eq!(cached.into_iter().collect::<BTreeSet<_>>(), expected);
    }
}

#[test]
fn index_visits_exactly_matching_chains_and_tag() {
    for seed in seeds(&[61, 62]) {
        let mut rng = Rng::new(seed);
        // Few buckets and a modest value range: bucket collisions, tag
        // collisions, repeated minimizers and overflow chains.
        for n in [1usize, 2, 3, 17, 64] {
            let mut idx = MinimizerIndex::with_buckets(n);
            let mut triples = Vec::new();
            for i in 0..3000u32 {
                let m = rng.below(400);
                let q = rng.below(1000) as u32;
                idx.insert(m, i, q);
                triples.push((m, i, q));
            }
            assert_eq!(idx.num_entries(), triples.len() as u64);
            assert_eq!(idx.num_buckets(), n);
            assert!(idx.overflow_blocks() > 0);
            let ms: Vec<u64> = (0..450).collect();
            check_candidates(&idx, &triples, &ms, &format!("seed={seed} n={n}"));
        }
    }
}

#[test]
fn index_grows_and_keeps_every_entry() {
    for seed in seeds(&[63]) {
        let mut rng = Rng::new(seed);
        // m is a function of (pid, q) here, as it is for real pieces.
        // Mostly distinct minimizers, some repeated (every 16th pid shares).
        let m_of = |p: u32, q: u32| ((p as u64 & !15) << 8) | (q as u64 % 50);
        let mut idx = MinimizerIndex::new(0);
        let mut triples = Vec::new();
        for i in 0..40_000u32 {
            let q = rng.below(100) as u32;
            idx.reserve(1, m_of);
            idx.insert(m_of(i, q), i, q);
            triples.push((m_of(i, q), i, q));
            // load never exceeds the limit after a reserve
            assert!(
                idx.num_entries() as f64
                    <= MAX_LOAD * (BLOCK_SLOTS * idx.num_buckets()) as f64 + 1.0
            );
        }
        assert!(idx.grows() > 5, "grew {} times", idx.grows());
        let ms: Vec<u64> = triples.iter().step_by(37).map(|t| t.0).collect();
        check_candidates(&idx, &triples, &ms, &format!("seed={seed} grown"));
        // Compaction (chains laid out contiguously) changes nothing visible.
        idx.shrink_to_fit();
        check_candidates(&idx, &triples, &ms, &format!("seed={seed} compacted"));
        // Two choices keep heads evenly filled: few overflow blocks.
        assert!(
            (idx.overflow_blocks() as f64) < 0.1 * idx.num_buckets() as f64,
            "overflow {} of {} buckets",
            idx.overflow_blocks(),
            idx.num_buckets()
        );
    }
}

#[test]
fn hot_minimizer_uses_overflow_chain() {
    let mut idx = MinimizerIndex::with_buckets(64);
    for i in 0..1000u32 {
        idx.insert(42, i, i);
    }
    let mut n = 0;
    idx.for_each_candidate(42, &mut ScanTrace::default(), |_, _| {
        n += 1;
        false
    });
    assert_eq!(n, 1000);
    idx.shrink_to_fit();
    let mut n2 = 0;
    idx.for_each_candidate(42, &mut ScanTrace::default(), |_, _| {
        n2 += 1;
        false
    });
    assert_eq!(n2, 1000);
    // 1000 entries: two full heads, the rest in full overflow blocks.
    assert_eq!(
        idx.overflow_blocks(),
        (1000 - 2 * BLOCK_SLOTS).div_ceil(BLOCK_SLOTS)
    );
    let mut per_bucket: HashMap<usize, usize> = HashMap::new();
    for (b, _) in idx.iter_entries_with_bucket() {
        *per_bucket.entry(b).or_default() += 1;
    }
    let s = idx.slot(42);
    assert!(per_bucket.keys().all(|&b| b == s.b1 || b == s.b2));
}

#[test]
fn index_early_stop() {
    let mut idx = MinimizerIndex::new(0);
    for i in 0..10 {
        idx.insert(5, i, i);
    }
    let mut n = 0;
    let mut tr = ScanTrace::default();
    assert!(idx.for_each_candidate(5, &mut tr, |_, _| {
        n += 1;
        n == 3
    }));
    assert_eq!(n, 3);
}

#[test]
fn drain_removes_exactly_one_minimizer_and_survives_growth() {
    for seed in seeds(&[64]) {
        let mut rng = Rng::new(seed);
        let m_of = |p: u32, _q: u32| (p as u64) % 300;
        let mut idx = MinimizerIndex::with_buckets(8);
        let mut triples = Vec::new();
        for i in 0..6000u32 {
            let q = rng.below(100) as u32;
            idx.reserve(1, m_of);
            idx.insert(m_of(i, q), i, q);
            triples.push((m_of(i, q), i, q));
        }
        let before = idx.num_entries();
        for hot in [7u64, 123, 299] {
            assert!(!idx.is_hot(hot));
            let expected = triples.iter().filter(|t| t.0 == hot).count();
            assert_eq!(idx.count_minimizer(hot, m_of), expected);
            let drained = idx.drain_minimizer(hot, m_of);
            assert_eq!(drained.len(), expected, "m={hot}");
            assert!(drained.iter().all(|&(p, q)| m_of(p, q) == hot));
            assert!(idx.is_hot(hot));
            assert_eq!(idx.count_minimizer(hot, m_of), 0);
            triples.retain(|t| t.0 != hot);
        }
        assert_eq!(idx.num_entries(), before - (6000 - triples.len()) as u64);
        let ms: Vec<u64> = (0..300).collect();
        check_candidates(&idx, &triples, &ms, &format!("seed={seed} drained"));
        // Growth keeps hot minimizers hot and the remaining entries intact.
        for i in 6000..20000u32 {
            let m = m_of(i, 0);
            if [7u64, 123, 299].contains(&m) {
                continue;
            }
            idx.reserve(1, m_of);
            idx.insert(m, i, 0);
            triples.push((m, i, 0));
        }
        assert!(idx.grows() > 0);
        assert!(idx.is_hot(7) && idx.is_hot(123) && idx.is_hot(299) && !idx.is_hot(8));
        check_candidates(&idx, &triples, &ms, &format!("seed={seed} drained+grown"));
        idx.shrink_to_fit();
        check_candidates(&idx, &triples, &ms, &format!("seed={seed} compacted"));
    }
}
