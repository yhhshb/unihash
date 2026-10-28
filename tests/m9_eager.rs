//! Eager run indexing against the deferred (original) insertion path: same
//! k-mer set and answers, no stored duplicates, k-mers visible while a read
//! is still open, and identical pieces whatever the chunking.

mod common;
use common::{check_against_oracle, seeds};
use std::collections::HashSet;
use unihash::codec::{canonical_kmers, revcomp_seq};
use unihash::rng::{repeat_rich_genome, sample_reads, Rng};
use unihash::{CandidateScan, Config, KmerTable, Oracle, RunIndexing};

fn configs() -> Vec<Config> {
    vec![
        Config::default(),
        Config {
            k: 21,
            l: 11,
            bucket_bits: 2,
            hot_threshold: 3,
            max_piece_len: 60,
            ..Config::default()
        },
        Config {
            k: 31,
            l: 15,
            hot_threshold: 1, // everything routes to layers 2/3
            scan: CandidateScan::Range,
            ..Config::default()
        },
        Config {
            k: 32,
            l: 19,
            l2: 32,
            hot_threshold: 2,
            bucket_bits: 0,
            ..Config::default()
        },
        Config {
            k: 7,
            l: 3,
            hot_threshold: 4,
            bucket_bits: 0,
            max_piece_len: 9,
            ..Config::default()
        },
        Config {
            k: 15,
            l: 9,
            hot_threshold: 0, // routing off
            ..Config::default()
        },
    ]
}

fn with(c: &Config, r: RunIndexing) -> Config {
    Config {
        run_indexing: r,
        ..c.clone()
    }
}

/// Insert `r`, cutting it into random chunks with probability 1/2.
fn insert(t: &mut KmerTable, r: &[u8], rng: &mut Rng) {
    if rng.chance(0.5) {
        let mut i = 0;
        while i < r.len() {
            let j = (i + rng.range(1, 50)).min(r.len());
            t.insert_chunk(&r[i..j]);
            i = j;
        }
        t.end_read();
    } else {
        t.insert_read(r);
    }
}

fn stored(t: &KmerTable) -> (HashSet<u64>, u64) {
    let mut s = HashSet::new();
    let mut n = 0;
    t.for_each_stored_kmer(|_, _, x| {
        s.insert(x);
        n += 1;
    });
    (s, n)
}

fn test_reads(rng: &mut Rng) -> Vec<Vec<u8>> {
    let genome = repeat_rich_genome(rng, 30_000);
    let mut reads = sample_reads(rng, &genome, 600, 20, 400, 0.01);
    // Intra-read repeats: stored twice by the deferred path, never by eager.
    reads.push(vec![b'A'; 200]);
    reads.push(b"ACGTTGCA".repeat(40));
    reads.push(b"ACGTNNACGTACGTNACGT".repeat(10));
    let mut tandem = rng.dna(37);
    tandem = tandem.repeat(9);
    reads.push(revcomp_seq(&tandem));
    reads.push(tandem);
    reads.push(genome[..5000].to_vec());
    reads
}

#[test]
fn eager_matches_deferred_and_oracle() {
    for seed in seeds(&[2001, 2002, 2003]) {
        for c in configs() {
            let ctx = format!("seed={seed} cfg={c:?}");
            let mut rng = Rng::new(seed);
            let reads = test_reads(&mut rng);
            let mut d = KmerTable::new(with(&c, RunIndexing::Deferred)).unwrap();
            let mut e = KmerTable::new(with(&c, RunIndexing::Eager)).unwrap();
            let mut o = Oracle::new(c.k);
            for r in &reads {
                insert(&mut d, r, &mut rng);
                insert(&mut e, r, &mut rng);
                o.insert_read(r);
            }
            for (t, name) in [(&d, "deferred"), (&e, "eager")] {
                let ctx = format!("{ctx} {name}");
                check_against_oracle(t, &o, &reads, &mut rng, 2000, &ctx);
                if let Err(err) = t.verify_index() {
                    panic!("{ctx}: {err}");
                }
                let (set, occ) = stored(t);
                assert_eq!(set.len(), o.len(), "{ctx}: distinct stored k-mers");
                assert_eq!(t.stats().stored_kmers, occ, "{ctx}: stored_kmers");
            }
            let (_, occ_e) = stored(&e);
            assert_eq!(occ_e, o.len() as u64, "{ctx}: eager stored a duplicate");
            assert_eq!(e.redundant_kmers(), 0, "{ctx}");
            // Every insertion lookup sees the same set in both modes, except
            // that deferred misses k-mers repeated inside its own open run.
            let (si, ei) = (d.insert_stats(), e.insert_stats());
            assert_eq!(si.kmers_processed, ei.kmers_processed, "{ctx}");
            assert_eq!(
                ei.present_hits,
                si.present_hits + d.redundant_kmers(),
                "{ctx}"
            );
            assert!(
                d.redundant_kmers() > 0,
                "{ctx}: test has no intra-run repeat"
            );
        }
    }
}

/// Without repeats or hot minimizers the two paths build the very same
/// pieces and index entries.
#[test]
fn eager_equals_deferred_without_repeats() {
    for seed in seeds(&[2011, 2012]) {
        let c = Config {
            hot_threshold: 0,
            max_piece_len: 500,
            ..Config::default()
        };
        let mut rng = Rng::new(seed);
        let genome = rng.dna(100_000);
        let reads = sample_reads(&mut rng, &genome, 1000, 50, 300, 0.0);
        let mut d = KmerTable::new(with(&c, RunIndexing::Deferred)).unwrap();
        let mut e = KmerTable::new(with(&c, RunIndexing::Eager)).unwrap();
        for r in &reads {
            insert(&mut d, r, &mut rng);
            insert(&mut e, r, &mut rng);
        }
        assert_eq!(d.redundant_kmers(), 0, "seed={seed}: random genome repeats");
        let pd: Vec<_> = d.pieces().iter().collect();
        let pe: Vec<_> = e.pieces().iter().collect();
        assert_eq!(pd, pe, "seed={seed}: pieces differ");
        let mut id: Vec<_> = d
            .index()
            .iter_entries()
            .map(|x| (x.piece_id, x.q))
            .collect();
        let mut ie: Vec<_> = e
            .index()
            .iter_entries()
            .map(|x| (x.piece_id, x.q))
            .collect();
        id.sort_unstable();
        ie.sort_unstable();
        assert_eq!(id, ie, "seed={seed}: index entries differ");
    }
}

/// Eager pieces do not depend on where the input is cut into chunks.
#[test]
fn eager_pieces_independent_of_chunking() {
    for seed in seeds(&[2021, 2022]) {
        for c in configs() {
            let c = with(&c, RunIndexing::Eager);
            let ctx = format!("seed={seed} cfg={c:?}");
            let mut rng = Rng::new(seed);
            let reads = test_reads(&mut rng);
            let mut whole = KmerTable::new(c.clone()).unwrap();
            let mut chunked = KmerTable::new(c.clone()).unwrap();
            for r in &reads {
                whole.insert_read(r);
                let mut i = 0;
                while i < r.len() {
                    let j = (i + rng.range(1, 7)).min(r.len());
                    chunked.insert_chunk(&r[i..j]);
                    i = j;
                }
                chunked.end_read();
            }
            let a: Vec<_> = whole.pieces().iter().collect();
            let b: Vec<_> = chunked.pieces().iter().collect();
            assert_eq!(a, b, "{ctx}");
        }
    }
}

/// With eager indexing, the k-mers of a read are queryable before
/// `end_read`, as soon as the chunk holding their last base was inserted.
#[test]
fn eager_kmers_visible_in_open_read() {
    for seed in seeds(&[2031, 2032]) {
        for c in configs() {
            let c = with(&c, RunIndexing::Eager);
            let ctx = format!("seed={seed} cfg={c:?}");
            let mut rng = Rng::new(seed);
            let reads = test_reads(&mut rng);
            let mut t = KmerTable::new(c.clone()).unwrap();
            let mut o = Oracle::new(c.k);
            for r in reads.iter().take(80) {
                let mut i = 0;
                while i < r.len() {
                    let j = (i + rng.range(1, 60)).min(r.len());
                    t.insert_chunk(&r[i..j]);
                    for x in canonical_kmers(&r[..j], c.k) {
                        assert!(t.contains(x), "{ctx}: k-mer of open read not visible");
                    }
                    i = j;
                }
                t.end_read();
                o.insert_read(r);
            }
            check_against_oracle(&t, &o, &reads[..80], &mut rng, 1000, &ctx);
        }
    }
}

/// Deletion and both defragmentations keep an eagerly built table exact and
/// duplicate-free.
#[test]
fn eager_with_deletion_and_defrag() {
    for seed in seeds(&[2041, 2042]) {
        for c in configs() {
            let c = with(&c, RunIndexing::Eager);
            let ctx = format!("seed={seed} cfg={c:?}");
            let mut rng = Rng::new(seed);
            let reads = test_reads(&mut rng);
            let mut t = KmerTable::new(c.clone()).unwrap();
            let mut o = Oracle::new(c.k);
            for (i, r) in reads.iter().enumerate() {
                insert(&mut t, r, &mut rng);
                o.insert_read(r);
                if i % 3 == 0 {
                    let s = &reads[rng.below(i as u64 + 1) as usize];
                    let a = rng.range(0, s.len() / 2);
                    let piece = &s[a..s.len().min(a + 80)];
                    t.delete_read(piece);
                    o.delete_read(piece);
                }
                match i % 97 {
                    50 => {
                        t.defrag_incremental();
                    }
                    96 => {
                        t.defrag();
                    }
                    _ => {}
                }
            }
            for r in &reads {
                for x in canonical_kmers(r, c.k) {
                    assert_eq!(t.contains(x), o.contains(x), "{ctx}");
                }
            }
            for &x in o.iter() {
                assert!(t.contains(x), "{ctx}: oracle k-mer missing");
            }
            if let Err(err) = t.verify_index() {
                panic!("{ctx}: {err}");
            }
            let (set, occ) = stored(&t);
            assert_eq!(set.len(), o.len(), "{ctx}");
            assert_eq!(occ, o.len() as u64, "{ctx}: duplicate after edits");
        }
    }
}
