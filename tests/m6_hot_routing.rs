//! Hot-minimizer routing (layers 2 and 3) and spilled index growth.

mod common;
use common::{check_against_oracle, check_invariants, seeds};
use unihash::codec::{canonical_kmers, kmer_mask};
use unihash::rng::{repeat_rich_genome, sample_reads, Rng};
use unihash::{CandidateScan, Config, KmerTable, Oracle};

fn build(c: &Config, reads: &[Vec<u8>], rng: &mut Rng) -> (KmerTable, Oracle) {
    let mut t = KmerTable::new(c.clone()).unwrap();
    let mut o = Oracle::new(c.k);
    for r in reads {
        // half of the reads through random chunks
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
        o.insert_read(r);
    }
    (t, o)
}

/// Stream every read through `query_stream` and compare with the oracle.
fn check_streaming(t: &KmerTable, o: &Oracle, reads: &[Vec<u8>], ctx: &str) {
    let mut sq = t.query_stream();
    for (ri, r) in reads.iter().enumerate() {
        let expect: Vec<u64> = canonical_kmers(r, t.k()).collect();
        let mut got = Vec::new();
        sq.push(r, |x, p| got.push((x, p)));
        sq.end();
        assert_eq!(got.len(), expect.len());
        for (&(x, p), &e) in got.iter().zip(&expect) {
            assert_eq!(x, e);
            assert_eq!(p, o.contains(x), "{ctx}: streaming read {ri}");
        }
    }
}

#[test]
fn routing_matches_oracle_with_low_thresholds() {
    for seed in seeds(&[901, 902, 903]) {
        let mut rng = Rng::new(seed);
        let genome = repeat_rich_genome(&mut rng, 30_000);
        let mut reads = sample_reads(&mut rng, &genome, 300, 40, 300, 0.01);
        // homopolymers and tandem repeats make ℓ2-minimizers hot as well
        reads.push(vec![b'A'; 400]);
        reads.push(b"ACGTTGCA".repeat(60));
        reads.push(b"AT".repeat(200));
        for &(k, l, l2) in &[(31, 17, 0), (31, 15, 21), (21, 11, 0), (32, 19, 32)] {
            for threshold in [1u32, 3, 16] {
                for scan in [CandidateScan::TieOffsets, CandidateScan::Range] {
                    let c = Config {
                        k,
                        l,
                        l2,
                        hot_threshold: threshold,
                        bucket_bits: 4,
                        scan,
                        ..Config::default()
                    };
                    let (t, o) = build(&c, &reads, &mut rng);
                    let ctx = format!("seed={seed} cfg={c:?}");
                    let st = t.stats();
                    assert!(st.insert.hot1 > 0, "{ctx}: no hot ℓ-minimizer");
                    // With ℓ2 = k an ℓ2-minimizer is the k-mer itself: it can
                    // only become hot through k-mers stored twice, which eager
                    // run indexing never does.
                    if threshold <= 3 && c.fallback_l() != Some(k) {
                        assert!(st.insert.hot2 > 0, "{ctx}: no hot ℓ2-minimizer");
                        assert!(st.layer3_kmers > 0, "{ctx}: layer 3 unused");
                    }
                    check_against_oracle(&t, &o, &reads, &mut rng, 5_000, &ctx);
                    check_invariants(&t, &o, &ctx);
                    check_streaming(&t, &o, &reads, &ctx);
                }
            }
        }
    }
}

#[test]
fn routing_small_k_exhaustive() {
    for seed in seeds(&[911, 912]) {
        let mut rng = Rng::new(seed);
        for &(k, l, l2) in &[(5, 2, 3), (6, 3, 5), (7, 3, 0), (4, 1, 4)] {
            for threshold in [1u32, 2, 5] {
                let c = Config {
                    k,
                    l,
                    l2,
                    hot_threshold: threshold,
                    bucket_bits: 0,
                    max_piece_len: k as u32 + 4,
                    ..Config::default()
                };
                let mut reads = Vec::new();
                for _ in 0..300 {
                    let len = rng.range(0, 30);
                    reads.push(match rng.below(3) {
                        0 => rng.dna_from(len, b"AT"),
                        1 => rng.dna_from(len, b"ACG"),
                        _ => rng.dna(len),
                    });
                }
                let (t, o) = build(&c, &reads, &mut rng);
                let ctx = format!("seed={seed} cfg={c:?}");
                if threshold == 1 {
                    let st = t.stats();
                    assert!(st.insert.hot1 > 0, "{ctx}: {st:?}");
                }
                for x in 0..=kmer_mask(k) {
                    assert_eq!(t.contains(x), o.contains(x), "{ctx}: k-mer {x:#x}");
                }
                check_invariants(&t, &o, &ctx);
                check_streaming(&t, &o, &reads, &ctx);
            }
        }
    }
}

#[test]
fn routing_disabled_is_identical_to_no_hot_layer() {
    let mut rng = Rng::new(921);
    let genome = repeat_rich_genome(&mut rng, 20_000);
    let reads = sample_reads(&mut rng, &genome, 200, 40, 200, 0.01);
    let c = Config {
        hot_threshold: 0,
        bucket_bits: 4,
        ..Config::default()
    };
    let (t, o) = build(&c, &reads, &mut rng);
    let st = t.stats();
    assert_eq!((st.insert.hot1, st.entries2, st.layer3_kmers), (0, 0, 0));
    check_against_oracle(&t, &o, &reads, &mut rng, 2_000, "disabled");
}

#[test]
fn disk_spilled_growth_matches_memory_spill() {
    let mut rng = Rng::new(931);
    let genome = repeat_rich_genome(&mut rng, 40_000);
    let reads = sample_reads(&mut rng, &genome, 400, 50, 300, 0.01);
    let dir = std::env::temp_dir().join(format!("unihash_spill_test_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let base = Config {
        bucket_bits: 0,
        hot_threshold: 8,
        ..Config::default()
    };
    let mut mem = KmerTable::new(base.clone()).unwrap();
    let mut disk = KmerTable::new(Config {
        spill_dir: Some(dir.clone()),
        ..base
    })
    .unwrap();
    let mut o = Oracle::new(31);
    for r in &reads {
        mem.insert_read(r);
        disk.insert_read(r);
        o.insert_read(r);
    }
    let (sm, sd) = (mem.stats(), disk.stats());
    assert!(sm.index_grows > 5);
    assert!(sm.spill.bytes_in_memory > 0 && sm.spill.bytes_to_disk == 0);
    assert!(sd.spill.bytes_to_disk > 0 && sd.spill.bytes_in_memory == 0);
    assert_eq!(sm.spill.bytes_in_memory, sd.spill.bytes_to_disk);
    // Same placement either way: identical entries bucket by bucket.
    let em: Vec<_> = mem.index().iter_entries_with_bucket().collect();
    let ed: Vec<_> = disk.index().iter_entries_with_bucket().collect();
    assert_eq!(em, ed);
    check_against_oracle(&disk, &o, &reads, &mut rng, 2_000, "disk spill");
    // Spill files are removed after each rebuild.
    assert_eq!(std::fs::read_dir(&dir).unwrap().count(), 0);
    std::fs::remove_dir_all(&dir).unwrap();
}
