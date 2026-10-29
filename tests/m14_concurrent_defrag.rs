//! Defragmentation of `ConcurrentKmerTable` in place (no conversion to
//! `KmerTable`): the same pieces as through a conversion, the same k-mer
//! set as the oracle, an exact index, and a table that keeps working with
//! concurrent writers and deleters afterwards.

mod common;
use common::seeds;
use unihash::codec::{canonical_kmers, kmer_mask};
use unihash::rng::{repeat_rich_genome, sample_reads, Rng};
use unihash::{CandidateScan, ConcurrentKmerTable, Config, KmerTable, Oracle};

fn configs() -> Vec<Config> {
    vec![
        Config::default(),
        Config {
            k: 21,
            l: 11,
            bucket_bits: 0,
            hot_threshold: 3,
            max_piece_len: 300,
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
            k: 15,
            l: 9,
            hot_threshold: 0,
            ..Config::default()
        },
    ]
}

/// Reads, and stretches of them to delete.
fn workload(rng: &mut Rng, k: usize) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    let g = repeat_rich_genome(rng, 15_000);
    let reads = sample_reads(rng, &g, 400, 40, 300, 0.01);
    let dels = reads
        .iter()
        .step_by(3)
        .filter(|r| r.len() > k + 5)
        .map(|r| {
            let a = rng.below((r.len() - k) as u64) as usize;
            let b = (a + k + rng.below(60) as usize).min(r.len());
            r[a..b].to_vec()
        })
        .collect();
    (reads, dels)
}

/// One writer, then sequential deletions: a deterministic layout.
fn build_one_writer(
    c: &Config,
    reads: &[Vec<u8>],
    dels: &[Vec<u8>],
    shards: usize,
) -> ConcurrentKmerTable {
    let t = ConcurrentKmerTable::with_shards(c.clone(), shards).unwrap();
    {
        let mut w = t.writer();
        for r in reads {
            w.insert_read(r);
        }
    }
    for d in dels {
        t.delete_read(d);
    }
    t
}

fn same_as_oracle(
    t: &ConcurrentKmerTable,
    o: &Oracle,
    reads: &[Vec<u8>],
    rng: &mut Rng,
    ctx: &str,
) {
    let k = t.k();
    for r in reads {
        for x in canonical_kmers(r, k) {
            assert_eq!(t.contains(x), o.contains(x), "{ctx}: read k-mer");
        }
    }
    for &x in o.iter() {
        assert!(t.contains(x), "{ctx}: oracle k-mer missing");
    }
    for _ in 0..3000 {
        let x = rng.next_u64() & kmer_mask(k);
        assert_eq!(t.contains(x), o.contains(x), "{ctx}: random k-mer");
    }
}

fn check_table(t: KmerTable, o: &Oracle, ctx: &str) {
    t.verify_index().unwrap_or_else(|e| panic!("{ctx}: {e}"));
    assert_eq!(t.redundant_kmers(), 0, "{ctx}: k-mer stored twice");
    assert_eq!(t.distinct_kmers(), o.len() as u64, "{ctx}");
}

#[test]
fn native_defrag_gives_the_conversion_pieces() {
    for seed in seeds(&[1401, 1402]) {
        for c in configs() {
            for shards in [1, 64] {
                let ctx = format!("seed={seed} shards={shards} cfg={c:?}");
                let mut rng = Rng::new(seed);
                let (reads, dels) = workload(&mut rng, c.k);
                let mut a = build_one_writer(&c, &reads, &dels, shards);
                let mut b = build_one_writer(&c, &reads, &dels, shards);
                let ra = a.defrag();
                let rb = b.with_table(KmerTable::defrag);
                assert_eq!(ra.joins, rb.joins, "{ctx}");
                assert_eq!(ra.pieces_after, rb.pieces_after, "{ctx}");
                assert_eq!(ra.array_bases_after, rb.array_bases_after, "{ctx}");
                assert_eq!(a.edit_stats().dead_bases, 0, "{ctx}");
                let (ta, tb) = (a.into_table(), b.into_table());
                assert_eq!(ta.pieces().len(), tb.pieces().len(), "{ctx}");
                for (i, (pa, pb)) in ta.pieces().iter().zip(tb.pieces().iter()).enumerate() {
                    assert_eq!(pa, pb, "{ctx}: piece {i} differs");
                }
                let mut o = Oracle::new(c.k);
                for r in &reads {
                    o.insert_read(r);
                }
                for d in &dels {
                    o.delete_read(d);
                }
                check_table(ta, &o, &ctx);
            }
        }
    }
}

#[test]
fn native_defrag_incremental_keeps_ids_and_bases() {
    for seed in seeds(&[1411, 1412]) {
        for c in configs() {
            let ctx = format!("seed={seed} cfg={c:?}");
            let mut rng = Rng::new(seed);
            let (reads, dels) = workload(&mut rng, c.k);
            let mut a = build_one_writer(&c, &reads, &dels, 16);
            let mut b = build_one_writer(&c, &reads, &dels, 16);
            let dead = a.edit_stats().dead_bases;
            assert!(dead > 0, "{ctx}: nothing to reclaim");
            let ra = a.defrag_incremental();
            b.with_table(KmerTable::defrag_incremental);
            assert_eq!(a.edit_stats().dead_bases, 0, "{ctx}");
            assert_eq!(
                a.edit_stats().shared_bases,
                b.edit_stats().shared_bases,
                "{ctx}"
            );
            assert!(ra.array_bases_after < ra.array_bases_before, "{ctx}");
            assert_eq!(ra.pieces_after, ra.pieces_before, "{ctx}");
            let (ta, tb) = (a.into_table(), b.into_table());
            assert_eq!(ta.pieces().len(), tb.pieces().len(), "{ctx}");
            for (i, (pa, pb)) in ta.pieces().iter().zip(tb.pieces().iter()).enumerate() {
                assert_eq!(pa, pb, "{ctx}: piece {i} differs");
            }
            let mut o = Oracle::new(c.k);
            for r in &reads {
                o.insert_read(r);
            }
            for d in &dels {
                o.delete_read(d);
            }
            check_table(ta, &o, &ctx);
        }
    }
}

/// Defragment between concurrent phases (several writers and deleters):
/// queries are right after each step, and the table keeps working.
#[test]
fn defrag_between_concurrent_phases() {
    for seed in seeds(&[1421, 1422]) {
        for c in configs() {
            for full in [true, false] {
                let ctx = format!("seed={seed} full={full} cfg={c:?}");
                let k = c.k;
                let mut rng = Rng::new(seed);
                let (reads, dels) = workload(&mut rng, k);
                let mut o = Oracle::new(k);
                let mut t = ConcurrentKmerTable::with_shards(c.clone(), 32).unwrap();
                let mut seen: Vec<Vec<u8>> = Vec::new();
                for phase in 0..3 {
                    let part: Vec<Vec<u8>> = reads.iter().skip(phase).step_by(3).cloned().collect();
                    let dpart: Vec<Vec<u8>> = dels.iter().skip(phase).step_by(3).cloned().collect();
                    // Writers first, then deleters (each phase's outcome is
                    // then deterministic as a set).
                    std::thread::scope(|s| {
                        for chunk in part.chunks(part.len().div_ceil(3)) {
                            let t = &t;
                            s.spawn(move || {
                                let mut w = t.writer();
                                for r in chunk {
                                    w.insert_read(r);
                                }
                            });
                        }
                    });
                    std::thread::scope(|s| {
                        for chunk in dpart.chunks(dpart.len().div_ceil(2).max(1)) {
                            let t = &t;
                            s.spawn(move || {
                                for d in chunk {
                                    t.delete_read(d);
                                }
                            });
                        }
                    });
                    for r in &part {
                        o.insert_read(r);
                    }
                    for d in &dpart {
                        o.delete_read(d);
                    }
                    seen.extend(part);
                    let pctx = format!("{ctx} phase={phase}");
                    if full {
                        t.defrag();
                    } else {
                        t.defrag_incremental();
                    }
                    assert_eq!(t.edit_stats().dead_bases, 0, "{pctx}");
                    same_as_oracle(&t, &o, &seen, &mut rng, &pctx);
                    // Deleted k-mers come back when inserted again.
                    if let Some(d) = dels.get(phase) {
                        t.writer().insert_read(d);
                        o.insert_read(d);
                        assert!(
                            canonical_kmers(d, k).all(|x| t.contains(x)),
                            "{pctx}: re-inserted stretch"
                        );
                    }
                }
                check_table(t.into_table(), &o, &ctx);
            }
        }
    }
}
