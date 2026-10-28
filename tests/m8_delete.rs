//! Deletion (single k-mers and stretches), alternating insertions and
//! deletions, and defragmentation, all checked against the HashSet oracle
//! and against an exact recomputation of the index (`verify_index`).

mod common;
use common::seeds;
use std::collections::HashSet;
use unihash::codec::{canonical_kmers, kmer_mask, revcomp_seq};
use unihash::rng::{repeat_rich_genome, sample_reads, Rng};
use unihash::{CandidateScan, Config, KmerTable, Oracle, RunIndexing};

/// Answers, stored k-mer set and index layers all agree with the oracle.
fn check(t: &KmerTable, o: &Oracle, probes: &[Vec<u8>], rng: &mut Rng, ctx: &str) {
    let k = t.k();
    for (i, r) in probes.iter().enumerate() {
        for x in canonical_kmers(r, k) {
            assert_eq!(
                t.contains(x),
                o.contains(x),
                "{ctx}: probe {i} k-mer {}",
                unihash::codec::decode_kmer(x, k)
            );
        }
        // streaming answers too
        let mut sq = t.query_stream();
        let mut got = Vec::new();
        sq.push(r, |x, p| got.push((x, p)));
        sq.end();
        for (x, p) in got {
            assert_eq!(p, o.contains(x), "{ctx}: streaming probe {i}");
        }
    }
    for &x in o.iter() {
        assert!(t.contains(x), "{ctx}: oracle k-mer missing");
    }
    for _ in 0..500 {
        let x = rng.next_u64() & kmer_mask(k);
        assert_eq!(t.contains(x), o.contains(x), "{ctx}: random k-mer");
    }
    let mut stored = HashSet::new();
    let mut occ = 0u64;
    t.for_each_stored_kmer(|_, _, x| {
        stored.insert(x);
        occ += 1;
    });
    assert_eq!(stored.len(), o.len(), "{ctx}: stored distinct != oracle");
    assert!(
        stored.iter().all(|x| o.contains(*x)),
        "{ctx}: stored k-mer not in oracle"
    );
    let st = t.stats();
    assert_eq!(st.stored_kmers, occ, "{ctx}: stored_kmers counter");
    if t.config().run_indexing == RunIndexing::Eager {
        assert_eq!(
            occ,
            o.len() as u64,
            "{ctx}: eager run indexing stored a duplicate"
        );
    }
    for (pid, p) in t.pieces().iter().enumerate() {
        assert!(
            p.is_empty() || p.len() >= k,
            "{ctx}: piece {pid} has {} bases",
            p.len()
        );
        assert!(
            p.len() <= t.config().max_piece_len as usize,
            "{ctx}: piece too long"
        );
    }
    if let Err(e) = t.verify_index() {
        panic!("{ctx}: {e}");
    }
}

/// A random substring of `r` (at least k bases when possible), maybe
/// reverse-complemented.
fn slice(rng: &mut Rng, r: &[u8], k: usize) -> Vec<u8> {
    if r.len() <= k {
        return r.to_vec();
    }
    let a = rng.range(0, r.len() - k);
    let b = rng.range(a + k, r.len());
    let s = r[a..b].to_vec();
    if rng.chance(0.5) {
        revcomp_seq(&s)
    } else {
        s
    }
}

fn configs() -> Vec<Config> {
    vec![
        Config {
            bucket_bits: 4,
            ..Config::default()
        },
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
            k: 21,
            l: 11,
            bucket_bits: 2,
            hot_threshold: 3,
            max_piece_len: 60,
            run_indexing: RunIndexing::Deferred,
            ..Config::default()
        },
        Config {
            k: 32,
            l: 19,
            l2: 32,
            hot_threshold: 2,
            bucket_bits: 0,
            run_indexing: RunIndexing::Deferred,
            ..Config::default()
        },
    ]
}

#[test]
fn delete_inserted_reads_and_stretches() {
    for seed in seeds(&[1001, 1002, 1003]) {
        for c in configs() {
            let ctx = format!("seed={seed} cfg={c:?}");
            let mut rng = Rng::new(seed);
            let genome = repeat_rich_genome(&mut rng, 20_000);
            let mut reads = sample_reads(&mut rng, &genome, 200, 40, 300, 0.01);
            reads.push(vec![b'A'; 120]);
            reads.push(b"ACGTTGCA".repeat(30));
            let mut t = KmerTable::new(c.clone()).unwrap();
            let mut o = Oracle::new(c.k);
            for r in &reads {
                t.insert_read(r);
                o.insert_read(r);
            }
            check(&t, &o, &reads, &mut rng, &ctx);
            // delete whole reads, sub-stretches (either strand), single k-mers
            for step in 0..120 {
                let r = reads[rng.below(reads.len() as u64) as usize].clone();
                match step % 3 {
                    0 => {
                        let removed = t.delete_read(&r);
                        let present = o.delete_read(&r);
                        assert!(removed >= present, "{ctx}: removed {removed} < {present}");
                    }
                    1 => {
                        let s = slice(&mut rng, &r, c.k);
                        t.delete_read(&s);
                        o.delete_read(&s);
                    }
                    _ => {
                        let ks: Vec<u64> = canonical_kmers(&r, c.k).collect();
                        if ks.is_empty() {
                            continue;
                        }
                        let x = ks[rng.below(ks.len() as u64) as usize];
                        assert_eq!(t.delete(x), o.delete(x), "{ctx}: delete single");
                        assert!(!t.contains(x));
                    }
                }
                if step % 20 == 19 {
                    check(&t, &o, &reads, &mut rng, &format!("{ctx} step {step}"));
                }
            }
            check(&t, &o, &reads, &mut rng, &ctx);
            // deleting everything empties the table
            for r in &reads {
                t.delete_read(r);
                o.delete_read(r);
            }
            assert!(o.is_empty());
            check(&t, &o, &reads, &mut rng, &ctx);
            assert_eq!(t.stats().stored_kmers, 0, "{ctx}");
            assert_eq!(
                t.index().num_entries() + t.index2().num_entries(),
                0,
                "{ctx}"
            );
            assert_eq!(t.layer3_len(), 0, "{ctx}");
            let rep = t.defrag();
            assert_eq!(rep.pieces_after, 0, "{ctx}");
            check(&t, &o, &reads, &mut rng, &ctx);
        }
    }
}

#[test]
fn alternating_insertions_deletions_and_defrag() {
    for seed in seeds(&[1011, 1012, 1013]) {
        for c in configs() {
            let ctx = format!("seed={seed} cfg={c:?}");
            let mut rng = Rng::new(seed);
            let genome = repeat_rich_genome(&mut rng, 15_000);
            let pool = sample_reads(&mut rng, &genome, 400, 30, 250, 0.01);
            let mut t = KmerTable::new(c.clone()).unwrap();
            let mut o = Oracle::new(c.k);
            for step in 0..600 {
                let r = &pool[rng.below(pool.len() as u64) as usize];
                match rng.below(10) {
                    0..=4 => {
                        // insert, sometimes chunked
                        if rng.chance(0.3) {
                            let mut i = 0;
                            while i < r.len() {
                                let j = (i + rng.range(1, 40)).min(r.len());
                                t.insert_chunk(&r[i..j]);
                                i = j;
                            }
                            t.end_read();
                        } else {
                            t.insert_read(r);
                        }
                        o.insert_read(r);
                    }
                    5..=7 => {
                        let s = slice(&mut rng, r, c.k);
                        t.delete_read(&s);
                        o.delete_read(&s);
                    }
                    8 => {
                        // delete then re-insert the same stretch (piece churn)
                        let s = slice(&mut rng, r, c.k);
                        t.delete_read(&s);
                        o.delete_read(&s);
                        t.insert_read(&s);
                        o.insert_read(&s);
                    }
                    _ => {
                        if rng.chance(0.5) {
                            // incremental: only around pieces changed since
                            let lim = [0, 100, usize::MAX][rng.below(3) as usize];
                            let rep = t.defrag_incremental_with(lim);
                            let after = t.fragmentation();
                            assert_eq!(after.dead_bases, 0, "{ctx}");
                            // No dead bases left: the array is the union of
                            // the piece views.
                            assert!(t.pieces().total_bases() <= after.live_bases, "{ctx}");
                            assert_eq!(
                                rep.pieces_before - rep.joins,
                                rep.pieces_after,
                                "{ctx}: {rep}"
                            );
                            if let Err(e) = t.verify_index() {
                                panic!("{ctx} step {step} after incremental defrag: {e}");
                            }
                        } else if rng.chance(0.3) {
                            let before = t.fragmentation();
                            let rep = t.defrag();
                            let after = t.fragmentation();
                            assert_eq!(after.dead_bases, 0, "{ctx}");
                            assert_eq!(after.empty_pieces, 0, "{ctx}");
                            assert!(rep.pieces_after <= before.pieces, "{ctx}: {rep}");
                            assert_eq!(
                                rep.pieces_before - rep.joins,
                                rep.pieces_after,
                                "{ctx}: {rep}"
                            );
                        }
                    }
                }
                if step % 50 == 49 {
                    check(&t, &o, &pool, &mut rng, &format!("{ctx} step {step}"));
                }
            }
            check(&t, &o, &pool, &mut rng, &ctx);
            t.defrag_incremental();
            check(
                &t,
                &o,
                &pool,
                &mut rng,
                &format!("{ctx} final incremental defrag"),
            );
            let e = t.stats().edits;
            assert!(
                e.splits > 0
                    && e.cuts > e.splits
                    && e.defrags > 0
                    && e.incremental_defrags > 0
                    && e.joins > 0,
                "{ctx}: {e:?}"
            );
            t.defrag();
            check(&t, &o, &pool, &mut rng, &format!("{ctx} final defrag"));
            // A second defrag right away finds nothing new to join.
            let again = t.defrag();
            assert_eq!(again.joins, 0, "{ctx}: {again}");
        }
    }
}

#[test]
fn small_k_exhaustive_with_edits() {
    for seed in seeds(&[1021, 1022]) {
        let mut rng = Rng::new(seed);
        for &(k, l, l2, thr) in &[
            (5, 2, 3, 2),
            (6, 3, 5, 1),
            (7, 3, 0, 0),
            (4, 1, 4, 3),
            (3, 2, 0, 64),
        ] {
            let c = Config {
                k,
                l,
                l2,
                hot_threshold: thr,
                bucket_bits: 0,
                max_piece_len: k as u32 + 5,
                ..Config::default()
            };
            let ctx = format!("seed={seed} cfg={c:?}");
            let mut t = KmerTable::new(c.clone()).unwrap();
            let mut o = Oracle::new(k);
            for step in 0..400 {
                let len = rng.range(0, 25);
                let s = match rng.below(3) {
                    0 => rng.dna_from(len, b"AT"),
                    1 => rng.dna_from(len, b"ACG"),
                    _ => rng.dna(len),
                };
                match rng.below(7) {
                    0..=2 => {
                        t.insert_read(&s);
                        o.insert_read(&s);
                    }
                    3..=5 => {
                        t.delete_read(&s);
                        o.delete_read(&s);
                    }
                    _ => {
                        if rng.chance(0.25) {
                            t.defrag();
                        } else if rng.chance(0.5) {
                            t.defrag_incremental_with([0, 100, usize::MAX][rng.below(3) as usize]);
                        }
                    }
                }
                if step % 40 == 39 {
                    for x in 0..=kmer_mask(k) {
                        assert_eq!(
                            t.contains(x),
                            o.contains(x),
                            "{ctx} step {step}: k-mer {x:#x}"
                        );
                    }
                    if let Err(e) = t.verify_index() {
                        panic!("{ctx} step {step}: {e}");
                    }
                }
            }
        }
    }
}

#[test]
fn deletion_cuts_keep_flanks_and_share_bases() {
    // One piece, delete a stretch from its middle: two pieces that keep the
    // k-1 flanking bases, no copy of the remaining bases.
    let c = Config {
        k: 11,
        l: 5,
        ..Config::default()
    };
    let mut rng = Rng::new(7);
    let seq = rng.dna(200);
    let mut t = KmerTable::new(c).unwrap();
    t.insert_read(&seq);
    assert_eq!(t.pieces().len(), 1);
    let bases = t.pieces().total_bases();
    // delete k-mers 50..=89 (a 50-base stretch holds 40 k-mers)
    let removed = t.delete_read(&seq[50..100]);
    assert_eq!(removed, 40);
    let ps: Vec<Vec<u8>> = t.pieces().iter().filter(|p| !p.is_empty()).collect();
    assert_eq!(ps.len(), 2);
    assert_eq!(ps[0], seq[..50 + 10].to_vec()); // k-mers 0..50 + flank
    assert_eq!(ps[1], seq[90..].to_vec()); // k-mers 90.. (flank included)
    assert_eq!(t.pieces().total_bases(), bases, "no copy on deletion");
    let f = t.fragmentation();
    assert_eq!(f.dead_bases, 90 - 60);
    // re-inserting the stretch bridges the two parts; defrag merges all three
    t.insert_read(&seq[50..100]);
    let rep = t.defrag();
    assert_eq!(rep.joins, 2, "{rep}");
    assert_eq!(rep.pieces_after, 1);
    assert_eq!(t.pieces().piece_seq(0), seq);
    // reverse-complement deletion of a single k-mer in the middle leaves two
    // parts sharing k-2 bases
    let rc = revcomp_seq(&seq[120..131]);
    assert_eq!(t.delete_read(&rc), 1);
    let f = t.fragmentation();
    assert_eq!((f.dead_bases, f.shared_bases), (0, 9));
    t.verify_index().unwrap();
}

#[test]
fn defrag_merges_overlapping_read_pieces() {
    // Overlapping reads leave pieces that extend each other by k-1 bases.
    let mut rng = Rng::new(42);
    let genome = rng.dna(30_000);
    let reads = sample_reads(&mut rng, &genome, 3_000, 100, 150, 0.0);
    let mut t = KmerTable::new(Config::default()).unwrap();
    let mut o = Oracle::new(31);
    for r in &reads {
        t.insert_read(r);
        o.insert_read(r);
    }
    let before = t.stats();
    let rep = t.defrag();
    let after = t.stats();
    assert!(rep.joins > 0);
    assert!(after.pieces * 4 < before.pieces, "{rep}");
    assert!(rep.table_bytes_after < rep.table_bytes_before, "{rep}");
    let mut r2 = Rng::new(1);
    check(&t, &o, &reads, &mut r2, "reads defrag");
}

#[test]
#[should_panic(expected = "a read is open")]
fn deletion_requires_closed_read() {
    let mut t = KmerTable::new(Config::default()).unwrap();
    t.insert_chunk(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT");
    t.delete_read(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT");
}

#[test]
fn incremental_defrag_compacts_and_merges() {
    for seed in seeds(&[1031, 1032, 1033]) {
        let c = Config {
            k: 15,
            l: 7,
            max_piece_len: 600,
            ..Config::default()
        };
        let ctx = format!("seed={seed}");
        let mut rng = Rng::new(seed);
        let genome = rng.dna(20_000);
        let mut t = KmerTable::new(c).unwrap();
        let mut o = Oracle::new(15);
        t.insert_read(&genome);
        o.insert_read(&genome);
        t.defrag();
        let start = t.fragmentation();
        // Delete and re-insert stretches (either strand): every cut is
        // bridged again.
        let mut cut = Vec::new();
        for _ in 0..60 {
            let a = rng.range(0, genome.len() - 200);
            let b = a + rng.range(15, 200);
            let s = genome[a..b].to_vec();
            let s = if rng.chance(0.5) { revcomp_seq(&s) } else { s };
            t.delete_read(&s);
            o.delete_read(&s);
            cut.push(s);
        }
        for s in &cut {
            t.insert_read(s);
            o.insert_read(s);
        }
        let mid = t.fragmentation();
        assert!(
            mid.pieces > start.pieces + 60 && mid.dead_bases > 0,
            "{ctx}: {mid:?}"
        );
        let mid_bytes = t.table_bytes();
        let mut t_all = t.clone();
        // Default: the base array is compacted; long pieces are not merged.
        t.defrag_incremental();
        let after = t.fragmentation();
        assert_eq!(after.dead_bases, 0, "{ctx}");
        assert!(t.table_bytes() < mid_bytes, "{ctx}");
        let mut r2 = Rng::new(seed);
        check(&t, &o, std::slice::from_ref(&genome), &mut r2, &ctx);
        // Without the length limit every bridge is merged back (chains are
        // capped at max_piece_len and may break elsewhere than the full
        // defrag did).
        let rep = t_all.defrag_incremental_with(usize::MAX);
        let after = t_all.fragmentation();
        assert!(rep.joins > 0, "{ctx}: {rep}");
        assert!(
            (after.pieces - start.pieces) * 4 <= mid.pieces - start.pieces,
            "{ctx}: start {} mid {} after {}",
            start.pieces,
            mid.pieces,
            after.pieces
        );
        check(&t_all, &o, std::slice::from_ref(&genome), &mut r2, &ctx);
        t.defrag();
        check(&t, &o, &[genome], &mut r2, &ctx);
    }
}

#[test]
fn incremental_defrag_merges_short_read_pieces() {
    let mut rng = Rng::new(43);
    let genome = rng.dna(30_000);
    let reads = sample_reads(&mut rng, &genome, 3_000, 100, 150, 0.0);
    let mut t = KmerTable::new(Config::default()).unwrap();
    let mut o = Oracle::new(31);
    for r in &reads {
        t.insert_read(r);
        o.insert_read(r);
    }
    let before = t.stats();
    let rep = t.defrag_incremental_with(100);
    let after = t.stats();
    let live = |s: &unihash::TableStats| s.pieces - s.empty_pieces;
    assert!(rep.joins > 0 && live(&after) * 2 < live(&before), "{rep}");
    assert!(rep.table_bytes_after < rep.table_bytes_before, "{rep}");
    let mut r2 = Rng::new(2);
    check(&t, &o, &reads, &mut r2, "reads incremental defrag");
}
