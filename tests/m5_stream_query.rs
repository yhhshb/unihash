//! In-order (streaming) queries must answer exactly like independent lookups.

mod common;
use common::seeds;
use unihash::codec::{canonical_kmers, revcomp_seq};
use unihash::rng::{repeat_rich_genome, sample_reads, Rng};
use unihash::{CandidateScan, Config, KmerTable, Oracle, StreamStats};

/// Stream `seq` through a fresh query in random chunks; check every answer
/// against the oracle and `contains`, and that the reported k-mers are the
/// read's canonical k-mers in order.
fn check_stream(
    t: &KmerTable,
    o: &Oracle,
    seqs: &[Vec<u8>],
    rng: &mut Rng,
    ctx: &str,
) -> StreamStats {
    let k = t.k();
    let mut sq = t.query_stream();
    for (si, s) in seqs.iter().enumerate() {
        let expected: Vec<u64> = canonical_kmers(s, k).collect();
        let mut got = Vec::with_capacity(expected.len());
        let mut i = 0;
        while i < s.len() {
            let n = rng.range(1, 80);
            let j = (i + n).min(s.len());
            sq.push(&s[i..j], |x, present| got.push((x, present)));
            i = j;
        }
        sq.end();
        assert_eq!(got.len(), expected.len(), "{ctx}: seq {si} k-mer count");
        for (p, (&(x, present), &e)) in got.iter().zip(&expected).enumerate() {
            assert_eq!(x, e, "{ctx}: seq {si} k-mer {p} identity");
            assert_eq!(
                present,
                o.contains(x),
                "{ctx}: seq {si} k-mer {p} vs oracle"
            );
            assert_eq!(
                present,
                t.contains(x),
                "{ctx}: seq {si} k-mer {p} vs contains"
            );
        }
    }
    *sq.stats()
}

fn query_sets(rng: &mut Rng, genome: &[u8], reads: &[Vec<u8>]) -> Vec<Vec<u8>> {
    let mut q: Vec<Vec<u8>> = reads.iter().take(60).cloned().collect();
    // reverse complements of inserted reads
    q.extend(reads.iter().skip(60).take(30).map(|r| revcomp_seq(r)));
    // reads with a higher error rate → alternating hits and misses
    q.extend(sample_reads(rng, genome, 60, 40, 300, 0.05));
    // long stretches of the genome (cross many piece boundaries)
    for _ in 0..5 {
        let s = rng.below((genome.len() - 2000) as u64) as usize;
        q.push(genome[s..s + 2000].to_vec());
    }
    // unrelated sequence and N runs
    q.push(rng.dna(1000));
    let mut with_n = genome[100..900].to_vec();
    with_n[300..305].fill(b'N');
    q.push(with_n);
    q
}

#[test]
fn streaming_matches_oracle() {
    for seed in seeds(&[801, 802, 803]) {
        let mut rng = Rng::new(seed);
        let genome = repeat_rich_genome(&mut rng, 20_000);
        let reads = sample_reads(&mut rng, &genome, 250, 40, 300, 0.01);
        let configs = [
            (31, 17, 8, 1 << 20, CandidateScan::TieOffsets),
            (31, 17, 8, 1 << 20, CandidateScan::Range),
            (31, 15, 0, 40, CandidateScan::TieOffsets), // heavy collisions + cap splits
            (32, 19, 6, 1 << 20, CandidateScan::TieOffsets),
            (32, 32, 2, 32, CandidateScan::Range),
            (21, 11, 6, 1 << 20, CandidateScan::TieOffsets),
        ];
        for &(k, l, bucket_bits, max_piece_len, scan) in &configs {
            let c = Config {
                k,
                l,
                bucket_bits,
                max_piece_len,
                scan,
                ..Config::default()
            };
            let mut t = KmerTable::new(c.clone()).unwrap();
            let mut o = Oracle::new(k);
            for r in &reads {
                t.insert_read(r);
                o.insert_read(r);
            }
            let qs = query_sets(&mut rng, &genome, &reads);
            let ctx = format!("seed={seed} cfg={c:?}");
            let st = check_stream(&t, &o, &qs, &mut rng, &ctx);
            // With one k-mer per piece (max_piece_len == k) nothing can extend.
            if max_piece_len > k as u32 {
                assert!(st.extensions > 0, "{ctx}: no extension used");
            }
            assert!(st.cache_reuses + st.extensions > 0);
            assert_eq!(
                st.queries,
                st.extensions + st.cache_reuses + st.bucket_scans + st.fallback_lookups,
                "{ctx}: every query is answered by exactly one path"
            );
        }
    }
}

#[test]
fn streaming_small_k_fuzz() {
    for seed in seeds(&[811, 812]) {
        let mut rng = Rng::new(seed);
        for &(k, l) in &[(5, 3), (4, 2), (3, 1), (2, 1), (6, 6), (7, 4)] {
            for scan in [CandidateScan::TieOffsets, CandidateScan::Range] {
                let c = Config {
                    k,
                    l,
                    bucket_bits: 2,
                    max_piece_len: k as u32 + 3,
                    scan,
                    ..Config::default()
                };
                let mut t = KmerTable::new(c.clone()).unwrap();
                let mut o = Oracle::new(k);
                let mut qs = Vec::new();
                for _ in 0..200 {
                    let len = rng.range(0, 30);
                    let r = match rng.below(3) {
                        0 => rng.dna_from(len, b"AT"),
                        1 => rng.dna_from(len, b"ACG"),
                        _ => rng.dna(len),
                    };
                    if rng.chance(0.5) {
                        t.insert_read(&r);
                        o.insert_read(&r);
                    }
                    qs.push(r);
                }
                let ctx = format!("seed={seed} cfg={c:?}");
                check_stream(&t, &o, &qs, &mut rng, &ctx);
            }
        }
    }
}

#[test]
fn streaming_reuses_work_on_long_hits() {
    let mut rng = Rng::new(821);
    let genome = rng.dna(50_000);
    // One piece for the whole genome (no cap split inside the query).
    let mut t = KmerTable::new(Config {
        bucket_bits: 12,
        max_piece_len: 1 << 20,
        ..Config::default()
    })
    .unwrap();
    t.insert_read(&genome);
    let mut sq = t.query_stream();
    let mut n = 0;
    sq.push(&genome[1000..11_000], |_, p| {
        assert!(p);
        n += 1;
    });
    let st = *sq.stats();
    assert_eq!(n, 10_000 - 30);
    // One index lookup for the first k-mer, extensions for all the others.
    assert_eq!(st.extensions, st.queries - 1);
    // Reverse complement: extension runs backwards along the piece.
    let mut sq = t.query_stream();
    sq.push(&revcomp_seq(&genome[2000..3000]), |_, p| assert!(p));
    assert_eq!(sq.stats().extensions, sq.stats().queries - 1);
    // Random sequence: misses; the bucket is scanned once per minimizer run.
    let mut sq = t.query_stream();
    sq.push(&rng.dna(10_000), |_, _| {});
    let st = sq.stats();
    assert!(st.cache_reuses > st.bucket_scans, "{st:?}");
}
