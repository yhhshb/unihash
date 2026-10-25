//! Milestone 4 — table insertion and query against the oracle.

mod common;
use common::{check_against_oracle, check_invariants, kmer, seeds};
use unihash::codec::{canonical_kmers, kmer_mask};
use unihash::rng::{repeat_rich_genome, sample_reads, Rng};
use unihash::{CandidateScan, Config, KmerTable, Oracle};

fn cfg(k: usize, l: usize) -> Config {
    Config {
        k,
        l,
        bucket_bits: 10,
        ..Config::default()
    }
}

/// Feed `read` through `insert_chunk` in random pieces (including empty and
/// 1-byte chunks), then `end_read`.
fn insert_in_random_chunks(t: &mut KmerTable, read: &[u8], rng: &mut Rng) {
    let mut i = 0;
    while i < read.len() {
        let n = match rng.below(4) {
            0 => 0,
            1 => 1,
            2 => rng.range(1, 7),
            _ => rng.range(1, 64),
        };
        let j = (i + n).min(read.len());
        t.insert_chunk(&read[i..j]);
        i = j;
    }
    t.end_read();
}

fn build(cfg: &Config, reads: &[Vec<u8>]) -> (KmerTable, Oracle) {
    let mut t = KmerTable::new(cfg.clone()).unwrap();
    let mut o = Oracle::new(cfg.k);
    for r in reads {
        t.insert_read(r);
        o.insert_read(r);
    }
    (t, o)
}

#[test]
fn hand_checked_overlapping_reads() {
    // k=5. r2 shares TCAGT and CAGTA with r1; its remaining k-mers are new
    // (none is the reverse complement of an r1 k-mer).
    let mut t = KmerTable::new(cfg(5, 3)).unwrap();
    t.insert_read(b"ACGGTCAGTA");
    t.insert_read(b"TCAGTACCAT");
    let p = t.pieces();
    assert_eq!(p.len(), 2);
    assert_eq!(p.piece_seq(0), b"ACGGTCAGTA");
    assert_eq!(p.piece_seq(1), b"AGTACCAT");
    let st = t.stats();
    assert_eq!(st.insert.present_hits, 2);
    assert_eq!(st.stored_kmers, 6 + 4);
    assert_eq!(t.redundant_kmers(), 0);
    for s in [
        "ACGGT", "CGGTC", "GGTCA", "GTCAG", "TCAGT", "CAGTA", "AGTAC", "GTACC", "TACCA", "ACCAT",
    ] {
        assert!(t.contains_seq(s.as_bytes()).unwrap(), "{s}");
    }
    // reverse complements are found too
    assert!(t.contains_seq(b"ATGGT").unwrap()); // rc(ACCAT)
    assert!(!t.contains_seq(b"AAAAA").unwrap());
    assert!(!t.contains_seq(b"CCATG").unwrap());
}

#[test]
fn hand_checked_reverse_complement_and_intra_run_repeat() {
    // r1 k-mers: AAAAA AAAAC AAACC AACCC ACCCC CCCCC.
    // r2 = CCCCCGGGGG: CCCCC present; CCCCG, CCCGG absent; CCGGG = rc(CCCGG)
    // and CGGGG = rc(CCCCG) are absent *when scanned* because the pending
    // run is not indexed yet (accepted redundancy); GGGGG = rc(CCCCC) is
    // present and closes the run.
    let mut t = KmerTable::new(cfg(5, 3)).unwrap();
    t.insert_read(b"AAAAACCCCC");
    t.insert_read(b"CCCCCGGGGG");
    assert_eq!(t.pieces().len(), 2);
    assert_eq!(t.pieces().piece_seq(1), b"CCCCGGGG");
    let st = t.stats();
    assert_eq!(st.insert.present_hits, 2);
    assert_eq!(st.stored_kmers, 10);
    assert_eq!(t.redundant_kmers(), 2);
    assert_eq!(t.distinct_kmers(), 8);
    for s in [
        "CCCCG", "CCCGG", "CCGGG", "CGGGG", "GGGGG", "TTTTT", "GGGTT",
    ] {
        assert!(t.contains_seq(s.as_bytes()).unwrap(), "{s}");
    }
    // a third insertion of the same read adds nothing
    t.insert_read(b"CCCCCGGGGG");
    assert_eq!(t.pieces().len(), 2);
}

#[test]
fn reads_shorter_than_k_produce_nothing() {
    let mut t = KmerTable::new(cfg(31, 17)).unwrap();
    t.insert_read(&[b'A'; 30]);
    t.insert_read(b"");
    t.insert_read(b"ACGTNACGTACGT");
    assert_eq!(t.pieces().len(), 0);
    assert_eq!(t.index().num_entries(), 0);
    assert_eq!(t.stats().insert.short_segments, 3);
}

fn oracle_configs() -> Vec<Config> {
    let mut v = Vec::new();
    for &(k, l) in &[(31, 17), (31, 15), (31, 19), (32, 19), (21, 11), (15, 7)] {
        for scan in [CandidateScan::TieOffsets, CandidateScan::Range] {
            v.push(Config {
                k,
                l,
                bucket_bits: 8,
                max_piece_len: 1 << 20,
                scan,
                ..Config::default()
            });
        }
    }
    // Force piece splitting and heavy bucket collisions.
    v.push(Config {
        k: 31,
        l: 17,
        bucket_bits: 0,
        max_piece_len: 40,
        scan: CandidateScan::Range,
        ..Config::default()
    });
    v.push(Config {
        k: 32,
        l: 32,
        bucket_bits: 2,
        max_piece_len: 32,
        scan: CandidateScan::TieOffsets,
        ..Config::default()
    });
    v
}

#[test]
fn property_against_oracle() {
    for seed in seeds(&[101, 102, 103]) {
        let mut rng = Rng::new(seed);
        let genome = repeat_rich_genome(&mut rng, 20_000);
        let mut reads = sample_reads(&mut rng, &genome, 250, 40, 300, 0.01);
        // a few exact duplicates and reverse complements of whole reads
        for i in 0..20 {
            reads.push(reads[i].clone());
            reads.push(unihash::codec::revcomp_seq(&reads[i + 20]));
        }
        for c in oracle_configs() {
            let (t, o) = build(&c, &reads);
            let ctx = format!("seed={seed} cfg={c:?}");
            let n = check_against_oracle(&t, &o, &reads, &mut rng, 20_000, &ctx);
            assert!(n > 0);
            check_invariants(&t, &o, &ctx);
        }
    }
}

#[test]
fn insertion_order_invariance() {
    for seed in seeds(&[201, 202]) {
        let mut rng = Rng::new(seed);
        let genome = repeat_rich_genome(&mut rng, 8_000);
        let reads = sample_reads(&mut rng, &genome, 200, 31, 150, 0.02);
        let c = Config {
            k: 31,
            l: 15,
            bucket_bits: 6,
            ..Config::default()
        };
        let (_, o) = build(&c, &reads);
        let mut sizes = Vec::new();
        for perm in 0..4 {
            let mut rs = reads.clone();
            rng.shuffle(&mut rs);
            let (t, o2) = build(&c, &rs);
            assert_eq!(o.len(), o2.len());
            let ctx = format!("seed={seed} perm={perm}");
            check_against_oracle(&t, &o, &reads, &mut rng, 5_000, &ctx);
            check_invariants(&t, &o, &ctx);
            sizes.push(t.distinct_kmers());
        }
        assert!(sizes.iter().all(|&s| s == o.len() as u64));
    }
}

/// Exhaustive check: every one of the 4^k k-mers is queried.
fn exhaustive(t: &KmerTable, o: &Oracle, ctx: &str) {
    let k = t.k();
    for x in 0..=kmer_mask(k) {
        assert_eq!(t.contains(x), o.contains(x), "{ctx}: k-mer {x:#x}");
    }
}

#[test]
fn small_k_fuzz() {
    let cases = [
        (5, 3),
        (5, 1),
        (5, 5),
        (4, 2),
        (3, 1),
        (2, 1),
        (2, 2),
        (6, 3),
        (7, 4),
    ];
    for seed in seeds(&[301, 302, 303, 304]) {
        let mut rng = Rng::new(seed);
        for &(k, l) in &cases {
            for chunked in [false, true] {
                for bucket_bits in [0u32, 3] {
                    for (mi, max_piece_len) in
                        [k as u32, k as u32 + 2, 1 << 20].into_iter().enumerate()
                    {
                        let scan = if mi % 2 == 0 {
                            CandidateScan::TieOffsets
                        } else {
                            CandidateScan::Range
                        };
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
                        let n_reads = rng.range(50, 400);
                        let mut reads = Vec::new();
                        for _ in 0..n_reads {
                            let len = rng.range(0, 25);
                            let mut r = match rng.below(4) {
                                0 => rng.dna_from(len, b"AT"),
                                1 => rng.dna_from(len, b"AC"),
                                _ => rng.dna(len),
                            };
                            if !r.is_empty() && rng.chance(0.1) {
                                let p = rng.below(r.len() as u64) as usize;
                                r[p] = b'N';
                            }
                            if chunked {
                                insert_in_random_chunks(&mut t, &r, &mut rng);
                            } else {
                                t.insert_read(&r);
                            }
                            o.insert_read(&r);
                            reads.push(r);
                        }
                        let ctx = format!("seed={seed} chunked={chunked} cfg={c:?}");
                        exhaustive(&t, &o, &ctx);
                        check_invariants(&t, &o, &ctx);
                    }
                }
            }
        }
    }
}

#[test]
fn k32_full_width() {
    for seed in seeds(&[401]) {
        let mut rng = Rng::new(seed);
        let genome = repeat_rich_genome(&mut rng, 10_000);
        let reads = sample_reads(&mut rng, &genome, 150, 32, 200, 0.01);
        for l in [8, 16, 19, 32] {
            let c = Config {
                k: 32,
                l,
                bucket_bits: 8,
                ..Config::default()
            };
            let (t, o) = build(&c, &reads);
            let ctx = format!("seed={seed} k=32 l={l}");
            check_against_oracle(&t, &o, &reads, &mut rng, 10_000, &ctx);
            check_invariants(&t, &o, &ctx);
            // all-T / all-A k-mers use every bit of the word
            assert_eq!(t.contains(u64::MAX), o.contains(u64::MAX));
            assert_eq!(t.contains(0), o.contains(0));
        }
    }
}

#[test]
fn long_read_split_at_cap() {
    let mut rng = Rng::new(501);
    let read = rng.dna(5_000);
    let c = Config {
        k: 31,
        l: 17,
        bucket_bits: 8,
        max_piece_len: 100,
        ..Config::default()
    };
    let (t, o) = build(&c, std::slice::from_ref(&read));
    let st = t.stats();
    assert!(st.insert.cap_splits > 40);
    // consecutive pieces overlap by exactly k-1 bases
    for i in 1..t.pieces().len() as u32 {
        let a = t.pieces().piece_seq(i - 1);
        let b = t.pieces().piece_seq(i);
        assert_eq!(&a[a.len() - 30..], &b[..30]);
    }
    check_against_oracle(&t, &o, &[read], &mut rng, 1_000, "cap");
    check_invariants(&t, &o, "cap");
}

#[test]
fn non_acgt_splits_reads() {
    let mut t = KmerTable::new(cfg(5, 3)).unwrap();
    t.insert_read(b"ACGTANNCCGTTRacgtt");
    // segments: ACGTA, CCGTT, acgtt
    assert_eq!(t.stats().insert.segments, 3);
    assert_eq!(t.pieces().len(), 3);
    assert!(t.contains(kmer("ACGTA")));
    assert!(t.contains(kmer("CCGTT")));
    assert!(t.contains(kmer("ACGTT")));
    // k-mers spanning the N/R characters are never stored
    assert!(!t.contains(kmer("CGTAC")));
    assert!(!t.contains(kmer("GTTAC")));
    let all: Vec<u64> = canonical_kmers(b"ACGTANNCCGTTRacgtt", 5).collect();
    assert_eq!(all.len(), 3);
    for x in all {
        assert!(t.contains(x));
    }
}

#[test]
fn chunked_insertion_is_identical_to_whole_reads() {
    for seed in seeds(&[601, 602, 603]) {
        let mut rng = Rng::new(seed);
        let genome = repeat_rich_genome(&mut rng, 30_000);
        let mut reads = sample_reads(&mut rng, &genome, 200, 20, 400, 0.01);
        // N runs at random places, including read ends.
        for r in reads.iter_mut() {
            for _ in 0..rng.below(3) {
                let p = rng.below(r.len() as u64) as usize;
                let n = rng.range(1, 5).min(r.len() - p);
                r[p..p + n].fill(b'N');
            }
        }
        reads.push(genome.clone()); // one long record
        for max_piece_len in [31u32, 64, 1 << 20] {
            let c = Config {
                k: 31,
                l: 17,
                bucket_bits: 8,
                max_piece_len,
                ..Config::default()
            };
            let mut whole = KmerTable::new(c.clone()).unwrap();
            let mut chunked = KmerTable::new(c.clone()).unwrap();
            for r in &reads {
                whole.insert_read(r);
                insert_in_random_chunks(&mut chunked, r, &mut rng);
            }
            let ctx = format!("seed={seed} max_piece_len={max_piece_len}");
            assert_eq!(whole.pieces().len(), chunked.pieces().len(), "{ctx}");
            for (a, b) in whole.pieces().iter().zip(chunked.pieces().iter()) {
                assert_eq!(a, b, "{ctx}: piece differs");
            }
            assert_eq!(whole.insert_stats(), chunked.insert_stats(), "{ctx}");
            let ea: Vec<_> = whole.index().iter_entries_with_bucket().collect();
            let eb: Vec<_> = chunked.index().iter_entries_with_bucket().collect();
            assert_eq!(ea, eb, "{ctx}: index differs");
        }
    }
}

#[test]
fn insertion_buffers_are_bounded_by_max_piece_len() {
    let mut rng = Rng::new(701);
    let max_piece_len = 4096u32;
    let c = Config {
        k: 31,
        l: 17,
        bucket_bits: 12,
        max_piece_len,
        ..Config::default()
    };
    let mut t = KmerTable::new(c).unwrap();
    let mut o = Oracle::new(31);
    // A 3 Mb "chromosome" streamed in 1000-byte chunks.
    let genome = rng.dna(3_000_000);
    for chunk in genome.chunks(1000) {
        t.insert_chunk(chunk);
        o.insert_chunk(chunk);
    }
    t.end_read();
    o.end_read();
    let bound = 2 * 14 * max_piece_len as usize + 1024;
    assert!(
        t.insert_buffer_bytes() <= bound,
        "insert buffers {} B exceed {bound} B",
        t.insert_buffer_bytes()
    );
    assert!(t.stats().insert.cap_splits > 500);
    assert_eq!(t.stats().insert.reads, 1);
    check_against_oracle(&t, &o, &[], &mut rng, 10_000, "bounded");
    check_invariants(&t, &o, "bounded");
}

#[test]
fn open_read_and_end_read_semantics() {
    let mut t = KmerTable::new(cfg(5, 3)).unwrap();
    t.end_read(); // no-op without an open read
    assert_eq!(t.stats().insert.reads, 0);
    t.insert_chunk(b"ACG");
    t.insert_chunk(b"TA"); // ACGTA spans the chunk boundary
                           // The pending run is not visible until it is emitted.
    assert!(!t.contains(kmer("ACGTA")));
    t.end_read();
    assert!(t.contains(kmer("ACGTA")));
    assert_eq!(t.stats().insert.reads, 1);
    // A new read does not continue the previous k-mer window.
    t.insert_chunk(b"CC");
    t.end_read();
    t.insert_chunk(b"GTT");
    t.end_read();
    assert!(!t.contains(kmer("ACCGT")));
    assert!(!t.contains(kmer("CCGTT")));
    assert_eq!(t.stats().insert.short_segments, 2);
}
