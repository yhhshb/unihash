//! Concurrent insertion and queries (`ConcurrentKmerTable`), checked against
//! the single-threaded `KmerTable` and the HashSet oracle: same k-mer set,
//! no k-mer stored twice whatever the interleaving, no false negative for a
//! completed insertion and no false positive at any time, exact index after
//! conversion, and round trips with deletion and defragmentation.
//!
//! `UNIHASH_STRESS=N` repeats the multi-threaded tests N times.

mod common;
use common::{check_against_oracle, seeds};
use std::collections::HashSet;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::Mutex;
use unihash::codec::{canonical_kmers, kmer_mask, revcomp_seq};
use unihash::rng::{repeat_rich_genome, sample_reads, Rng};
use unihash::{CandidateScan, ConcurrentKmerTable, Config, KmerTable, Oracle, RunIndexing};

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

fn stress() -> usize {
    std::env::var("UNIHASH_STRESS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(1)
}

fn test_reads(rng: &mut Rng, genome_len: usize, n: usize) -> Vec<Vec<u8>> {
    let genome = repeat_rich_genome(rng, genome_len);
    let mut reads = sample_reads(rng, &genome, n, 20, 400, 0.01);
    reads.push(vec![b'A'; 200]);
    reads.push(b"ACGTTGCA".repeat(40));
    reads.push(b"ACGTNNACGTACGTNACGT".repeat(10));
    let tandem = rng.dna(37).repeat(9);
    reads.push(revcomp_seq(&tandem));
    reads.push(tandem);
    reads.push(genome[..genome_len.min(5000)].to_vec());
    reads
}

/// Insert `r` through `w`, cut into random chunks half of the time.
fn insert(w: &mut unihash::Writer<'_>, r: &[u8], rng: &mut Rng) {
    if rng.chance(0.5) {
        let mut i = 0;
        while i < r.len() {
            let j = (i + rng.range(1, 50)).min(r.len());
            w.insert_chunk(&r[i..j]);
            i = j;
        }
        w.end_read();
    } else {
        w.insert_read(r);
    }
}

/// The converted table is exact, duplicate-free (unless it was partly built
/// with deferred run indexing) and its index exact.
fn check_table(t: &KmerTable, o: &Oracle, probes: &[Vec<u8>], rng: &mut Rng, ctx: &str) {
    check_table_dups(t, o, probes, rng, ctx, false)
}

fn check_table_dups(
    t: &KmerTable,
    o: &Oracle,
    probes: &[Vec<u8>],
    rng: &mut Rng,
    ctx: &str,
    dups_allowed: bool,
) {
    for r in probes {
        for x in canonical_kmers(r, t.k()) {
            assert_eq!(t.contains(x), o.contains(x), "{ctx}: probe k-mer");
        }
    }
    for &x in o.iter() {
        assert!(t.contains(x), "{ctx}: oracle k-mer missing");
    }
    for _ in 0..2000 {
        let x = rng.next_u64() & kmer_mask(t.k());
        assert_eq!(t.contains(x), o.contains(x), "{ctx}: random k-mer");
    }
    let mut set = HashSet::new();
    let mut occ = 0u64;
    t.for_each_stored_kmer(|_, _, x| {
        set.insert(x);
        occ += 1;
    });
    assert_eq!(set.len(), o.len(), "{ctx}: distinct stored k-mers");
    if !dups_allowed {
        assert_eq!(occ, o.len() as u64, "{ctx}: a k-mer is stored twice");
    }
    assert_eq!(t.stats().stored_kmers, occ, "{ctx}: stored_kmers counter");
    for (pid, p) in t.pieces().iter().enumerate() {
        assert!(
            p.is_empty() || (p.len() >= t.k() && p.len() <= t.config().max_piece_len as usize),
            "{ctx}: piece {pid} has {} bases",
            p.len()
        );
    }
    if let Err(e) = t.verify_index() {
        panic!("{ctx}: {e}");
    }
}

#[test]
fn one_writer_matches_single_threaded_table() {
    for seed in seeds(&[3001, 3002]) {
        for c in configs() {
            for shards in [1, 7, 1024] {
                let ctx = format!("seed={seed} shards={shards} cfg={c:?}");
                let mut rng = Rng::new(seed);
                let reads = test_reads(&mut rng, 30_000, 500);
                let ct = ConcurrentKmerTable::with_shards(c.clone(), shards).unwrap();
                let mut single = KmerTable::new(c.clone()).unwrap();
                let mut o = Oracle::new(c.k);
                {
                    let mut w = ct.writer();
                    for r in &reads {
                        insert(&mut w, r, &mut rng);
                        single.insert_read(r);
                        o.insert_read(r);
                    }
                }
                for r in &reads {
                    for x in canonical_kmers(r, c.k) {
                        assert!(ct.contains(x), "{ctx}: concurrent table misses a k-mer");
                    }
                }
                let st = ct.stats();
                assert_eq!(st.stored_kmers, o.len() as u64, "{ctx}");
                let a = single.insert_stats();
                assert_eq!(st.insert.kmers_processed, a.kmers_processed, "{ctx}");
                assert_eq!(st.insert.present_hits, a.present_hits, "{ctx}");
                let t = ct.into_table();
                check_table(&t, &o, &reads, &mut rng, &ctx);
                check_against_oracle(&t, &o, &reads, &mut rng, 1000, &ctx);
            }
        }
    }
}

/// Every thread inserts every read (in its own order and chunking), so the
/// same new k-mers race in all threads.
#[test]
fn many_writers_store_each_kmer_once() {
    for _ in 0..stress() {
        for seed in seeds(&[3011, 3012, 3013]) {
            for c in configs() {
                // shards == 0: sized with `with_capacity` (256 shards here,
                // pre-sized, so the indexes hardly grow).
                for (shards, threads) in [(1, 4), (3, 8), (64, 6), (0, 8)] {
                    let ctx = format!("seed={seed} shards={shards} threads={threads} cfg={c:?}");
                    let mut rng = Rng::new(seed);
                    let reads = test_reads(&mut rng, 20_000, 300);
                    let ct = if shards == 0 {
                        ConcurrentKmerTable::with_capacity(c.clone(), 40_000).unwrap()
                    } else {
                        ConcurrentKmerTable::with_shards(c.clone(), shards).unwrap()
                    };
                    std::thread::scope(|s| {
                        for th in 0..threads {
                            let (ct, reads) = (&ct, &reads);
                            s.spawn(move || {
                                let mut rng = Rng::new(seed * 100 + th as u64);
                                let mut order: Vec<usize> = (0..reads.len()).collect();
                                rng.shuffle(&mut order);
                                let mut w = ct.writer();
                                for i in order {
                                    insert(&mut w, &reads[i], &mut rng);
                                }
                            });
                        }
                    });
                    let mut o = Oracle::new(c.k);
                    for r in &reads {
                        o.insert_read(r);
                    }
                    for &x in o.iter() {
                        assert!(ct.contains(x), "{ctx}: missing before conversion");
                    }
                    let st = ct.stats();
                    assert_eq!(st.stored_kmers, o.len() as u64, "{ctx}: duplicates");
                    assert_eq!(st.insert.reads, (threads * reads.len()) as u64, "{ctx}");
                    let t = ct.into_table();
                    check_table(&t, &o, &reads, &mut rng, &ctx);
                }
            }
        }
    }
}

/// Readers query while writers insert: a read whose insertion has completed
/// is always found entirely, and nothing outside the final set is ever
/// reported.
#[test]
fn readers_during_insertion() {
    for _ in 0..stress() {
        for seed in seeds(&[3021, 3022]) {
            for c in configs() {
                let ctx = format!("seed={seed} cfg={c:?}");
                let mut rng = Rng::new(seed);
                let reads = test_reads(&mut rng, 30_000, 600);
                let mut all = Oracle::new(c.k);
                for r in &reads {
                    all.insert_read(r);
                }
                let ct = ConcurrentKmerTable::with_shards(c.clone(), 16).unwrap();
                let done: Vec<AtomicBool> = reads.iter().map(|_| AtomicBool::new(false)).collect();
                let next = AtomicUsize::new(0);
                let writers_left = AtomicUsize::new(3);
                let checked = Mutex::new(0u64);
                std::thread::scope(|s| {
                    for th in 0..3 {
                        let (ct, reads, done, next, left) =
                            (&ct, &reads, &done, &next, &writers_left);
                        s.spawn(move || {
                            let mut rng = Rng::new(seed + th);
                            let mut w = ct.writer();
                            loop {
                                let i = next.fetch_add(1, Ordering::Relaxed);
                                if i >= reads.len() {
                                    break;
                                }
                                insert(&mut w, &reads[i], &mut rng);
                                done[i].store(true, Ordering::Release);
                            }
                            drop(w);
                            left.fetch_sub(1, Ordering::Release);
                        });
                    }
                    for th in 0..2 {
                        let (ct, reads, done, left, all, checked, ctx) =
                            (&ct, &reads, &done, &writers_left, &all, &checked, &ctx);
                        s.spawn(move || {
                            let mut rng = Rng::new(seed * 7 + th);
                            let mut n = 0u64;
                            while left.load(Ordering::Acquire) > 0 {
                                let i = rng.below(reads.len() as u64) as usize;
                                if done[i].load(Ordering::Acquire) {
                                    for x in canonical_kmers(&reads[i], c.k) {
                                        assert!(ct.contains(x), "{ctx}: completed read not found");
                                        n += 1;
                                    }
                                }
                                for _ in 0..20 {
                                    let x = rng.next_u64() & kmer_mask(c.k);
                                    if ct.contains(x) {
                                        assert!(all.contains(x), "{ctx}: false positive");
                                    }
                                }
                            }
                            *checked.lock().unwrap() += n;
                        });
                    }
                });
                let t = ct.into_table();
                check_table(&t, &all, &reads, &mut rng, &ctx);
            }
        }
    }
}

/// `KmerTable` → concurrent → parallel inserts → `KmerTable` → deletions and
/// defragmentation → concurrent again → parallel inserts → `KmerTable`.
#[test]
fn round_trips_with_deletion_and_defrag() {
    for seed in seeds(&[3031, 3032]) {
        for c in configs() {
            for mode in [RunIndexing::Eager, RunIndexing::Deferred] {
                let c = Config {
                    run_indexing: mode,
                    ..c.clone()
                };
                let ctx = format!("seed={seed} cfg={c:?}");
                let mut rng = Rng::new(seed);
                let reads = test_reads(&mut rng, 25_000, 600);
                let (a, rest) = reads.split_at(200);
                let (b, d) = rest.split_at(200);
                let mut o = Oracle::new(c.k);
                let mut t = KmerTable::new(c.clone()).unwrap();
                for r in a {
                    t.insert_read(r);
                    o.insert_read(r);
                }
                let parallel = |t: KmerTable, part: &[Vec<u8>], shards: usize| {
                    let ct = ConcurrentKmerTable::from_table_with_shards(t, shards);
                    std::thread::scope(|s| {
                        for chunk in part.chunks(part.len().div_ceil(3)) {
                            let ct = &ct;
                            s.spawn(move || {
                                let mut w = ct.writer();
                                for r in chunk {
                                    w.insert_read(r);
                                }
                            });
                        }
                    });
                    ct.into_table()
                };
                t = parallel(t, b, 5);
                for r in b {
                    o.insert_read(r);
                }
                // Deferred tables may hold intra-run duplicates from part a.
                if mode == RunIndexing::Eager {
                    check_table(&t, &o, &reads, &mut rng, &format!("{ctx} after first"));
                }
                for r in reads.iter().step_by(3) {
                    let s = &r[r.len() / 4..r.len() * 3 / 4];
                    t.delete_read(s);
                    o.delete_read(s);
                }
                t.defrag_incremental();
                assert!(t.verify_index().is_ok(), "{ctx}");
                t.defrag();
                t = parallel(t, d, 64);
                for r in d {
                    o.insert_read(r);
                }
                t.defrag_incremental();
                check_table_dups(
                    &t,
                    &o,
                    &reads,
                    &mut rng,
                    &format!("{ctx} final"),
                    mode == RunIndexing::Deferred,
                );
                let st = t.stats();
                assert!(st.insert.reads >= reads.len() as u64, "{ctx}");
            }
        }
    }
}

/// The k-mers of a writer's open read are visible to other threads after
/// each `insert_chunk`.
#[test]
fn open_read_is_visible_after_each_chunk() {
    for seed in seeds(&[3041]) {
        for c in configs() {
            let ctx = format!("seed={seed} cfg={c:?}");
            let mut rng = Rng::new(seed);
            let reads = test_reads(&mut rng, 10_000, 60);
            let ct = ConcurrentKmerTable::with_shards(c.clone(), 4).unwrap();
            let mut w = ct.writer();
            for r in &reads {
                let mut i = 0;
                while i < r.len() {
                    let j = (i + rng.range(1, 60)).min(r.len());
                    w.insert_chunk(&r[i..j]);
                    let prefix = &r[..j];
                    std::thread::scope(|s| {
                        let (ct, ctx) = (&ct, &ctx);
                        s.spawn(move || {
                            for x in canonical_kmers(prefix, c.k) {
                                assert!(ct.contains(x), "{ctx}: open read k-mer not visible");
                            }
                        });
                    });
                    i = j;
                }
                w.end_read();
            }
        }
    }
}
