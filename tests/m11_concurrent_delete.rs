//! Concurrent deletion (`ConcurrentKmerTable::delete_read` / `delete`),
//! alone and racing with writers, readers and other deleters, checked
//! against the HashSet oracle: exact final set, no k-mer stored twice, exact
//! index after conversion, no false negative for a k-mer that is never
//! deleted and no false positive at any moment.
//!
//! `UNIHASH_STRESS=N` repeats the multi-threaded tests N times.

mod common;
use common::seeds;
use std::collections::HashSet;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
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

fn kmers(reads: &[Vec<u8>], k: usize) -> HashSet<u64> {
    reads.iter().flat_map(|r| canonical_kmers(r, k)).collect()
}

/// A random substring of `r` (at least k bases when possible), maybe
/// reverse-complemented.
fn slice(rng: &mut Rng, r: &[u8], k: usize) -> Vec<u8> {
    if r.len() <= k {
        return r.to_vec();
    }
    let a = rng.below((r.len() - k) as u64) as usize;
    let b = (a + k + rng.below(120) as usize).min(r.len());
    let s = r[a..b].to_vec();
    if rng.chance(0.5) {
        revcomp_seq(&s)
    } else {
        s
    }
}

/// The table holds exactly `want` (k-mers in `ambiguous` may be either way),
/// each once, with an exact index.
fn check(t: &KmerTable, want: &HashSet<u64>, ambiguous: &HashSet<u64>, dups: bool, ctx: &str) {
    let mut seen = HashSet::new();
    let mut occ = 0u64;
    t.for_each_stored_kmer(|_, _, x| {
        assert!(
            want.contains(&x) || ambiguous.contains(&x),
            "{ctx}: stored k-mer should be absent"
        );
        seen.insert(x);
        occ += 1;
    });
    for x in want {
        if !ambiguous.contains(x) {
            assert!(t.contains(*x), "{ctx}: k-mer missing");
        }
    }
    if !dups {
        assert_eq!(occ, seen.len() as u64, "{ctx}: a k-mer is stored twice");
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

/// One thread: the concurrent table follows the oracle through random
/// insertions and deletions, in both run-indexing modes of the starting
/// table (deferred tables hold intra-run duplicates).
#[test]
fn deletion_matches_oracle_single_thread() {
    for seed in seeds(&[4001, 4002]) {
        for c in configs() {
            for mode in [RunIndexing::Eager, RunIndexing::Deferred] {
                let c = Config {
                    run_indexing: mode,
                    ..c.clone()
                };
                let ctx = format!("seed={seed} cfg={c:?}");
                let mut rng = Rng::new(seed);
                let genome = repeat_rich_genome(&mut rng, 15_000);
                let mut pool = sample_reads(&mut rng, &genome, 300, 30, 250, 0.01);
                pool.push(b"ACGTTGCA".repeat(30));
                pool.push(vec![b'A'; 150]);
                // Start from a KmerTable (deferred mode keeps duplicates).
                let mut st = KmerTable::new(c.clone()).unwrap();
                let mut o = Oracle::new(c.k);
                for r in &pool[..100] {
                    st.insert_read(r);
                    o.insert_read(r);
                }
                let ct = ConcurrentKmerTable::from_table_with_shards(st, 7);
                {
                    let mut w = ct.writer();
                    for step in 0..500 {
                        let r = &pool[rng.below(pool.len() as u64) as usize];
                        match rng.below(10) {
                            0..=3 => {
                                w.insert_read(r);
                                o.insert_read(r);
                            }
                            4..=7 => {
                                let s = slice(&mut rng, r, c.k);
                                ct.delete_read(&s);
                                o.delete_read(&s);
                            }
                            _ => {
                                let x = canonical_kmers(r, c.k).next();
                                if let Some(x) = x {
                                    assert_eq!(ct.delete(x), o.delete(x), "{ctx} step {step}");
                                }
                            }
                        }
                        if step % 100 == 99 {
                            for _ in 0..3000 {
                                let x = rng.next_u64() & kmer_mask(c.k);
                                assert_eq!(
                                    ct.contains(x),
                                    o.contains(x),
                                    "{ctx} step {step}: random"
                                );
                            }
                            for r in &pool {
                                for x in canonical_kmers(r, c.k) {
                                    assert_eq!(ct.contains(x), o.contains(x), "{ctx} step {step}");
                                }
                            }
                        }
                    }
                }
                let e = ct.edit_stats();
                assert!(e.cuts > 0 && e.splits > 0, "{ctx}: {e:?}");
                let want: HashSet<u64> = o.iter().copied().collect();
                let t = ct.into_table();
                check(
                    &t,
                    &want,
                    &HashSet::new(),
                    mode == RunIndexing::Deferred,
                    &ctx,
                );
            }
        }
    }
}

/// Deleters (overlapping stretches and single k-mers), writers (new reads)
/// and a reader run at the same time. A k-mer that is never deleted is
/// always found; nothing outside the inserted k-mers is ever found; the end
/// state is exact.
#[test]
fn deleters_writers_and_readers_race() {
    for _ in 0..stress() {
        for seed in seeds(&[4011, 4012]) {
            for c in configs() {
                for shards in [3, 64] {
                    let ctx = format!("seed={seed} shards={shards} cfg={c:?}");
                    let mut rng = Rng::new(seed);
                    let k = c.k;
                    let ga = repeat_rich_genome(&mut rng, 20_000);
                    let a_reads = sample_reads(&mut rng, &ga, 300, 40, 300, 0.01);
                    // Repeat-rich too: minimizers keep becoming hot while
                    // deletions run.
                    let gn = repeat_rich_genome(&mut rng, 20_000);
                    let n_reads = sample_reads(&mut rng, &gn, 200, 40, 300, 0.01);
                    // Stretches to delete: three overlapping lists.
                    let mut dels: Vec<Vec<Vec<u8>>> = vec![Vec::new(); 3];
                    for (i, r) in a_reads.iter().enumerate().filter(|(i, _)| i % 3 == 0) {
                        let s = slice(&mut rng, r, k);
                        dels[i % 3].push(s.clone());
                        dels[(i + 1) % 3].push(s);
                    }
                    let a_k = kmers(&a_reads, k);
                    let n_k = kmers(&n_reads, k);
                    let d_k: HashSet<u64> = dels.iter().flat_map(|v| kmers(v, k)).collect();
                    let ambiguous: HashSet<u64> = d_k.intersection(&n_k).copied().collect();
                    let mut want: HashSet<u64> = a_k.difference(&d_k).copied().collect();
                    want.extend(n_k.difference(&d_k));
                    let stable: Vec<u64> = a_k
                        .iter()
                        .filter(|x| !d_k.contains(x))
                        .copied()
                        .take(5000)
                        .collect();
                    let all: HashSet<u64> = a_k.union(&n_k).copied().collect();

                    let ct = ConcurrentKmerTable::with_shards(c.clone(), shards).unwrap();
                    std::thread::scope(|s| {
                        for part in a_reads.chunks(a_reads.len().div_ceil(3)) {
                            let ct = &ct;
                            s.spawn(move || {
                                let mut w = ct.writer();
                                for r in part {
                                    w.insert_read(r);
                                }
                            });
                        }
                    });
                    let busy = AtomicUsize::new(dels.len() + 2);
                    let done = AtomicBool::new(false);
                    std::thread::scope(|s| {
                        for (i, list) in dels.iter().enumerate() {
                            let (ct, busy) = (&ct, &busy);
                            s.spawn(move || {
                                for (j, d) in list.iter().enumerate() {
                                    if (i + j) % 4 == 0 {
                                        for x in canonical_kmers(d, k).collect::<Vec<_>>() {
                                            ct.delete(x);
                                        }
                                    } else {
                                        ct.delete_read(d);
                                    }
                                }
                                busy.fetch_sub(1, Ordering::Release);
                            });
                        }
                        for part in n_reads.chunks(n_reads.len().div_ceil(2)) {
                            let (ct, busy) = (&ct, &busy);
                            s.spawn(move || {
                                let mut w = ct.writer();
                                for r in part {
                                    w.insert_read(r);
                                }
                                drop(w);
                                busy.fetch_sub(1, Ordering::Release);
                            });
                        }
                        let (ct, busy, done, stable, all, ctx) =
                            (&ct, &busy, &done, &stable, &all, &ctx);
                        s.spawn(move || {
                            let mut rng = Rng::new(seed ^ 77);
                            while busy.load(Ordering::Acquire) > 0 {
                                for &x in stable.iter().step_by(7) {
                                    assert!(ct.contains(x), "{ctx}: k-mer never deleted not found");
                                }
                                for _ in 0..200 {
                                    let x = rng.next_u64() & kmer_mask(k);
                                    if ct.contains(x) {
                                        let c = unihash::codec::canonicalize(x, k);
                                        assert!(all.contains(&c), "{ctx}: false positive");
                                    }
                                }
                            }
                            done.store(true, Ordering::Release);
                        });
                    });
                    assert!(done.load(Ordering::Acquire));
                    let e = ct.edit_stats();
                    assert!(e.cuts > 0, "{ctx}");
                    let t = ct.into_table();
                    check(&t, &want, &ambiguous, false, &ctx);
                }
            }
        }
    }
}

/// A deletion that hits a writer's open run waits for the run to close.
#[test]
fn deletion_waits_for_open_runs() {
    for seed in seeds(&[4021]) {
        let c = Config::default();
        let k = c.k;
        let mut rng = Rng::new(seed);
        let read = rng.dna(3000);
        let stretch = read[1000..1200].to_vec();
        let ct = ConcurrentKmerTable::with_shards(c, 16).unwrap();
        let fed = AtomicBool::new(false);
        let closed = AtomicBool::new(false);
        std::thread::scope(|s| {
            let (ct, read, fed, closed) = (&ct, &read, &fed, &closed);
            s.spawn(move || {
                let mut w = ct.writer();
                w.insert_chunk(&read[..2000]);
                fed.store(true, Ordering::Release);
                std::thread::sleep(std::time::Duration::from_millis(200));
                w.insert_chunk(&read[2000..]);
                closed.store(true, Ordering::Release);
                w.end_read();
            });
            let stretch = &stretch;
            s.spawn(move || {
                while !fed.load(Ordering::Acquire) {
                    std::thread::yield_now();
                }
                let n = ct.delete_read(stretch);
                // The run was still open: this returns only after the
                // writer went on (its run closes at end_read).
                assert!(closed.load(Ordering::Acquire), "cut an open run");
                assert_eq!(n, (stretch.len() - k + 1) as u64);
            });
        });
        for x in canonical_kmers(&stretch, k) {
            assert!(!ct.contains(x));
        }
        let mut want = kmers(std::slice::from_ref(&read), k);
        for x in canonical_kmers(&stretch, k) {
            want.remove(&x);
        }
        let t = ct.into_table();
        check(&t, &want, &HashSet::new(), false, "open run");
    }
}

/// Defragmentation through the exclusive wrappers, between concurrent
/// phases of insertion and deletion.
#[test]
fn maintenance_between_concurrent_phases() {
    for seed in seeds(&[4031, 4032]) {
        for c in configs() {
            let ctx = format!("seed={seed} cfg={c:?}");
            let mut rng = Rng::new(seed);
            let genome = repeat_rich_genome(&mut rng, 20_000);
            let reads = sample_reads(&mut rng, &genome, 600, 40, 300, 0.01);
            let mut o = Oracle::new(c.k);
            let mut ct = ConcurrentKmerTable::with_shards(c.clone(), 32).unwrap();
            for phase in 0..3 {
                let part = &reads[phase * 200..(phase + 1) * 200];
                let dels: Vec<Vec<u8>> = reads[..(phase + 1) * 200]
                    .iter()
                    .step_by(5)
                    .map(|r| slice(&mut rng, r, c.k))
                    .collect();
                std::thread::scope(|s| {
                    for chunk in part.chunks(70) {
                        let ct = &ct;
                        s.spawn(move || {
                            let mut w = ct.writer();
                            for r in chunk {
                                w.insert_read(r);
                            }
                        });
                    }
                });
                for r in part {
                    o.insert_read(r);
                }
                std::thread::scope(|s| {
                    for chunk in dels.chunks(15) {
                        let ct = &ct;
                        s.spawn(move || {
                            for d in chunk {
                                ct.delete_read(d);
                            }
                        });
                    }
                });
                for d in &dels {
                    o.delete_read(d);
                }
                let e = ct.edit_stats();
                assert!(e.dead_bases > 0 || e.cuts == 0, "{ctx}");
                if phase == 0 {
                    ct.defrag_incremental();
                    assert_eq!(ct.edit_stats().dead_bases, 0, "{ctx}");
                } else {
                    ct.defrag();
                }
                for r in &reads[..(phase + 1) * 200] {
                    for x in canonical_kmers(r, c.k) {
                        assert_eq!(ct.contains(x), o.contains(x), "{ctx} phase {phase}");
                    }
                }
            }
            let want: HashSet<u64> = o.iter().copied().collect();
            let t = ct.into_table();
            check(&t, &want, &HashSet::new(), false, &ctx);
        }
    }
}
