//! Streaming queries on `ConcurrentKmerTable`: every answer equals
//! `contains` (on a quiet table), stays exact while the table changes
//! between and during pushes, and the callback may write to the table.
//!
//! `UNIHASH_STRESS=N` repeats the multi-threaded test N times.

mod common;
use common::seeds;
use std::collections::HashSet;
use std::sync::atomic::{AtomicUsize, Ordering};
use unihash::codec::{canonical_kmers, canonicalize, encode_kmer};
use unihash::rng::{repeat_rich_genome, sample_reads, Rng};
use unihash::{CandidateScan, ConcurrentKmerTable, Config};

fn configs() -> Vec<Config> {
    vec![
        Config::default(),
        Config {
            k: 21,
            l: 11,
            bucket_bits: 0,
            hot_threshold: 3,
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
            l: 7,
            hot_threshold: 0,
            max_piece_len: 40,
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

/// Streamed answers of `r` (pushed in random chunks), in order.
fn stream(t: &ConcurrentKmerTable, r: &[u8], rng: &mut Rng) -> Vec<(u64, bool)> {
    let mut sq = t.query_stream();
    let mut out = Vec::new();
    let mut s = 0;
    while s < r.len() {
        let e = (s + rng.range(1, 120)).min(r.len());
        sq.push(&r[s..e], |x, p| out.push((x, p)));
        s = e;
    }
    out
}

#[test]
fn answers_equal_contains() {
    for seed in seeds(&[1501, 1502]) {
        for c in configs() {
            let ctx = format!("seed={seed} cfg={c:?}");
            let k = c.k;
            let mut rng = Rng::new(seed);
            let g = repeat_rich_genome(&mut rng, 20_000);
            let reads = sample_reads(&mut rng, &g, 300, 40, 300, 0.01);
            let t = ConcurrentKmerTable::with_shards(c.clone(), 16).unwrap();
            std::thread::scope(|s| {
                for part in reads.chunks(reads.len().div_ceil(3)) {
                    let t = &t;
                    s.spawn(move || {
                        let mut w = t.writer();
                        for r in part {
                            w.insert_read(r);
                        }
                    });
                }
            });
            // Inserted reads, held-out reads with errors, N runs.
            let mut queries = reads.clone();
            queries.extend(sample_reads(&mut rng, &g, 100, 40, 300, 0.03));
            for q in queries.iter_mut().skip(250) {
                let p = rng.below(q.len() as u64) as usize;
                q[p] = b'N';
            }
            let mut sq = t.query_stream();
            for q in &queries {
                let want: Vec<(u64, bool)> =
                    canonical_kmers(q, k).map(|x| (x, t.contains(x))).collect();
                let got = stream(&t, q, &mut rng);
                assert_eq!(got, want, "{ctx}");
                // One query object for all sequences, too.
                let mut got = Vec::new();
                sq.push(q, |x, p| got.push((x, p)));
                sq.end();
                assert_eq!(got, want, "{ctx}");
            }
            let s = sq.stats();
            assert!(s.extensions > 0, "{ctx}: no extension");
            if c.hot_threshold != 1 {
                assert!(s.cache_reuses > 0, "{ctx}: no cache reuse");
            }
        }
    }
}

/// The same query object, its cursor and candidates kept across pushes,
/// while the table changes between pushes.
#[test]
fn answers_stay_exact_across_changes_between_pushes() {
    for seed in seeds(&[1511]) {
        for c in configs() {
            let k = c.k;
            let mut rng = Rng::new(seed);
            for trial in 0..60 {
                let ctx = format!("seed={seed} trial={trial} cfg={c:?}");
                let mut t = ConcurrentKmerTable::with_shards(c.clone(), 8).unwrap();
                let r = rng.dna(240);
                let mut sq = t.query_stream();
                // Absent, then inserted between two pushes of the same read.
                let mut got = Vec::new();
                sq.push(&r[..120], |_, p| got.push(p));
                assert!(got.iter().all(|&p| !p), "{ctx}");
                t.writer().insert_read(&r);
                got.clear();
                sq.push(&r[120..], |_, p| got.push(p));
                assert!(got.iter().all(|&p| p), "{ctx}: inserted k-mers not seen");
                sq.end();
                // Present, then the next k-mer deleted between two pushes.
                sq.push(&r[..150], |_, _| {});
                let next = canonicalize(encode_kmer(&r[151 - k..151]).unwrap(), k);
                assert!(t.delete(next), "{ctx}");
                let mut first = None;
                sq.push(&r[150..], |x, p| {
                    first.get_or_insert((x, p));
                });
                assert_eq!(first, Some((next, false)), "{ctx}: deleted k-mer seen");
                sq.end();
                drop(sq);
                // And after defragmentation.
                if trial % 2 == 0 {
                    t.defrag();
                } else {
                    t.defrag_incremental();
                }
                let want: Vec<(u64, bool)> =
                    canonical_kmers(&r, k).map(|x| (x, t.contains(x))).collect();
                assert_eq!(stream(&t, &r, &mut rng), want, "{ctx}");
            }
        }
    }
}

/// Streaming readers while writers insert and deleters delete: a k-mer
/// inserted before and never deleted is always found; nothing outside the
/// inserted k-mers is ever found.
#[test]
fn streams_during_insertion_and_deletion() {
    for _ in 0..stress() {
        for seed in seeds(&[1521]) {
            for c in configs() {
                let ctx = format!("seed={seed} cfg={c:?}");
                let k = c.k;
                let mut rng = Rng::new(seed);
                let g = repeat_rich_genome(&mut rng, 20_000);
                let a = sample_reads(&mut rng, &g, 200, 40, 300, 0.01);
                let gb = repeat_rich_genome(&mut rng, 60_000);
                let b = sample_reads(&mut rng, &gb, 3000, 40, 300, 0.01);
                // Delete every other read of a; the rest is stable.
                let dels: Vec<&Vec<u8>> = a.iter().step_by(2).collect();
                let del_k: HashSet<u64> = dels.iter().flat_map(|r| canonical_kmers(r, k)).collect();
                let b_k: HashSet<u64> = b.iter().flat_map(|r| canonical_kmers(r, k)).collect();
                let all: HashSet<u64> = a
                    .iter()
                    .chain(&b)
                    .flat_map(|r| canonical_kmers(r, k))
                    .collect();
                let t = ConcurrentKmerTable::with_shards(c.clone(), 16).unwrap();
                {
                    let mut w = t.writer();
                    for r in &a {
                        w.insert_read(r);
                    }
                }
                let busy = AtomicUsize::new(3);
                std::thread::scope(|s| {
                    for part in b.chunks(b.len().div_ceil(2)) {
                        let (t, busy) = (&t, &busy);
                        s.spawn(move || {
                            let mut w = t.writer();
                            for r in part {
                                w.insert_read(r);
                            }
                            drop(w);
                            busy.fetch_sub(1, Ordering::Release);
                        });
                    }
                    {
                        let (t, busy, dels) = (&t, &busy, &dels);
                        s.spawn(move || {
                            for d in dels {
                                t.delete_read(d);
                            }
                            busy.fetch_sub(1, Ordering::Release);
                        });
                    }
                    for i in 0..2 {
                        let (t, busy, a, all, del_k, b_k, ctx) =
                            (&t, &busy, &a, &all, &del_k, &b_k, &ctx);
                        s.spawn(move || {
                            let mut rng = Rng::new(seed + i);
                            let mut sq = t.query_stream();
                            let mut rounds = 0;
                            while busy.load(Ordering::Acquire) > 0 || rounds < 20 {
                                rounds += 1;
                                let r = &a[rng.below(a.len() as u64) as usize];
                                let mut s0 = 0;
                                while s0 < r.len() {
                                    let e = (s0 + rng.range(1, 100)).min(r.len());
                                    sq.push(&r[s0..e], |x, p| {
                                        if !del_k.contains(&x) && !b_k.contains(&x) {
                                            assert!(p, "{ctx}: stable k-mer not found");
                                        }
                                        if p {
                                            assert!(all.contains(&x), "{ctx}: false positive");
                                        }
                                    });
                                    s0 = e;
                                }
                                sq.end();
                                let q = rng.dna(150);
                                sq.push(&q, |x, p| {
                                    assert!(!p || all.contains(&x), "{ctx}: false positive");
                                });
                                sq.end();
                            }
                        });
                    }
                });
            }
        }
    }
}

/// The callback runs with no lock held: it may insert into the table.
#[test]
fn callback_may_write() {
    let c = Config::default();
    let k = c.k;
    let mut rng = Rng::new(1531);
    let t = ConcurrentKmerTable::with_shards(c, 4).unwrap();
    let mut w = t.writer();
    let mut sq = t.query_stream();
    for _ in 0..50 {
        let r = rng.dna(200);
        let mut missing = Vec::new();
        sq.push(&r, |x, p| {
            if !p {
                missing.push(x);
                w.insert_read(&r[..k]); // writes during the callback
            }
        });
        sq.end();
        w.insert_read(&r);
        assert!(canonical_kmers(&r, k).all(|x| t.contains(x)));
    }
}
