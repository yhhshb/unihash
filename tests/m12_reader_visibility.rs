//! Visibility for readers of `ConcurrentKmerTable`: a k-mer that is present
//! stays visible while deletions move it to another piece id (every cut of a
//! piece's first k-mer files all the remaining ones under a new id) and while
//! writers make its minimizer hot (its entries move to layer 2 or 3).
//!
//! `UNIHASH_STRESS=N` repeats every test N times.

mod common;
use common::seeds;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use unihash::codec::{canonical_kmers, canonicalize, encode_kmer};
use unihash::rng::Rng;
use unihash::{CandidateScan, ConcurrentKmerTable, Config};

fn stress() -> usize {
    std::env::var("UNIHASH_STRESS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(1)
}

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
    ]
}

/// One deleter removes the k-mers of a piece from the front, one at a time:
/// every cut files all the remaining k-mers under a new piece id and empties
/// the old one. Readers check the k-mers not deleted yet the whole time.
#[test]
fn kmers_moved_by_cuts_stay_visible() {
    for _ in 0..stress() {
        for seed in seeds(&[5101, 5102]) {
            for c in configs() {
                for shards in [4, 256] {
                    let ctx = format!("seed={seed} shards={shards} cfg={c:?}");
                    let k = c.k;
                    let mut rng = Rng::new(seed);
                    let seq = rng.dna(1200);
                    let xs: Vec<u64> = (0..=seq.len() - k)
                        .map(|i| canonicalize(encode_kmer(&seq[i..i + k]).unwrap(), k))
                        .collect();
                    let ct = ConcurrentKmerTable::with_shards(c.clone(), shards).unwrap();
                    ct.writer().insert_read(&seq);
                    // K-mers 0..progress are deleted or being deleted.
                    let progress = AtomicUsize::new(0);
                    let done = AtomicBool::new(false);
                    std::thread::scope(|s| {
                        for _ in 0..2 {
                            let (ct, xs, progress, done, ctx) = (&ct, &xs, &progress, &done, &ctx);
                            s.spawn(move || {
                                while !done.load(Ordering::Acquire) {
                                    let p = progress.load(Ordering::Acquire);
                                    for (j, &x) in xs.iter().enumerate().skip(p + 1) {
                                        if !ct.contains(x) {
                                            let p2 = progress.load(Ordering::Acquire);
                                            assert!(
                                                j <= p2,
                                                "{ctx}: k-mer {j} not found (deleting {p2})"
                                            );
                                        }
                                    }
                                }
                            });
                        }
                        for (i, &x) in xs.iter().enumerate() {
                            progress.store(i, Ordering::Release);
                            assert!(ct.delete(x), "{ctx}: k-mer {i} not deleted");
                        }
                        done.store(true, Ordering::Release);
                    });
                    let st = ct.stats();
                    assert_eq!(st.live_bases, 0, "{ctx}");
                    assert_eq!(
                        (st.entries, st.entries2, st.layer3_kmers),
                        (0, 0, 0),
                        "{ctx}"
                    );
                }
            }
        }
    }
}

/// Writers insert mutated copies of a sequence inserted before: the new
/// k-mers share its minimizers, which keep becoming hot (their entries move
/// to layer 2, then 3) while readers check the sequence's k-mers.
#[test]
fn kmers_stay_visible_through_promotions() {
    for _ in 0..stress() {
        for seed in seeds(&[5201]) {
            for c in configs() {
                let ctx = format!("seed={seed} cfg={c:?}");
                let k = c.k;
                let mut rng = Rng::new(seed);
                let base = rng.dna(3000);
                let copies: Vec<Vec<u8>> = (0..300)
                    .map(|_| {
                        let mut r = base.clone();
                        for b in r.iter_mut() {
                            if rng.chance(0.02) {
                                *b = rng.base();
                            }
                        }
                        r
                    })
                    .collect();
                let ct = ConcurrentKmerTable::with_shards(c.clone(), 16).unwrap();
                ct.writer().insert_read(&base);
                let stable: Vec<u64> = canonical_kmers(&base, k).collect();
                let busy = AtomicUsize::new(2);
                std::thread::scope(|s| {
                    for part in copies.chunks(copies.len().div_ceil(2)) {
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
                    let (ct, busy, stable, ctx) = (&ct, &busy, &stable, &ctx);
                    s.spawn(move || {
                        while busy.load(Ordering::Acquire) > 0 {
                            for &x in stable {
                                assert!(ct.contains(x), "{ctx}: inserted k-mer not found");
                            }
                        }
                    });
                });
                if c.hot_threshold > 0 {
                    assert!(ct.stats().hot1 > 0, "{ctx}: no promotion happened");
                }
                for r in std::iter::once(&base).chain(&copies) {
                    for x in canonical_kmers(r, k) {
                        assert!(ct.contains(x), "{ctx}: final state");
                    }
                }
            }
        }
    }
}
