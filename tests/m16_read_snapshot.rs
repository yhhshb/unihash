//! `ReadSnapshot`: lookups and streaming queries through the snapshot (no
//! locks) give exactly the answers of the locked queries, from several
//! threads, on every index layer; writers wait until the snapshot is dropped.

mod common;
use common::seeds;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Duration;
use unihash::codec::{canonical_kmers, kmer_mask};
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

fn build(c: &Config, reads: &[Vec<u8>]) -> ConcurrentKmerTable {
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
    t
}

#[test]
fn snapshot_answers_equal_locked_answers() {
    for seed in seeds(&[1601, 1602]) {
        for c in configs() {
            let ctx = format!("seed={seed} cfg={c:?}");
            let k = c.k;
            let mut rng = Rng::new(seed);
            let g = repeat_rich_genome(&mut rng, 20_000);
            let reads = sample_reads(&mut rng, &g, 300, 40, 300, 0.01);
            let queries = sample_reads(&mut rng, &g, 100, 40, 300, 0.05);
            let t = build(&c, &reads);
            // Locked answers first, then the snapshot's from 3 threads.
            let mut kmers: Vec<u64> = reads
                .iter()
                .chain(&queries)
                .flat_map(|r| canonical_kmers(r, k).collect::<Vec<_>>())
                .collect();
            kmers.extend((0..2000).map(|_| rng.next_u64() & kmer_mask(k)));
            let locked: Vec<bool> = kmers.iter().map(|&x| t.contains(x)).collect();
            let locked_stream: Vec<Vec<(u64, bool)>> = queries
                .iter()
                .map(|r| {
                    let mut sq = t.query_stream();
                    let mut out = Vec::new();
                    sq.push(r, |x, p| out.push((x, p)));
                    out
                })
                .collect();
            assert!(
                locked.iter().any(|&p| p) && locked.iter().any(|&p| !p),
                "{ctx}"
            );
            let snap = t.read_snapshot();
            std::thread::scope(|s| {
                for part in 0..3 {
                    let (snap, kmers, locked, queries, locked_stream, ctx) =
                        (&snap, &kmers, &locked, &queries, &locked_stream, &ctx);
                    s.spawn(move || {
                        for (i, &x) in kmers.iter().enumerate().skip(part).step_by(3) {
                            assert_eq!(snap.contains(x), locked[i], "{ctx}: k-mer {i}");
                        }
                        for (i, r) in queries.iter().enumerate().skip(part).step_by(3) {
                            // Pushed in random chunks: the cache and the
                            // extension must not depend on chunk borders.
                            let mut rng = Rng::new(i as u64);
                            let mut sq = snap.query_stream();
                            let mut out = Vec::new();
                            let mut s = 0;
                            while s < r.len() {
                                let e = (s + rng.range(1, 120)).min(r.len());
                                sq.push(&r[s..e], |x, p| out.push((x, p)));
                                s = e;
                            }
                            assert_eq!(out, locked_stream[i], "{ctx}: streamed query {i}");
                        }
                    });
                }
            });
        }
    }
}

#[test]
fn writers_wait_for_the_snapshot() {
    let c = Config::default();
    let mut rng = Rng::new(1603);
    let g = repeat_rich_genome(&mut rng, 5_000);
    let reads = sample_reads(&mut rng, &g, 25, 60, 200, 0.0);
    let t = build(&c, &reads);
    let new_read: Vec<u8> = (0..200).map(|_| b"ACGT"[rng.below(4) as usize]).collect();
    let new_read = &new_read;
    let x = canonical_kmers(new_read, c.k)
        .find(|&x| !t.contains(x))
        .expect("a k-mer not inserted yet");
    let done = AtomicBool::new(false);
    std::thread::scope(|s| {
        let snap = t.read_snapshot();
        let (t, done) = (&t, &done);
        let h = s.spawn(move || {
            let mut w = t.writer();
            w.insert_read(new_read);
            done.store(true, Ordering::SeqCst);
        });
        std::thread::sleep(Duration::from_millis(200));
        assert!(
            !done.load(Ordering::SeqCst),
            "a writer ran while a snapshot was held"
        );
        assert!(
            !snap.contains(x),
            "the snapshot sees a k-mer inserted after it"
        );
        drop(snap);
        h.join().unwrap();
    });
    assert!(t.contains(x));
}
