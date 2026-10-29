//! Insertion shortcuts: a k-mer right after a present one is checked by
//! extending the previous hit by one base, and the candidates of a
//! minimizer are reused while successive k-mers share it. The table must be
//! exactly the oracle's set, with every k-mer stored once (eager run
//! indexing), whatever the configuration, the chunking, and the edits made
//! between reads.

mod common;
use common::{check_against_oracle, check_invariants, seeds};
use unihash::rng::{repeat_rich_genome, sample_reads, Rng};
use unihash::{CandidateScan, Config, KmerTable, Oracle, RunIndexing};

fn configs() -> Vec<Config> {
    let base = vec![
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
            bucket_bits: 2,
            hot_threshold: 0, // routing off: long chains
            max_piece_len: 40,
            ..Config::default()
        },
    ];
    let mut v = Vec::new();
    for c in base {
        for ri in [RunIndexing::Eager, RunIndexing::Deferred] {
            v.push(Config {
                run_indexing: ri,
                ..c.clone()
            });
        }
    }
    v
}

/// Reads with repeats inside runs (tandem repeats, homopolymers) as well as
/// sampled reads with errors.
fn reads(rng: &mut Rng) -> Vec<Vec<u8>> {
    let g = repeat_rich_genome(rng, 20_000);
    let mut r = sample_reads(rng, &g, 500, 40, 300, 0.01);
    for _ in 0..20 {
        let ulen = rng.range(1, 6);
        let unit = rng.dna(ulen);
        let n = rng.range(30, 200);
        let mut t: Vec<u8> = unit.iter().copied().cycle().take(n).collect();
        let flank = rng.dna(20);
        t.extend_from_slice(&flank);
        r.push(t);
    }
    rng.shuffle(&mut r);
    r
}

fn insert_chunked(t: &mut KmerTable, o: &mut Oracle, r: &[u8], rng: &mut Rng) {
    let mut s = 0;
    while s < r.len() {
        let e = (s + rng.range(1, 80)).min(r.len());
        t.insert_chunk(&r[s..e]);
        o.insert_chunk(&r[s..e]);
        s = e;
    }
    t.end_read();
    o.end_read();
}

/// Every k-mer of every read, every oracle k-mer and random k-mers: the
/// table answers as the oracle (k-mers may have been deleted).
fn same_as_oracle(t: &KmerTable, o: &Oracle, rs: &[Vec<u8>], rng: &mut Rng, ctx: &str) {
    let k = t.k();
    for r in rs {
        for x in unihash::codec::canonical_kmers(r, k) {
            assert_eq!(t.contains(x), o.contains(x), "{ctx}: read k-mer");
        }
    }
    for &x in o.iter() {
        assert!(t.contains(x), "{ctx}: oracle k-mer missing");
    }
    for _ in 0..2000 {
        let x = rng.next_u64() & unihash::codec::kmer_mask(k);
        assert_eq!(t.contains(x), o.contains(x), "{ctx}: random k-mer");
    }
}

#[test]
fn shortcuts_match_oracle_and_are_used() {
    for seed in seeds(&[1301, 1302]) {
        for c in configs() {
            let ctx = format!("seed={seed} cfg={c:?}");
            let mut rng = Rng::new(seed);
            let rs = reads(&mut rng);
            let mut t = KmerTable::new(c.clone()).unwrap();
            let mut o = Oracle::new(c.k);
            for (i, r) in rs.iter().enumerate() {
                if i % 3 == 0 {
                    insert_chunked(&mut t, &mut o, r, &mut rng);
                } else {
                    t.insert_read(r);
                    o.insert_read(r);
                }
            }
            check_against_oracle(&t, &o, &rs, &mut rng, 2000, &ctx);
            check_invariants(&t, &o, &ctx);
            t.verify_index().unwrap_or_else(|e| panic!("{ctx}: {e}"));
            if c.run_indexing == RunIndexing::Eager {
                assert_eq!(t.redundant_kmers(), 0, "{ctx}: k-mer stored twice");
            }
            let s = t.insert_stats();
            assert!(s.extensions > 0, "{ctx}: no hit extension");
            assert!(s.cache_reuses > 0, "{ctx}: no candidate reuse");
            assert_eq!(
                s.lookup.queries + s.extensions,
                s.kmers_processed,
                "{ctx}: every k-mer is extended or looked up"
            );
        }
    }
}

/// Deletions and defragmentation between reads change pieces and index
/// entries under the minimizers a read just used; the next read (often a
/// suffix of a previous one, so its first minimizers are the cached ones)
/// must not reuse anything stale.
#[test]
fn edits_between_reads_leave_nothing_stale() {
    for seed in seeds(&[1311, 1312, 1313]) {
        for c in configs() {
            let ctx = format!("seed={seed} cfg={c:?}");
            let k = c.k;
            let mut rng = Rng::new(seed);
            let rs = reads(&mut rng);
            let mut t = KmerTable::new(c.clone()).unwrap();
            let mut o = Oracle::new(k);
            let mut prev: Option<Vec<u8>> = None;
            for (i, r) in rs.iter().take(300).enumerate() {
                insert_chunked(&mut t, &mut o, r, &mut rng);
                if r.len() > k + 10 {
                    // Delete a stretch inside the read just inserted (its
                    // piece is cut, items move to a new id) ...
                    let a = rng.below((r.len() - k) as u64) as usize;
                    let b = (a + k + rng.below(20) as usize).min(r.len());
                    t.delete_read(&r[a..b]);
                    o.delete_read(&r[a..b]);
                    if i % 7 == 0 {
                        t.defrag();
                    } else if i % 5 == 0 {
                        t.defrag_incremental();
                    }
                    // ... then insert its tail and a tail of the previous read.
                    let tail = &r[r.len().saturating_sub(k + 8)..];
                    t.insert_read(tail);
                    o.insert_read(tail);
                    if let Some(p) = &prev {
                        let tail = &p[p.len().saturating_sub(k + 20)..];
                        insert_chunked(&mut t, &mut o, tail, &mut rng);
                    }
                }
                prev = Some(r.clone());
            }
            same_as_oracle(&t, &o, &rs, &mut rng, &ctx);
            t.verify_index().unwrap_or_else(|e| panic!("{ctx}: {e}"));
            // Eager: no k-mer in two different pieces (with deferred run
            // indexing, a k-mer repeated inside a run is stored twice in its
            // piece, and a cut may separate the copies).
            let eager = c.run_indexing == RunIndexing::Eager;
            let mut owner = std::collections::HashMap::new();
            t.for_each_stored_kmer(|pid, _, x| {
                let p = *owner.entry(x).or_insert(pid);
                assert!(!eager || p == pid, "{ctx}: k-mer stored in two pieces");
            });
            assert_eq!(owner.len(), o.len(), "{ctx}");
        }
    }
}

// ----------------------------------------------------------------------
// Concurrent writers
// ----------------------------------------------------------------------

use unihash::codec::{canonicalize, encode_kmer};
use unihash::ConcurrentKmerTable;

fn ckmer(s: &[u8], k: usize) -> u64 {
    canonicalize(encode_kmer(s).unwrap(), k)
}

/// A writer's extension cursor outlives a chunk; meanwhile another thread
/// deletes the next k-mer of the cursor's piece. The writer must see the
/// piece's new view and insert that k-mer again.
#[test]
fn writer_extension_sees_concurrent_deletions() {
    for seed in seeds(&[1321]) {
        for c in [Config::default(), configs()[2].clone()] {
            let k = c.k;
            let mut rng = Rng::new(seed);
            for trial in 0..100 {
                let ctx = format!("seed={seed} trial={trial} cfg={c:?}");
                let t = ConcurrentKmerTable::with_shards(c.clone(), 16).unwrap();
                let r = rng.dna(300);
                t.writer().insert_read(&r);
                let mut w = t.writer();
                w.insert_chunk(&r[50..150]); // all present: cursor in r's piece
                let next = ckmer(&r[151 - k..151], k);
                assert!(t.delete(next), "{ctx}");
                w.insert_chunk(&r[150..250]);
                w.end_read();
                drop(w);
                assert!(t.contains(next), "{ctx}: re-inserted k-mer missing");
                let tt = t.into_table();
                assert_eq!(tt.redundant_kmers(), 0, "{ctx}");
                tt.verify_index().unwrap_or_else(|e| panic!("{ctx}: {e}"));
            }
        }
    }
}

/// A writer's cached candidates outlive a chunk (it releases its shard
/// lock); meanwhile another writer stores the next k-mer under the same
/// minimizer. The first writer must find it, not store it a second time.
#[test]
fn writer_candidates_see_other_writers() {
    for seed in seeds(&[1331]) {
        let c = Config::default();
        let k = c.k;
        let mut rng = Rng::new(seed);
        let t = ConcurrentKmerTable::with_shards(c.clone(), 4).unwrap();
        let mut stats = 0;
        for _ in 0..300 {
            let r = rng.dna(200);
            let mut w = t.writer();
            w.insert_chunk(&r[..100]); // all absent: candidates of the last minimizer cached
            t.writer().insert_read(&r[100 - k + 1..140]); // stores the next k-mers
            w.insert_chunk(&r[100..]);
            w.end_read();
            stats += 1;
        }
        let tt = t.into_table();
        assert_eq!(tt.redundant_kmers(), 0, "seed={seed}: k-mer stored twice");
        tt.verify_index().unwrap_or_else(|e| panic!("{e}"));
        assert!(stats > 0);
    }
}

/// One writer goes through the same shortcuts as `KmerTable` (same pieces,
/// same counters).
#[test]
fn one_writer_takes_the_same_shortcuts() {
    for seed in seeds(&[1341]) {
        for c in configs() {
            if c.run_indexing != RunIndexing::Eager {
                continue;
            }
            let ctx = format!("seed={seed} cfg={c:?}");
            let mut rng = Rng::new(seed);
            let rs = reads(&mut rng);
            let mut st = KmerTable::new(c.clone()).unwrap();
            let ct = ConcurrentKmerTable::with_shards(c.clone(), 1).unwrap();
            {
                let mut w = ct.writer();
                for r in &rs {
                    st.insert_read(r);
                    w.insert_read(r);
                }
            }
            let ct = ct.into_table();
            let (a, b) = (st.insert_stats(), ct.insert_stats());
            assert_eq!(a.present_hits, b.present_hits, "{ctx}");
            assert_eq!(a.extensions, b.extensions, "{ctx}");
            assert!(b.extensions > 0 && b.cache_reuses > 0, "{ctx}");
            assert_eq!(st.stats().stored_kmers, ct.stats().stored_kmers, "{ctx}");
            assert_eq!(ct.redundant_kmers(), 0, "{ctx}");
        }
    }
}
