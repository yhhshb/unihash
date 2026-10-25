//! Churn benchmark: alternating deletions and re-insertions, under several
//! defragmentation policies, to choose when to call `KmerTable::defrag`.
//!
//! ```text
//! unihash-churn [INPUT...] [--synthetic-repeats LEN] [--max-records N] [--unit BASES] [--ops N] [--lag N]
//!               [--policy never|every:N|overhead:X|inc:X|inc:X+full:Y]...
//!               [--check] [--seed S]
//!               [--max-piece-len N] [--k K] [--l L]
//! ```
//!
//! The table is built from the inputs and defragmented once (the starting
//! point of every policy). Then each operation deletes one *unit* (a random
//! window of `--unit` bases of a record, or a whole record when records are
//! shorter, i.e. reads) and re-inserts the unit deleted `--lag` operations
//! earlier, so the table size stays steady while its pieces fragment. Each
//! policy runs the same operation sequence from a fresh table.

use std::time::{Duration, Instant};
use unihash::codec::canonical_kmers;
use unihash::rng::Rng;
use unihash::{Config, KmerTable, Oracle};

#[derive(Clone, Copy, Debug)]
enum Policy {
    Never,
    Every(u64),
    /// Full `defrag` when `overhead(k) > X`.
    Overhead(f64),
    /// `defrag_incremental` when `incremental_overhead() > X`; full
    /// `defrag` when `piece_overhead(k)` exceeds Y (if given).
    Incremental(f64, Option<f64>),
}

impl std::fmt::Display for Policy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Policy::Never => write!(f, "never"),
            Policy::Every(n) => write!(f, "every:{n}"),
            Policy::Overhead(x) => write!(f, "overhead:{x}"),
            Policy::Incremental(x, None) => write!(f, "inc:{x}"),
            Policy::Incremental(x, Some(y)) => write!(f, "inc:{x}+full:{y}"),
        }
    }
}

fn parse_policy(s: &str) -> Result<Policy, String> {
    if s == "never" {
        return Ok(Policy::Never);
    }
    if let Some(n) = s.strip_prefix("every:") {
        return n
            .parse()
            .map(Policy::Every)
            .map_err(|e| format!("{s}: {e}"));
    }
    if let Some(rest) = s.strip_prefix("inc:") {
        let (x, y) = match rest.split_once("+full:") {
            Some((x, y)) => (x, Some(y.parse::<f64>().map_err(|e| format!("{s}: {e}"))?)),
            None => (rest, None),
        };
        let x = x.parse::<f64>().map_err(|e| format!("{s}: {e}"))?;
        return Ok(Policy::Incremental(x, y));
    }
    if let Some(x) = s.strip_prefix("overhead:") {
        return x
            .parse()
            .map(Policy::Overhead)
            .map_err(|e| format!("{s}: {e}"));
    }
    Err(format!(
        "unknown policy {s} (never | every:N | overhead:X | inc:X | inc:X+full:Y)"
    ))
}

struct Args {
    inputs: Vec<String>,
    max_records: usize,
    unit: usize,
    ops: u64,
    lag: usize,
    policies: Vec<Policy>,
    check: bool,
    seed: u64,
    synthetic: usize,
    cfg: Config,
}

fn parse_args() -> Args {
    let mut a = Args {
        inputs: Vec::new(),
        max_records: usize::MAX,
        unit: 1000,
        ops: 20_000,
        lag: 1000,
        policies: Vec::new(),
        check: false,
        seed: 1,
        synthetic: 0,
        cfg: Config::default(),
    };
    let mut it = std::env::args().skip(1);
    let die = |m: String| -> ! {
        eprintln!("unihash-churn: {m}");
        std::process::exit(2)
    };
    while let Some(x) = it.next() {
        let mut val = || {
            it.next()
                .unwrap_or_else(|| die(format!("{x} needs a value")))
        };
        match x.as_str() {
            "--max-records" => {
                a.max_records = val().parse().unwrap_or_else(|e| die(format!("{e}")))
            }
            "--unit" => a.unit = val().parse().unwrap_or_else(|e| die(format!("{e}"))),
            "--ops" => a.ops = val().parse().unwrap_or_else(|e| die(format!("{e}"))),
            "--lag" => a.lag = val().parse().unwrap_or_else(|e| die(format!("{e}"))),
            "--seed" => a.seed = val().parse().unwrap_or_else(|e| die(format!("{e}"))),
            "--policy" => a
                .policies
                .push(parse_policy(&val()).unwrap_or_else(|e| die(e))),
            "--k" => a.cfg.k = val().parse().unwrap_or_else(|e| die(format!("{e}"))),
            "--l" => a.cfg.l = val().parse().unwrap_or_else(|e| die(format!("{e}"))),
            "--max-piece-len" => {
                a.cfg.max_piece_len = val().parse().unwrap_or_else(|e| die(format!("{e}")))
            }
            "--check" => a.check = true,
            "--synthetic-repeats" => {
                a.synthetic = val().parse().unwrap_or_else(|e| die(format!("{e}")))
            }
            "-h" | "--help" => {
                println!(
                    "{}",
                    include_str!("unihash_churn.rs")
                        .lines()
                        .take(15)
                        .collect::<Vec<_>>()
                        .join("\n")
                );
                std::process::exit(0);
            }
            _ if x.starts_with('-') => die(format!("unknown option {x}")),
            _ => a.inputs.push(x),
        }
    }
    if a.inputs.is_empty() && a.synthetic == 0 {
        die("no input".into());
    }
    if a.policies.is_empty() {
        a.policies = vec![
            Policy::Never,
            Policy::Overhead(0.02),
            Policy::Overhead(0.05),
            Policy::Overhead(0.10),
            Policy::Overhead(0.20),
            Policy::Overhead(0.40),
        ];
    }
    a
}

fn load(a: &Args) -> Vec<Vec<u8>> {
    let mut recs = Vec::new();
    if a.synthetic > 0 {
        recs.push(unihash::rng::repeat_rich_genome(
            &mut Rng::new(a.seed ^ 77),
            a.synthetic,
        ));
    }
    for p in &a.inputs {
        let mut r = unihash::fastx::open(std::path::Path::new(p))
            .unwrap_or_else(|e| panic!("cannot open {p}: {e}"));
        let mut seq = Vec::new();
        while recs.len() < a.max_records && r.next_record(&mut seq).expect("read error") {
            recs.push(std::mem::take(&mut seq));
        }
    }
    recs
}

/// The unit deleted by operation `i`: (record, start, end).
fn units(recs: &[Vec<u8>], unit: usize, n: usize, rng: &mut Rng) -> Vec<(usize, usize, usize)> {
    let total: usize = recs.iter().map(|r| r.len()).sum();
    let cum: Vec<usize> = recs
        .iter()
        .scan(0usize, |s, r| {
            *s += r.len();
            Some(*s)
        })
        .collect();
    (0..n)
        .map(|_| {
            let g = rng.below(total as u64) as usize;
            let ri = cum.partition_point(|&c| c <= g);
            let r = &recs[ri];
            if r.len() <= unit {
                (ri, 0, r.len())
            } else {
                let a = rng.range(0, r.len() - unit);
                (ri, a, a + unit)
            }
        })
        .collect()
}

/// Exact set equality with the oracle — every stored k-mer is in the
/// oracle, every oracle k-mer is found by the table — plus the query sample
/// and an exact recomputation of the index. Exits with status 1 on failure.
fn check_sets(t: &KmerTable, o: &Oracle, sample: &[u64], what: &str) {
    let mut extra = 0u64; // stored in the table, not in the oracle
    t.for_each_stored_kmer(|_, _, x| extra += !o.contains(x) as u64);
    let mut missing = 0u64; // in the oracle, not found by the table
    for &x in o.iter() {
        missing += !t.contains(x) as u64;
    }
    let mut wrong = 0u64;
    for &x in sample {
        wrong += (t.contains(x) != o.contains(x)) as u64;
    }
    let index = t.verify_index();
    println!(
        "  check vs oracle ({what}): {} k-mers in both, {extra} only in unihash, {missing} only in the HashSet, {wrong} wrong sample answers, index {}",
        o.len() as u64 - missing,
        match &index {
            Ok(()) => "exact".to_string(),
            Err(e) => e.clone(),
        }
    );
    if extra + missing + wrong > 0 || index.is_err() {
        std::process::exit(1);
    }
}

fn mb(b: f64) -> String {
    format!("{:.2} MB", b / 1e6)
}

fn main() {
    let a = parse_args();
    let recs = load(&a);
    let bases: usize = recs.iter().map(|r| r.len()).sum();
    let mut rng = Rng::new(a.seed);
    let ops = units(&recs, a.unit, a.ops as usize, &mut rng);
    println!("== unihash-churn ==");
    println!(
        "input: {:?}, {} records, {bases} bases",
        a.inputs,
        recs.len()
    );
    println!(
        "operations: {} × (delete one unit of ≤ {} bases, re-insert the unit deleted {} operations earlier)",
        a.ops, a.unit, a.lag
    );

    // Initial build + defrag, timed once.
    let t0 = Instant::now();
    let mut base = KmerTable::new(a.cfg.clone()).unwrap();
    for r in &recs {
        base.insert_read(r);
    }
    let t_build = t0.elapsed();
    base.shrink_to_fit();
    let b0 = base.table_bytes();
    let p0 = base.stats().pieces;
    let t1 = Instant::now();
    let rep = base.defrag();
    let t_defrag0 = t1.elapsed();
    base.shrink_to_fit();
    println!(
        "build {:.2} s; table {} in {p0} pieces; initial defrag {:.3} s: {rep}",
        t_build.as_secs_f64(),
        mb(b0 as f64),
        t_defrag0.as_secs_f64()
    );
    let start_bytes = base.table_bytes() as f64;
    println!(
        "starting table {} ({:.3} B/k-mer)\n",
        mb(start_bytes),
        start_bytes / base.stats().stored_kmers as f64
    );

    // Query sample: k-mers of 3000 units (some deleted, most present).
    let mut qrng = Rng::new(a.seed ^ 0xabcdef);
    let qunits = units(&recs, a.unit, 3000, &mut qrng);
    let mut qkmers: Vec<u64> = qunits
        .iter()
        .flat_map(|&(ri, s, e)| canonical_kmers(&recs[ri][s..e], a.cfg.k).collect::<Vec<_>>())
        .collect();
    qrng.shuffle(&mut qkmers);

    println!(
        "{:<20} {:>9} {:>9} {:>9} {:>10} {:>9} {:>11} {:>11} {:>11} {:>9} {:>9}",
        "policy",
        "full+inc",
        "ops s",
        "defrag s",
        "total s",
        "defrag %",
        "mean table",
        "max table",
        "end table",
        "random",
        "stream"
    );
    for &pol in &a.policies {
        let mut t = base.clone();
        let mut oracle = a.check.then(|| {
            let mut o = Oracle::new(a.cfg.k);
            for r in &recs {
                o.insert_read(r);
            }
            o
        });
        let mut queue = std::collections::VecDeque::new();
        let (mut t_ops, mut t_def) = (Duration::ZERO, Duration::ZERO);
        let mut defrags = 0u64;
        let mut inc_defrags = 0u64;
        let (mut sum_b, mut n_b, mut max_b) = (0f64, 0u64, 0f64);
        let sample_every = (a.ops / 200).max(1);
        for (i, &(ri, s, e)) in ops.iter().enumerate() {
            let u = &recs[ri][s..e];
            let t1 = Instant::now();
            t.delete_read(u);
            queue.push_back((ri, s, e));
            let back = if queue.len() > a.lag {
                queue.pop_front()
            } else {
                None
            };
            if let Some((rj, s2, e2)) = back {
                t.insert_read(&recs[rj][s2..e2]);
            }
            t_ops += t1.elapsed();
            if let Some(o) = oracle.as_mut() {
                o.delete_read(u);
                if let Some((rj, s2, e2)) = back {
                    o.insert_read(&recs[rj][s2..e2]);
                }
            }
            let i = i as u64 + 1;
            // 1 = full defrag, 2 = incremental
            let action = match pol {
                Policy::Never => 0,
                Policy::Every(n) => (i % n == 0) as u8,
                Policy::Overhead(x) => (t.fragmentation().overhead(a.cfg.k) > x) as u8,
                Policy::Incremental(x, y) => {
                    let f = t.fragmentation();
                    if y.is_some_and(|y| f.piece_overhead(a.cfg.k) > y) {
                        1
                    } else if f.incremental_overhead() > x {
                        2
                    } else {
                        0
                    }
                }
            };
            if action > 0 {
                let t2 = Instant::now();
                if action == 1 {
                    t.defrag();
                    defrags += 1;
                } else {
                    t.defrag_incremental();
                    inc_defrags += 1;
                }
                t_def += t2.elapsed();
            }
            if i % sample_every == 0 {
                let b = t.table_bytes() as f64;
                sum_b += b;
                n_b += 1;
                max_b = max_b.max(b);
            }
        }
        let end_b = t.table_bytes() as f64;
        // Query speed in the final (possibly fragmented) state.
        let tq = Instant::now();
        let mut hits = 0u64;
        for &x in &qkmers {
            hits += t.contains(x) as u64;
        }
        let rnd = qkmers.len() as f64 / tq.elapsed().as_secs_f64() / 1e6;
        let ts = Instant::now();
        let mut n_str = 0u64;
        let mut sq = t.query_stream();
        for &(ri, s, e) in &qunits {
            sq.push(&recs[ri][s..e], |_, p| {
                n_str += 1;
                hits += p as u64;
            });
            sq.end();
        }
        let strm = n_str as f64 / ts.elapsed().as_secs_f64() / 1e6;
        let total = t_ops + t_def;
        println!(
            "{:<20} {:>9} {:>9.2} {:>9.2} {:>10.2} {:>8.1}% {:>11} {:>11} {:>11} {:>7.2} M {:>7.2} M",
            pol.to_string(),
            format!("{defrags}+{inc_defrags}"),
            t_ops.as_secs_f64(),
            t_def.as_secs_f64(),
            total.as_secs_f64(),
            100.0 * t_def.as_secs_f64() / total.as_secs_f64(),
            mb(sum_b / n_b.max(1) as f64),
            mb(max_b),
            mb(end_b),
            rnd,
            strm
        );
        if let Some(o) = &oracle {
            check_sets(&t, o, &qkmers, "final state");
        }
        if matches!(pol, Policy::Never) {
            let f = t.fragmentation();
            println!(
                "  never-defragmented end state: {} pieces (start {}), {} dead bases, {} B index slack, {}",
                f.pieces,
                base.stats().pieces,
                f.dead_bases,
                f.index_slack_bytes,
                t.stats().edits
            );
            let t3 = Instant::now();
            let r = t.defrag();
            println!(
                "  one defrag at the end: {:.3} s, {r}",
                t3.elapsed().as_secs_f64()
            );
            let reclaimed = 1.0 - r.table_bytes_after as f64 / r.table_bytes_before as f64;
            println!(
                "  overhead estimate {:.1} % vs reclaimed {:.1} %",
                100.0 * f.overhead(a.cfg.k),
                100.0 * reclaimed
            );
            if let Some(o) = &oracle {
                check_sets(&t, o, &qkmers, "after the final defrag");
            }
        }
        let _ = hits;
    }
}
