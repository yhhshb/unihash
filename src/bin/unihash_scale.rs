//! unihash-scale — thread scaling and input-size scaling of insertion and
//! queries.
//!
//! The input is loaded into memory first (long records are cut into
//! 64 kbp work items overlapping by k − 1 bases, which keeps every k-mer),
//! so only table work is timed. For each thread count the tool builds a
//! [`ConcurrentKmerTable`] with that many writers (work items handed out
//! dynamically), queries it with as many threads, converts it with
//! `into_table`, and queries the single-threaded [`KmerTable`] lock-free with
//! as many threads. The single-threaded `KmerTable` build is the reference
//! for speed-ups and, with `--check`, for exact set equality. With
//! `--sizes`, the same measures are taken on growing prefixes of the input.

use std::process::exit;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{Duration, Instant};
use unihash::codec::{canonical_kmers, kmer_mask};
use unihash::rng::{repeat_rich_genome, sample_reads, Rng};
use unihash::{fastx, ConcurrentKmerTable, Config, KmerTable};

#[global_allocator]
static GLOBAL: counting::Counting = counting::Counting;

#[path = "common/counting.rs"]
mod counting;

const USAGE: &str = "\
unihash-scale — thread and input-size scaling of insertion and queries

USAGE:
    unihash-scale [OPTIONS] [INPUT.fa|.fq[.gz] ...]

TABLE:
    --k N, --l N, --bucket-bits N, --max-piece-len N, --hot-threshold N
                          as in unihash-bench
    --shards N            layer-1 index shards of the concurrent table (default 1024)
    --capacity            size the concurrent table for the distinct k-mers of the
                          input (known from the single-threaded build) with
                          ConcurrentKmerTable::with_capacity: shards scale with
                          the table and never grow (overrides --shards)

SCALING:
    --threads LIST        thread counts, e.g. 1,2,4,8 (default: powers of two up
                          to the available cores, and the core count)
    --sizes LIST          also measure growing prefixes of the input, as
                          fractions, e.g. 0.125,0.25,0.5,1
    --size-threads N      threads used for --sizes (default: largest --threads)
    --reps N              repetitions of each build; the median time is shown (default 1)
    --check               check each multi-threaded table against the
                          single-threaded one (same k-mers, none stored twice)

INPUT:
    --max-records N       use the first N records of each input file
    --synthetic LEN       uniform random genome of LEN bases
    --synthetic-repeats LEN
                          repeat-rich synthetic genome of LEN bases
    --reads N             sample N reads (100-150 b, 0.5 % errors) from the
                          synthetic genome instead of inserting it whole
    --seed N              RNG seed (default 42)

QUERIES:
    --query-kmers N       positive k-mers sampled from the input (default 10000000)

DELETION:
    --delete-frac F       also time concurrent deletion: build the table, then
                          delete every (1/F)-th work item (whole reads, or
                          64 kbp windows) with as many threads, against
                          KmerTable::delete_read on one thread (default 0: off)
    --neg N               random k-mers (default 10000000)
";

/// Work item length for long records.
const ITEM: usize = 1 << 16;
/// Items handed to a writer at a time (fewer when there are few items).
const BATCH: usize = 32;

struct Args {
    cfg: Config,
    shards: usize,
    capacity: bool,
    threads: Vec<usize>,
    sizes: Vec<f64>,
    size_threads: Option<usize>,
    reps: usize,
    check: bool,
    inputs: Vec<String>,
    max_records: Option<usize>,
    synthetic: Option<(usize, bool)>,
    reads: Option<usize>,
    seed: u64,
    query_kmers: usize,
    neg: usize,
    delete_frac: f64,
}

fn die(msg: impl std::fmt::Display) -> ! {
    eprintln!("error: {msg}\n\nRun with --help for usage.");
    exit(2)
}

fn parse_args() -> Args {
    let cores = std::thread::available_parallelism().map_or(1, |n| n.get());
    let mut threads: Vec<usize> = std::iter::successors(Some(1usize), |&t| Some(t * 2))
        .take_while(|&t| t <= cores)
        .collect();
    if *threads.last().unwrap() != cores {
        threads.push(cores);
    }
    let mut a = Args {
        cfg: Config::default(),
        shards: unihash::concurrent::DEFAULT_SHARDS,
        capacity: false,
        threads,
        sizes: Vec::new(),
        size_threads: None,
        reps: 1,
        check: false,
        inputs: Vec::new(),
        max_records: None,
        synthetic: None,
        reads: None,
        seed: 42,
        query_kmers: 10_000_000,
        neg: 10_000_000,
        delete_frac: 0.0,
    };
    let mut it = std::env::args().skip(1);
    fn val(it: &mut impl Iterator<Item = String>, f: &str) -> String {
        it.next()
            .unwrap_or_else(|| die(format!("{f} needs a value")))
    }
    fn num<T: std::str::FromStr>(it: &mut impl Iterator<Item = String>, f: &str) -> T {
        let v = val(it, f);
        v.replace('_', "")
            .parse()
            .unwrap_or_else(|_| die(format!("invalid value '{v}' for {f}")))
    }
    fn list<T: std::str::FromStr>(it: &mut impl Iterator<Item = String>, f: &str) -> Vec<T> {
        let v = val(it, f);
        v.split(',')
            .map(|x| {
                x.trim()
                    .parse()
                    .unwrap_or_else(|_| die(format!("invalid value '{x}' in {f}")))
            })
            .collect()
    }
    while let Some(x) = it.next() {
        match x.as_str() {
            "-h" | "--help" => {
                print!("{USAGE}");
                exit(0)
            }
            "--k" => a.cfg.k = num(&mut it, "--k"),
            "--l" => a.cfg.l = num(&mut it, "--l"),
            "--bucket-bits" => a.cfg.bucket_bits = num(&mut it, "--bucket-bits"),
            "--max-piece-len" => a.cfg.max_piece_len = num(&mut it, "--max-piece-len"),
            "--hot-threshold" => a.cfg.hot_threshold = num(&mut it, "--hot-threshold"),
            "--shards" => a.shards = num(&mut it, "--shards"),
            "--capacity" => a.capacity = true,
            "--threads" => a.threads = list(&mut it, "--threads"),
            "--sizes" => a.sizes = list(&mut it, "--sizes"),
            "--size-threads" => a.size_threads = Some(num(&mut it, "--size-threads")),
            "--reps" => a.reps = num::<usize>(&mut it, "--reps").max(1),
            "--check" => a.check = true,
            "--max-records" => a.max_records = Some(num(&mut it, "--max-records")),
            "--synthetic" => a.synthetic = Some((num(&mut it, "--synthetic"), false)),
            "--synthetic-repeats" => {
                a.synthetic = Some((num(&mut it, "--synthetic-repeats"), true))
            }
            "--reads" => a.reads = Some(num(&mut it, "--reads")),
            "--seed" => a.seed = num(&mut it, "--seed"),
            "--query-kmers" => a.query_kmers = num(&mut it, "--query-kmers"),
            "--neg" => a.neg = num(&mut it, "--neg"),
            "--delete-frac" => a.delete_frac = num(&mut it, "--delete-frac"),
            _ if x.starts_with('-') => die(format!("unknown option {x}")),
            _ => a.inputs.push(x),
        }
    }
    if a.inputs.is_empty() && a.synthetic.is_none() {
        die("no input: give FASTA/FASTQ files or --synthetic/--synthetic-repeats");
    }
    if let Err(e) = a.cfg.validate() {
        die(e);
    }
    if a.threads.is_empty() || a.threads.contains(&0) {
        die("--threads needs positive counts");
    }
    if a.sizes.iter().any(|&f| !(f > 0.0 && f <= 1.0)) {
        die("--sizes fractions must be in (0, 1]");
    }
    a
}

/// Load the input as work items.
fn load(a: &Args) -> Vec<Vec<u8>> {
    let k = a.cfg.k;
    let mut items = Vec::new();
    let push = |rec: &[u8], items: &mut Vec<Vec<u8>>| {
        if rec.len() <= ITEM {
            items.push(rec.to_vec());
            return;
        }
        let mut s = 0;
        while s + k <= rec.len() {
            let e = (s + ITEM).min(rec.len());
            items.push(rec[s..e].to_vec());
            if e == rec.len() {
                break;
            }
            s = e + 1 - k;
        }
    };
    if let Some((len, rep)) = a.synthetic {
        let mut rng = Rng::new(a.seed);
        let g = if rep {
            repeat_rich_genome(&mut rng, len)
        } else {
            rng.dna(len)
        };
        match a.reads {
            Some(n) => items.extend(sample_reads(&mut rng, &g, n, 100, 150, 0.005)),
            None => push(&g, &mut items),
        }
    }
    for path in &a.inputs {
        let mut r = fastx::open(std::path::Path::new(path)).unwrap_or_else(|e| die(e));
        let mut seq = Vec::new();
        let mut n = 0;
        while r.next_record(&mut seq).unwrap_or_else(|e| die(e)) {
            push(&seq, &mut items);
            n += 1;
            if a.max_records.is_some_and(|m| n >= m) {
                break;
            }
        }
    }
    items
}

fn secs(d: Duration) -> f64 {
    d.as_secs_f64()
}

fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[v.len() / 2]
}

fn single_build(cfg: &Config, items: &[Vec<u8>]) -> (KmerTable, Duration, usize) {
    let base = counting::current();
    counting::reset_peak();
    let t0 = Instant::now();
    let mut t = KmerTable::new(cfg.clone()).unwrap();
    for r in items {
        t.insert_read(r);
    }
    let d = t0.elapsed();
    let peak = counting::peak() - base;
    t.shrink_to_fit();
    (t, d, peak)
}

struct Built {
    table: ConcurrentKmerTable,
    time: Duration,
    peak: usize,
    heap: usize,
    /// Writer context switches: voluntary (slept, e.g. on a lock) and
    /// involuntary (preempted); minor page faults of the process.
    csw: (u64, u64),
    faults: u64,
    busy: String,
}

/// (voluntary, involuntary) context switches of the calling thread.
fn thread_csw() -> (u64, u64) {
    let s = std::fs::read_to_string("/proc/thread-self/status").unwrap_or_default();
    let get = |key: &str| {
        s.lines()
            .find(|l| l.starts_with(key))
            .and_then(|l| l.split_whitespace().nth(1))
            .and_then(|v| v.parse().ok())
            .unwrap_or(0)
    };
    (
        get("voluntary_ctxt_switches:"),
        get("nonvoluntary_ctxt_switches:"),
    )
}

/// Minor page faults of the process (live and finished threads).
fn minor_faults() -> u64 {
    let s = std::fs::read_to_string("/proc/self/stat").unwrap_or_default();
    // Fields after the parenthesised command name; minflt is field 10.
    s.rsplit_once(')')
        .and_then(|(_, rest)| rest.split_whitespace().nth(7))
        .and_then(|v| v.parse().ok())
        .unwrap_or(0)
}

/// CPU seconds used by the process so far (user + system; 10 ms ticks).
fn cpu_seconds() -> f64 {
    let s = std::fs::read_to_string("/proc/self/stat").unwrap_or_default();
    let f: Vec<u64> = s
        .rsplit_once(')')
        .map(|(_, rest)| {
            rest.split_whitespace()
                .skip(11)
                .take(2)
                .filter_map(|v| v.parse().ok())
                .collect()
        })
        .unwrap_or_default();
    f.iter().sum::<u64>() as f64 / 100.0
}

/// Share of `threads × wall` the process spent on a CPU: below 100 % the
/// threads waited (locks, sleeps) or were preempted by other processes.
fn busy(cpu0: f64, wall: Duration, threads: usize) -> String {
    format!(
        "CPU busy {:.0} % of {threads} thread(s)",
        100.0 * (cpu_seconds() - cpu0) / (wall.as_secs_f64() * threads as f64)
    )
}

fn concurrent_build(
    cfg: &Config,
    shards: usize,
    capacity: Option<u64>,
    items: &[Vec<u8>],
    threads: usize,
) -> Built {
    let base = counting::current();
    counting::reset_peak();
    let faults0 = minor_faults();
    let cpu0 = cpu_seconds();
    let t0 = Instant::now();
    let table = match capacity {
        Some(n) => ConcurrentKmerTable::with_capacity(cfg.clone(), n).unwrap(),
        None => ConcurrentKmerTable::with_shards(cfg.clone(), shards).unwrap(),
    };
    let csw = std::sync::Mutex::new((0u64, 0u64));
    let next = AtomicUsize::new(0);
    let batch = (items.len() / (threads * 16)).clamp(1, BATCH);
    std::thread::scope(|s| {
        for _ in 0..threads {
            let (table, next, csw) = (&table, &next, &csw);
            s.spawn(move || {
                let mut w = table.writer();
                loop {
                    let i = next.fetch_add(batch, Ordering::Relaxed);
                    if i >= items.len() {
                        break;
                    }
                    for r in &items[i..(i + batch).min(items.len())] {
                        w.insert_read(r);
                    }
                }
                drop(w);
                let (v, n) = thread_csw();
                let mut c = csw.lock().unwrap();
                c.0 += v;
                c.1 += n;
            });
        }
    });
    let time = t0.elapsed();
    Built {
        table,
        time,
        peak: counting::peak() - base,
        heap: counting::current() - base,
        csw: csw.into_inner().unwrap(),
        faults: minor_faults() - faults0,
        busy: busy(cpu0, time, threads),
    }
}

/// Run `f` over `n` parts in `threads` threads; returns (hits, time).
fn parallel_count<F: Fn(usize, usize) -> u64 + Sync>(
    n: usize,
    threads: usize,
    f: F,
) -> (u64, Duration) {
    let t0 = Instant::now();
    let hits = std::thread::scope(|s| {
        let per = n.div_ceil(threads);
        let hs: Vec<_> = (0..threads)
            .map(|i| {
                let f = &f;
                s.spawn(move || f((i * per).min(n), ((i + 1) * per).min(n)))
            })
            .collect();
        hs.into_iter().map(|h| h.join().unwrap()).sum()
    });
    (hits, t0.elapsed())
}

struct Queries {
    pos: Vec<u64>,
    neg: Vec<u64>,
    /// Sequences for the streaming query (sampled items).
    seqs: Vec<Vec<u8>>,
    seq_kmers: usize,
}

fn make_queries(a: &Args, items: &[Vec<u8>]) -> Queries {
    let k = a.cfg.k;
    let total: usize = items.iter().map(|r| r.len().saturating_sub(k - 1)).sum();
    let stride = (total / a.query_kmers.max(1)).max(1);
    let mut pos = Vec::with_capacity(a.query_kmers.min(total));
    let mut i = 0usize;
    for r in items {
        for x in canonical_kmers(r, k) {
            if i % stride == 0 {
                pos.push(x);
            }
            i += 1;
        }
    }
    let mut rng = Rng::new(a.seed ^ 0x5eed);
    rng.shuffle(&mut pos);
    let neg: Vec<u64> = (0..a.neg).map(|_| rng.next_u64() & kmer_mask(k)).collect();
    // Streaming: every s-th item, about query_kmers k-mers in all.
    let per_item = total / items.len().max(1);
    let want = (a.query_kmers / per_item.max(1)).max(1);
    let step = (items.len() / want).max(1);
    let seqs: Vec<Vec<u8>> = items.iter().step_by(step).cloned().collect();
    let seq_kmers = seqs.iter().map(|r| r.len().saturating_sub(k - 1)).sum();
    Queries {
        pos,
        neg,
        seqs,
        seq_kmers,
    }
}

fn rate(n: usize, d: Duration) -> f64 {
    n as f64 / secs(d) / 1e6
}

fn check_equal(st: &KmerTable, mt: &KmerTable, what: &str) -> bool {
    // `st` stores each of its k-mers once (eager run indexing). If `mt`
    // contains all of them and stores as many k-mers, it stores exactly
    // these, each once.
    let mut missing = 0u64;
    st.for_each_stored_kmer(|_, _, x| missing += !mt.contains(x) as u64);
    let (a, b) = (st.stats().stored_kmers, mt.stats().stored_kmers);
    let index = mt.verify_index();
    let ok = missing == 0 && a == b && index.is_ok();
    println!(
        "  check {what}: {b} stored k-mers vs {a} single-threaded, {missing} missing, index {} → {}",
        match &index {
            Ok(()) => "exact".to_string(),
            Err(e) => format!("WRONG ({e})"),
        },
        if ok { "OK" } else { "MISMATCH" }
    );
    ok
}

fn main() {
    let a = parse_args();
    let k = a.cfg.k;
    println!("== unihash-scale ==");
    println!(
        "config: {:?}, shards {}, cores {}",
        a.cfg,
        a.shards,
        std::thread::available_parallelism().map_or(1, |n| n.get())
    );
    println!("{}", cpu_note());
    let t0 = Instant::now();
    let items = load(&a);
    let bases: usize = items.iter().map(Vec::len).sum();
    println!(
        "input: {} items, {bases} bases (loaded in {:.1} s)",
        items.len(),
        secs(t0.elapsed())
    );

    // Single-threaded reference.
    let mut st_times = Vec::new();
    let mut st = None;
    for _ in 0..a.reps {
        let (t, d, peak) = single_build(&a.cfg, &items);
        st_times.push(secs(d));
        st = Some((t, peak));
    }
    let (st, st_peak) = st.unwrap();
    let st_time = median(st_times);
    // k-mers looked up during insertion (windows without non-ACGT bases).
    let kmers = st.insert_stats().kmers_processed as usize;
    let distinct = st.stats().stored_kmers as f64;
    let st_bytes = {
        let s = st.stats();
        s.piece_bytes_reserved
            + s.index_bytes_reserved
            + s.index2_bytes_reserved
            + s.layer3_bytes_reserved
            + s.dirty_bytes_reserved
    };
    println!(
        "single-threaded KmerTable: {kmers} k-mers in {:.3} s, {:.2} M k-mers/s, {} distinct k-mers, {:.3} B/k-mer (peak {:.3})",
        st_time,
        kmers as f64 / st_time / 1e6,
        distinct,
        st_bytes as f64 / distinct,
        st_peak as f64 / distinct
    );
    let q = make_queries(&a, &items);
    println!(
        "queries: {} positive, {} random, streaming over {} items ({} k-mers)",
        q.pos.len(),
        q.neg.len(),
        q.seqs.len(),
        q.seq_kmers
    );

    let sizing = if a.capacity {
        format!(
            "with_capacity({distinct}): {} shards",
            ConcurrentKmerTable::with_capacity(a.cfg.clone(), distinct as u64)
                .unwrap()
                .stats()
                .shards
        )
    } else {
        format!("{} shards", a.shards)
    };
    println!("\n-- thread scaling: insertion (ConcurrentKmerTable, {sizing}) --");
    println!("threads | insert s | M k-mers/s | speed-up vs 1 thread | vs KmerTable | efficiency | B/k-mer (heap) | peak B/k-mer | into_table s | B/k-mer after | chunk splits");
    let mut base_rate = None;
    let mut all_ok = true;
    let mut locked = Vec::new();
    for &th in &a.threads {
        let mut times = Vec::new();
        let mut last = None;
        for _ in 0..a.reps {
            let b = concurrent_build(
                &a.cfg,
                a.shards,
                a.capacity.then_some(distinct as u64),
                &items,
                th,
            );
            times.push(secs(b.time));
            last = Some(b);
        }
        let b = last.unwrap();
        let time = median(times);
        let r = kmers as f64 / time / 1e6;
        let r1 = *base_rate.get_or_insert(r);
        let cs = b.table.stats();
        // Locked queries on the concurrent table.
        let (hp, dp) = parallel_count(q.pos.len(), th, |s, e| {
            q.pos[s..e].iter().filter(|&&x| b.table.contains(x)).count() as u64
        });
        let (_, dn) = parallel_count(q.neg.len(), th, |s, e| {
            q.neg[s..e].iter().filter(|&&x| b.table.contains(x)).count() as u64
        });
        assert_eq!(hp as usize, q.pos.len(), "positive query missed");
        let (hs, ds) = parallel_count(q.seqs.len(), th, |s, e| {
            let mut sq = b.table.query_stream();
            let mut h = 0u64;
            for r in &q.seqs[s..e] {
                sq.push(r, |_, p| h += p as u64);
                sq.end();
            }
            h
        });
        assert_eq!(
            hs as usize,
            q.seq_kmers - count_invalid(&q.seqs, a.cfg.k),
            "concurrent streaming missed"
        );
        locked.push((
            th,
            rate(q.pos.len(), dp),
            rate(q.neg.len(), dn),
            rate(q.seq_kmers, ds),
        ));
        let t1 = Instant::now();
        let mt = b.table.into_table();
        let conv = t1.elapsed();
        let ms = mt.stats();
        let mt_bytes = ms.piece_bytes_reserved
            + ms.index_bytes_reserved
            + ms.index2_bytes_reserved
            + ms.layer3_bytes_reserved
            + ms.dirty_bytes_reserved;
        println!(
            "{th:7} | {time:8.3} | {r:10.2} | {:20.2} | {:12.2} | {:9.0} % | {:14.3} | {:12.3} | {:12.3} | {:13.3} | {}",
            r / r1,
            st_time / time,
            100.0 * r / r1 / th as f64,
            b.heap as f64 / distinct,
            b.peak as f64 / distinct,
            secs(conv),
            mt_bytes as f64 / distinct,
            cs.insert.chunk_splits
        );
        let (vol, invol) = b.csw;
        println!(
            "  {}; writer context switches {vol} voluntary (slept), {invol} involuntary (preempted); {} minor page faults; {}; {}",
            cs.contention, b.faults, b.busy, cpu_note()
        );
        if a.check {
            all_ok &= check_equal(&st, &mt, &format!("{th} threads"));
        }
    }

    if a.delete_frac > 0.0 {
        all_ok &= deletion_scaling(&a, &items, &st, distinct);
    }

    println!("\n-- thread scaling: queries (M queries/s; KmerTable is lock-free, ConcurrentKmerTable takes shard read locks) --");
    println!("threads | KmerTable positive | random | streaming | Concurrent positive | random | streaming");
    let mut sq_hits = 0;
    for (i, &th) in a.threads.iter().enumerate() {
        let (hp, dp) = parallel_count(q.pos.len(), th, |s, e| {
            q.pos[s..e].iter().filter(|&&x| st.contains(x)).count() as u64
        });
        let (_, dn) = parallel_count(q.neg.len(), th, |s, e| {
            q.neg[s..e].iter().filter(|&&x| st.contains(x)).count() as u64
        });
        let (hs, ds) = parallel_count(q.seqs.len(), th, |s, e| {
            let mut sq = st.query_stream();
            let mut h = 0u64;
            for r in &q.seqs[s..e] {
                sq.push(r, |_, p| h += p as u64);
                sq.end();
            }
            h
        });
        assert_eq!(hp as usize, q.pos.len(), "positive query missed");
        sq_hits = hs;
        let (_, lp, ln, ls) = locked[i];
        println!(
            "{th:7} | {:18.2} | {:6.2} | {:9.2} | {:19.2} | {:6.2} | {:9.2}",
            rate(q.pos.len(), dp),
            rate(q.neg.len(), dn),
            rate(q.seq_kmers, ds),
            lp,
            ln,
            ls
        );
    }
    assert_eq!(
        sq_hits as usize,
        q.seq_kmers - count_invalid(&q.seqs, k),
        "streaming missed"
    );

    if !a.sizes.is_empty() {
        let th = a
            .size_threads
            .unwrap_or_else(|| *a.threads.iter().max().unwrap());
        println!("\n-- input-size scaling (prefixes of the input; {th} threads for the concurrent build) --");
        println!("fraction | distinct k-mers | KmerTable M k-mers/s | B/k-mer | peak B/k-mer | {th}-thread M k-mers/s | speed-up | queries positive / random (M/s, 1 thread, KmerTable)");
        for &f in &a.sizes {
            let n = ((items.len() as f64 * f).ceil() as usize).clamp(1, items.len());
            let part = &items[..n];
            let (t, d, peak) = single_build(&a.cfg, part);
            let km = t.insert_stats().kmers_processed as usize;
            let s = t.stats();
            let dist = s.stored_kmers as f64;
            let bytes = s.piece_bytes_reserved
                + s.index_bytes_reserved
                + s.index2_bytes_reserved
                + s.layer3_bytes_reserved
                + s.dirty_bytes_reserved;
            let b = concurrent_build(
                &a.cfg,
                a.shards,
                a.capacity.then_some(dist as u64),
                part,
                th,
            );
            drop(b.table);
            let pq = make_queries(&a, part);
            let t1 = Instant::now();
            let hp = pq.pos.iter().filter(|&&x| t.contains(x)).count();
            let dp = t1.elapsed();
            assert_eq!(hp, pq.pos.len());
            let t2 = Instant::now();
            let _ = pq.neg.iter().filter(|&&x| t.contains(x)).count();
            let dn = t2.elapsed();
            println!(
                "{f:8} | {:15} | {:20.2} | {:7.3} | {:12.3} | {:20.2} | {:8.2} | {:.2} / {:.2}",
                s.stored_kmers,
                km as f64 / secs(d) / 1e6,
                bytes as f64 / dist,
                peak as f64 / dist,
                km as f64 / secs(b.time) / 1e6,
                secs(d) / secs(b.time),
                rate(pq.pos.len(), dp),
                rate(pq.neg.len(), dn)
            );
        }
    }
    if a.check {
        println!(
            "\ncheck: {}",
            if all_ok {
                "all multi-threaded tables equal the single-threaded one"
            } else {
                "MISMATCH"
            }
        );
        if !all_ok {
            exit(1);
        }
    }
}

/// Concurrent deletion of every (1/F)-th item, for each thread count, against
/// single-threaded `KmerTable::delete_read`. Returns whether every check
/// passed.
fn deletion_scaling(a: &Args, items: &[Vec<u8>], st: &KmerTable, distinct: f64) -> bool {
    let step = (1.0 / a.delete_frac).round().max(1.0) as usize;
    let dels: Vec<&Vec<u8>> = items.iter().step_by(step).collect();
    let dkmers: usize = dels
        .iter()
        .map(|r| canonical_kmers(r, a.cfg.k).count())
        .sum();
    let mut reference = st.clone();
    let t0 = Instant::now();
    let mut removed = 0u64;
    for r in &dels {
        removed += reference.delete_read(r);
    }
    let st_time = secs(t0.elapsed());
    println!(
        "\n-- thread scaling: deletion (every {step}th item: {} items, {dkmers} k-mers looked up, {removed} removed) --",
        dels.len()
    );
    println!(
        "KmerTable::delete_read, 1 thread: {st_time:.3} s, {:.2} M k-mers/s",
        dkmers as f64 / st_time / 1e6
    );
    println!("threads | delete s | M k-mers/s | speed-up vs 1 thread | vs KmerTable | removed | dead bases after | B/k-mer after into_table + defrag_incremental");
    let mut ok = true;
    let mut r1 = None;
    for &th in &a.threads {
        let b = concurrent_build(
            &a.cfg,
            a.shards,
            a.capacity.then_some(distinct as u64),
            items,
            th,
        );
        let before = b.table.stats();
        let next = AtomicUsize::new(0);
        let batch = (dels.len() / (th * 16)).clamp(1, BATCH);
        let cpu0 = cpu_seconds();
        let t0 = Instant::now();
        let n: u64 = std::thread::scope(|s| {
            let hs: Vec<_> = (0..th)
                .map(|_| {
                    let (t, next, dels) = (&b.table, &next, &dels);
                    s.spawn(move || {
                        let mut n = 0;
                        loop {
                            let i = next.fetch_add(batch, Ordering::Relaxed);
                            if i >= dels.len() {
                                break n;
                            }
                            for r in &dels[i..(i + batch).min(dels.len())] {
                                n += t.delete_read(r);
                            }
                        }
                    })
                })
                .collect();
            hs.into_iter().map(|h| h.join().unwrap()).sum()
        });
        let wall = t0.elapsed();
        let busy_del = busy(cpu0, wall, th);
        let time = secs(wall);
        let after = b.table.stats();
        let rate = dkmers as f64 / time / 1e6;
        let base = *r1.get_or_insert(rate);
        let c = after.contention.since(&before.contention);
        let dead = after.edits.dead_bases;
        let mut mt = b.table.into_table();
        mt.defrag_incremental();
        println!(
            "{th:7} | {time:8.3} | {rate:10.2} | {:20.2} | {:12.2} | {n:7} | {dead:16} | {:.3}",
            rate / base,
            st_time / time,
            mt.table_bytes() as f64 / mt.stats().stored_kmers.max(1) as f64
        );
        println!("  {c}; {busy_del}; {}", cpu_note());
        if a.check {
            ok &= check_equal(&reference, &mt, &format!("{th} threads after deletion"));
        }
    }
    ok
}

/// Speed of one thread on a fixed integer loop (M iterations/s, ~20 ms).
fn cpu_speed() -> f64 {
    let n = 4_000_000u64;
    let t0 = Instant::now();
    let mut x = 0x9e37_79b9_7f4a_7c15u64;
    for i in 0..n {
        x ^= x >> 31;
        x = x.wrapping_mul(0xbf58_476d_1ce4_e5b9).wrapping_add(i);
    }
    std::hint::black_box(x);
    n as f64 / t0.elapsed().as_secs_f64() / 1e6
}

static CPU_BASE: std::sync::OnceLock<f64> = std::sync::OnceLock::new();

/// One-thread CPU speed now, relative to the start of the run: a value well
/// below 1 means the machine was throttled or busy when the row above was
/// measured.
fn cpu_note() -> String {
    let now = cpu_speed();
    let base = *CPU_BASE.get_or_init(|| now);
    format!("CPU check {:.2}× of start", now / base)
}

/// k-mers of `seqs` that contain a non-ACGT character (not answered).
fn count_invalid(seqs: &[Vec<u8>], k: usize) -> usize {
    let all: usize = seqs.iter().map(|r| r.len().saturating_sub(k - 1)).sum();
    let valid: usize = seqs.iter().map(|r| canonical_kmers(r, k).count()).sum();
    all - valid
}
