//! Benchmark and memory-accounting driver.
//!
//! 1. Inserts the input (FASTA/FASTQ files, optionally gzipped, or a
//!    synthetic genome) into the table **alone** and reports its exact heap
//!    (accounted from capacities and cross-checked with a counting
//!    allocator), its peak heap during insertion, and RSS — before the
//!    `HashSet<u64>` oracle is built in a second pass.
//! 2. Builds query workloads that cover each query source **uniformly**
//!    (windows of consecutive bases spread over the whole input, or whole
//!    records sampled at a fixed stride for reads).
//! 3. Times every access pattern on **the same k-mers, measured the same
//!    way**, for the table and the oracle:
//!    pre-parsed k-mers in sequence order and shuffled, and the sequence
//!    itself (parsing only, independent lookups, streaming query).
//! 4. Random negatives, minimizer skew, peak RSS.
//!
//! Run `unihash-bench --help` for options.

use std::collections::HashMap;
use std::path::PathBuf;
use std::process::exit;
use std::time::{Duration, Instant};

use unihash::codec::{decode_kmer, kmer_mask, KmerStream};
use unihash::index::{BLOCK_BYTES, BLOCK_SLOTS};
use unihash::minimizer::lmer_of_order_value;
use unihash::rng::{repeat_rich_genome, sample_reads, Rng};
use unihash::{
    fastx, CandidateScan, Config, KmerTable, Oracle, QueryStats, QueryTrace, RunIndexing,
};

// ---------------------------------------------------------------------------
// Counting allocator: exact live and peak heap bytes of this process.
// ---------------------------------------------------------------------------

#[path = "common/counting.rs"]
mod counting;

#[global_allocator]
static GLOBAL: counting::Counting = counting::Counting;

const USAGE: &str = "\
unihash-bench — benchmark the unihash k-mer presence table

USAGE:
    unihash-bench [OPTIONS] [INPUT.fa|.fq[.gz] ...]

TABLE:
    --k N                 k-mer length (default 31)
    --l N                 minimizer length (default 17)
    --bucket-bits N       log2 of the initial bucket count; the index grows
                          automatically (default 10)
    --max-piece-len N     maximum piece length in bases (default 4096)
    --scan ties|range     candidate check per index entry (default ties)
    --run-indexing deferred|eager
                          index each new k-mer as it is appended (eager,
                          default) or a run of new k-mers when it ends
    --hot-threshold N     index entries that make a minimizer hot and route
                          its k-mers to the fallback layers (default 64; 0 = off)
    --l2 N                fallback minimizer length (default 0 = min(k, ℓ+8))
    --spill-dir DIR       spill the old index to temporary files in DIR while
                          the index grows (default: compact in-memory list)
    --defrag              defragment after insertion (merge pieces that
                          overlap by k-1) before measuring memory and queries

INPUT (files and/or synthetic data):
    --max-records N       use the first N records of each input file
    --synthetic LEN       uniform random genome of LEN bases
    --synthetic-repeats LEN
                          repeat-rich synthetic genome of LEN bases
    --reads N             sample N reads from the synthetic genome instead
                          of inserting it whole
    --read-len A-B        read length range (default 100-150)
    --error-rate E        per-base substitution rate for sampled reads (default 0.005)
    --seed N              RNG seed (default 42)

QUERIES (each source is sampled uniformly; all patterns use the same k-mers):
    --query-kmers N       target k-mers per query workload (default 20000000)
    --window N            bases per window for long records (default 1000);
                          records shorter than this are sampled whole
    --holdout             also query the next --max-records records of each
                          input file (held out from insertion)
    --query FILE          additional query set (e.g. a low-hit read set). Repeatable.
    --neg N               random negative queries (default 1000000)
    --no-oracle           skip the hash-set baseline
    --top N               print the N hottest minimizers (default 10)

COMPARISON:
    --dump DIR            write the insertion inputs, every query workload
                          (window sequences, in-order and shuffled k-mers) and
                          the negatives to DIR, with expected hit counts, for
                          tools/kache/kache_bench (see tools/kache_compare.sh)
";

struct Args {
    cfg: Config,
    inputs: Vec<PathBuf>,
    queries: Vec<PathBuf>,
    holdout: bool,
    max_records: u64,
    synthetic: Option<(usize, bool)>,
    reads: Option<usize>,
    read_len: (usize, usize),
    error_rate: f64,
    seed: u64,
    query_kmers: u64,
    window: usize,
    neg: u64,
    oracle: bool,
    top: usize,
    dump: Option<PathBuf>,
    defrag: bool,
}

fn die(msg: impl std::fmt::Display) -> ! {
    eprintln!("error: {msg}\n\nRun with --help for usage.");
    exit(2)
}

fn parse_args() -> Args {
    parse_args_from(std::env::args().skip(1)).unwrap_or_else(|e| die(e))
}

/// Next value of `flag`, verbatim (paths and names keep their underscores).
fn string(it: &mut impl Iterator<Item = String>, flag: &str) -> Result<String, String> {
    it.next().ok_or_else(|| format!("{flag} needs a value"))
}

/// Next value of `flag` parsed as a number; `_` digit separators are
/// allowed (`--neg 1_000_000`).
fn num<T: std::str::FromStr>(
    it: &mut impl Iterator<Item = String>,
    flag: &str,
) -> Result<T, String> {
    let v = string(it, flag)?;
    v.replace('_', "")
        .parse()
        .map_err(|_| format!("invalid value '{v}' for {flag}"))
}

fn parse_args_from(mut it: impl Iterator<Item = String>) -> Result<Args, String> {
    let mut a = Args {
        cfg: Config {
            bucket_bits: 10,
            ..Config::default()
        },
        inputs: vec![],
        queries: vec![],
        holdout: false,
        max_records: u64::MAX,
        synthetic: None,
        reads: None,
        read_len: (100, 150),
        error_rate: 0.005,
        seed: 42,
        query_kmers: 20_000_000,
        window: 1000,
        neg: 1_000_000,
        oracle: true,
        top: 10,
        dump: None,
        defrag: false,
    };
    while let Some(arg) = it.next() {
        match arg.as_str() {
            "-h" | "--help" => {
                print!("{USAGE}");
                exit(0)
            }
            "--k" => a.cfg.k = num(&mut it, "--k")?,
            "--l" => a.cfg.l = num(&mut it, "--l")?,
            "--scan" => a.cfg.scan = string(&mut it, "--scan")?.parse::<CandidateScan>()?,
            "--run-indexing" => {
                a.cfg.run_indexing = string(&mut it, "--run-indexing")?.parse::<RunIndexing>()?
            }
            "--bucket-bits" => a.cfg.bucket_bits = num(&mut it, "--bucket-bits")?,
            "--hot-threshold" => a.cfg.hot_threshold = num(&mut it, "--hot-threshold")?,
            "--l2" => a.cfg.l2 = num(&mut it, "--l2")?,
            "--spill-dir" => {
                let d = PathBuf::from(string(&mut it, "--spill-dir")?);
                if !d.is_dir() {
                    return Err(format!("spill directory {} does not exist", d.display()));
                }
                a.cfg.spill_dir = Some(d);
            }
            "--max-piece-len" => a.cfg.max_piece_len = num(&mut it, "--max-piece-len")?,
            "--max-records" => a.max_records = num(&mut it, "--max-records")?,
            "--synthetic" => a.synthetic = Some((num(&mut it, "--synthetic")?, false)),
            "--synthetic-repeats" => {
                a.synthetic = Some((num(&mut it, "--synthetic-repeats")?, true))
            }
            "--reads" => a.reads = Some(num(&mut it, "--reads")?),
            "--read-len" => {
                let v = string(&mut it, "--read-len")?;
                let (x, y) = v.split_once('-').unwrap_or((&v, &v));
                let bad = || format!("invalid value '{v}' for --read-len (expected A-B)");
                a.read_len = (x.parse().map_err(|_| bad())?, y.parse().map_err(|_| bad())?);
            }
            "--error-rate" => a.error_rate = num(&mut it, "--error-rate")?,
            "--seed" => a.seed = num(&mut it, "--seed")?,
            "--query" => a.queries.push(PathBuf::from(string(&mut it, "--query")?)),
            "--query-kmers" | "--pos-limit" => a.query_kmers = num(&mut it, "--query-kmers")?,
            "--window" => a.window = num(&mut it, "--window")?,
            "--holdout" => a.holdout = true,
            "--neg" => a.neg = num(&mut it, "--neg")?,
            "--no-oracle" => a.oracle = false,
            "--defrag" => a.defrag = true,
            "--dump" => a.dump = Some(PathBuf::from(string(&mut it, "--dump")?)),
            "--top" => a.top = num(&mut it, "--top")?,
            s if s.starts_with("--") => return Err(format!("unknown option {s}")),
            _ => a.inputs.push(PathBuf::from(arg)),
        }
    }
    if a.inputs.is_empty() && a.synthetic.is_none() {
        return Err("no input: give FASTA/FASTQ files or --synthetic/--synthetic-repeats".into());
    }
    if a.holdout && (a.inputs.is_empty() || a.max_records == u64::MAX) {
        return Err("--holdout needs input files and --max-records".into());
    }
    if a.window <= a.cfg.k {
        return Err(format!("--window must be larger than k = {}", a.cfg.k));
    }
    a.cfg.validate().map_err(|e| e.to_string())?;
    // Fail now rather than after a long insertion pass.
    for p in a.inputs.iter().chain(&a.queries) {
        if p.as_os_str() != "-" && !p.is_file() {
            return Err(format!("input file {} does not exist", p.display()));
        }
    }
    Ok(a)
}

// ---------------------------------------------------------------------------
// Input sources and replayable buffers
// ---------------------------------------------------------------------------

const BATCH_BYTES: usize = 4 << 20;

/// A buffer of input events (sequence chunks and end-of-record marks) that
/// can be replayed into several consumers.
#[derive(Default)]
struct Batch {
    data: Vec<u8>,
    /// `Some(end)`: chunk `data[prev_end..end]`; `None`: end of record.
    events: Vec<Option<usize>>,
}

impl Batch {
    fn push_chunk(&mut self, c: &[u8]) {
        if !c.is_empty() {
            self.data.extend_from_slice(c);
            self.events.push(Some(self.data.len()));
        }
    }
    fn push_end(&mut self) {
        self.events.push(None);
    }
    fn bytes(&self) -> usize {
        self.data.len()
    }
    fn clear(&mut self) {
        self.data.clear();
        self.events.clear();
    }
    fn replay(&self, mut f: impl FnMut(Option<&[u8]>)) {
        let mut start = 0;
        for e in &self.events {
            match *e {
                Some(end) => {
                    f(Some(&self.data[start..end]));
                    start = end;
                }
                None => f(None),
            }
        }
    }
}

/// A source of records that can be read several times.
enum Source {
    Memory(Vec<Vec<u8>>),
    /// Files; per file, `skip` records are skipped and at most `max` used.
    Files {
        paths: Vec<PathBuf>,
        skip: u64,
        max: u64,
    },
}

impl Source {
    /// Stream every record as chunks: `Some(chunk)` for sequence data and
    /// `None` at the end of each record. File records are never buffered
    /// whole, so a chromosome costs one input buffer, not its length.
    fn for_each_chunk(&self, mut f: impl FnMut(Option<&[u8]>)) {
        match self {
            Source::Memory(v) => {
                for r in v {
                    f(Some(r));
                    f(None);
                }
            }
            Source::Files { paths, skip, max } => {
                for p in paths {
                    let mut rd = fastx::open(p).unwrap_or_else(|e| die(e));
                    let mut skipped = 0u64;
                    while skipped < *skip
                        && rd.next_record_chunks(|_| {}).unwrap_or_else(|e| die(e))
                    {
                        skipped += 1;
                    }
                    let mut n = 0u64;
                    while n < *max
                        && rd
                            .next_record_chunks(|c| f(Some(c)))
                            .unwrap_or_else(|e| die(e))
                    {
                        f(None);
                        n += 1;
                    }
                }
            }
        }
    }

    /// (records, bases) in one pass.
    fn count(&self) -> (u64, u64) {
        let (mut r, mut b) = (0u64, 0u64);
        self.for_each_chunk(|c| match c {
            Some(c) => b += c.len() as u64,
            None => r += 1,
        });
        (r, b)
    }
}

// ---------------------------------------------------------------------------
// Uniform query workloads
// ---------------------------------------------------------------------------

/// Query workload: sequence windows and, in the same order, their k-mers.
struct Workload {
    seq: Batch,
    kmers: Vec<u64>,
    windows: u64,
    /// Distance between window starts (genome mode), or 0 (record mode).
    #[cfg_attr(not(test), allow(dead_code))]
    stride: u64,
    description: String,
}

/// Sample about `target` k-mers from `src`, spread uniformly over it.
///
/// * Records shorter than `window` (reads): every `s`-th record, whole.
/// * Longer records (genomes): windows of `window` consecutive bases whose
///   starts are a fixed stride apart over the concatenated input; a window
///   ends early at the end of its record.
///
/// Each window is a separate sequence (the streaming query restarts cold).
fn build_workload(src: &Source, k: usize, target: u64, window: usize) -> Workload {
    let (records, bases) = src.count();
    let mean = bases as f64 / records.max(1) as f64;
    let mut seq = Batch::default();
    let mut windows = 0u64;
    let mut stride_used = 0u64;
    let description;
    if mean <= window as f64 {
        let per = (mean - k as f64 + 1.0).max(1.0);
        let every = ((records as f64 * per / target as f64).floor() as u64).max(1);
        let mut r = 0u64;
        let mut take = false;
        src.for_each_chunk(|c| match c {
            Some(c) => {
                if r % every == 0 {
                    seq.push_chunk(c);
                    take = true;
                }
            }
            None => {
                if take {
                    seq.push_end();
                    windows += 1;
                    take = false;
                }
                r += 1;
            }
        });
        description = format!("1 record in {every} of {records} (mean length {mean:.0}), whole");
    } else {
        let per = (window - k + 1) as f64;
        let n_windows = (target as f64 / per).ceil().max(1.0);
        let stride = ((bases as f64 / n_windows).floor() as u64).max(window as u64);
        stride_used = stride;
        let mut pos = 0u64; // bases seen so far
        let mut next_start = 0u64;
        let mut remaining = 0usize; // > 0 while a window is open
        src.for_each_chunk(|c| match c {
            Some(mut c) => {
                while !c.is_empty() {
                    if remaining > 0 {
                        let take = remaining.min(c.len());
                        seq.push_chunk(&c[..take]);
                        pos += take as u64;
                        remaining -= take;
                        c = &c[take..];
                        if remaining == 0 {
                            seq.push_end();
                            windows += 1;
                            next_start += stride;
                        }
                    } else if pos + c.len() as u64 <= next_start {
                        pos += c.len() as u64;
                        c = &[];
                    } else {
                        let skip = (next_start - pos) as usize;
                        pos += skip as u64;
                        c = &c[skip..];
                        remaining = window;
                    }
                }
            }
            None => {
                if remaining > 0 {
                    // Window cut short by the end of its record.
                    seq.push_end();
                    windows += 1;
                    remaining = 0;
                    next_start += stride;
                }
                // Keep starts ahead of the current position.
                while next_start < pos {
                    next_start += stride;
                }
            }
        });
        description = format!(
            "windows of ≤{window} bases every {stride} bases over {bases} bases ({records} records)"
        );
    }
    let mut kmers = Vec::new();
    let mut ks = KmerStream::new(k);
    seq.replay(|c| match c {
        Some(c) => ks.push(c, |x| kmers.push(x)),
        None => ks.end(),
    });
    Workload {
        seq,
        kmers,
        windows,
        stride: stride_used,
        description,
    }
}

// ---------------------------------------------------------------------------
// Measurements
// ---------------------------------------------------------------------------

fn rate(n: u64, d: Duration) -> String {
    let s = d.as_secs_f64();
    if s == 0.0 {
        return "n/a".into();
    }
    let r = n as f64 / s;
    if r >= 1e6 {
        format!("{:.2} M/s", r / 1e6)
    } else {
        format!("{:.1} k/s", r / 1e3)
    }
}

fn mib(b: usize) -> String {
    format!("{:.1} MiB", b as f64 / (1 << 20) as f64)
}

/// (VmRSS, VmHWM) in bytes.
fn rss() -> (usize, usize) {
    let s = std::fs::read_to_string("/proc/self/status").unwrap_or_default();
    let get = |key: &str| -> usize {
        s.lines()
            .find(|l| l.starts_with(key))
            .and_then(|l| l.split_whitespace().nth(1))
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(0)
            * 1024
    };
    (get("VmRSS:"), get("VmHWM:"))
}

fn timed<T>(f: impl FnOnce() -> T) -> (T, Duration) {
    let t0 = Instant::now();
    let v = f();
    (v, t0.elapsed())
}

fn count_mismatches(a: &[bool], b: &[bool]) -> u64 {
    assert_eq!(a.len(), b.len());
    a.iter().zip(b).filter(|(x, y)| x != y).count() as u64
}

/// Time every access pattern on the same k-mers. Returns the mismatches.
fn write_u64s(path: &std::path::Path, xs: &[u64]) {
    let mut bytes = Vec::with_capacity(xs.len() * 8);
    for x in xs {
        bytes.extend_from_slice(&x.to_le_bytes());
    }
    std::fs::write(path, bytes).unwrap_or_else(|e| die(format!("{}: {e}", path.display())));
}

/// Write a query workload for the kache-hash driver: `q<id>.fa` (one record
/// per window), `q<id>.inorder.u64` and `q<id>.shuffled.u64` (canonical
/// 2-bit k-mers, little-endian), and a manifest line with the expected hits.
fn dump_workload(
    dir: &std::path::Path,
    id: usize,
    label: &str,
    w: &Workload,
    shuffled: &[u64],
    hits: u64,
) {
    use std::io::Write;
    let mut fa = Vec::with_capacity(w.seq.bytes() + 16 * w.windows as usize);
    let mut n = 0u64;
    let mut open = false;
    w.seq.replay(|c| match c {
        Some(c) => {
            if !open {
                fa.extend_from_slice(format!(">w{n}\n").as_bytes());
                open = true;
            }
            fa.extend_from_slice(c);
        }
        None => {
            if open {
                fa.push(b'\n');
                open = false;
            }
            n += 1;
        }
    });
    std::fs::write(dir.join(format!("q{id}.fa")), fa).unwrap_or_else(|e| die(e));
    write_u64s(&dir.join(format!("q{id}.inorder.u64")), &w.kmers);
    write_u64s(&dir.join(format!("q{id}.shuffled.u64")), shuffled);
    let mut m = std::fs::OpenOptions::new()
        .append(true)
        .open(dir.join("manifest.tsv"))
        .unwrap_or_else(|e| die(e));
    writeln!(m, "query\t{id}\t{label}\t{}\t{hits}", w.kmers.len()).unwrap_or_else(|e| die(e));
}

fn measure_queries(
    label: &str,
    w: &Workload,
    table: &KmerTable,
    oracle: Option<&Oracle>,
    rng: &mut Rng,
    dump: Option<(&std::path::Path, usize)>,
) -> u64 {
    let k = table.k();
    let n = w.kmers.len() as u64;
    let mut shuffled = w.kmers.clone();
    rng.shuffle(&mut shuffled);

    // Pre-parsed k-mers: in sequence order and shuffled.
    let lookups = |xs: &[u64]| {
        let mut st = QueryStats::default();
        let mut ans = Vec::with_capacity(xs.len());
        for &x in xs {
            let mut tr = QueryTrace::default();
            let hit = table.contains_traced(x, &mut tr);
            st.record(hit, &tr);
            ans.push(hit);
        }
        (st, ans)
    };
    let ((st_in, a_in), t_in) = timed(|| lookups(&w.kmers));
    let ((st_rnd, a_rnd), t_rnd) = timed(|| lookups(&shuffled));
    let oracle_lookups =
        |o: &Oracle, xs: &[u64]| timed(|| xs.iter().map(|&x| o.contains(x)).collect::<Vec<_>>());
    let or_in = oracle.map(|o| oracle_lookups(o, &w.kmers));
    let or_rnd = oracle.map(|o| oracle_lookups(o, &shuffled));

    // From the sequence (parsing included), same k-mers in the same order.
    let (parse_sum, t_parse) = timed(|| {
        let mut ks = KmerStream::new(k);
        let mut acc = 0u64;
        w.seq.replay(|c| match c {
            Some(c) => ks.push(c, |x| acc = acc.wrapping_add(x)),
            None => ks.end(),
        });
        acc
    });
    std::hint::black_box(parse_sum);
    let (a_seq, t_seq) = timed(|| {
        let mut ks = KmerStream::new(k);
        let mut ans = Vec::with_capacity(n as usize);
        w.seq.replay(|c| match c {
            Some(c) => ks.push(c, |x| ans.push(table.contains(x))),
            None => ks.end(),
        });
        ans
    });
    let mut sq = table.query_stream();
    let (a_str, t_str) = timed(|| {
        let mut ans = Vec::with_capacity(n as usize);
        w.seq.replay(|c| match c {
            Some(c) => sq.push(c, |_, p| ans.push(p)),
            None => sq.end(),
        });
        ans
    });
    let st_str = *sq.stats();
    let or_seq = oracle.map(|o| {
        timed(|| {
            let mut ks = KmerStream::new(k);
            let mut ans = Vec::with_capacity(n as usize);
            w.seq.replay(|c| match c {
                Some(c) => ks.push(c, |x| ans.push(o.contains(x))),
                None => ks.end(),
            });
            ans
        })
    });

    // Correctness: everything against the oracle, or against each other.
    let truth_in: &[bool] = or_in.as_ref().map_or(&a_in, |(v, _)| v);
    let mm_in = count_mismatches(&a_in, truth_in);
    let mm_rnd = match &or_rnd {
        Some((v, _)) => count_mismatches(&a_rnd, v),
        None => 0,
    };
    let mm_seq = count_mismatches(&a_seq, truth_in);
    let mm_str = count_mismatches(&a_str, truth_in);

    println!("\n-- queries: {label} --");
    println!(
        "workload: {n} k-mers from {} windows; {}; hit rate {:.4}",
        w.windows,
        w.description,
        st_in.hits as f64 / n.max(1) as f64
    );
    println!("pattern (same k-mers)        table                         oracle HashSet<u64>");
    let row = |tag: &str, t: Duration, o: Option<Duration>| {
        let o = o.map_or("n/a".to_string(), |o| {
            format!("{:>8.3} s  {:>11}", o.as_secs_f64(), rate(n, o))
        });
        println!(
            "{tag:<26} {:>8.3} s  {:>11}      {o}",
            t.as_secs_f64(),
            rate(n, t)
        );
    };
    row("kmers-in-order", t_in, or_in.as_ref().map(|x| x.1));
    row("kmers-random-order", t_rnd, or_rnd.as_ref().map(|x| x.1));
    row("seq-parse-only", t_parse, None);
    row("seq-independent", t_seq, or_seq.as_ref().map(|x| x.1));
    row("seq-streaming", t_str, None);
    println!("k-mers in order (pre-parsed):\n{st_in}");
    println!("k-mers in random order (pre-parsed):\n{st_rnd}");
    println!("streaming:\n{st_str}");
    if let Some((dir, id)) = dump {
        let hits = truth_in.iter().filter(|&&h| h).count() as u64;
        dump_workload(dir, id, label, w, &shuffled, hits);
    }
    let total = mm_in + mm_rnd + mm_seq + mm_str;
    println!(
        "  MISMATCHES vs {}       {total} (in-order {mm_in}, random {mm_rnd}, seq {mm_seq}, streaming {mm_str})",
        if oracle.is_some() { "oracle" } else { "table " }
    );
    total
}

/// Start a dump directory: `manifest.tsv` with k and the insertion inputs
/// (file paths with their record cap, or `insert.fa` for in-memory input).
fn dump_inputs(dir: &std::path::Path, k: usize, source: &Source) {
    use std::io::Write;
    std::fs::create_dir_all(dir).unwrap_or_else(|e| die(e));
    let mut m = std::fs::File::create(dir.join("manifest.tsv")).unwrap_or_else(|e| die(e));
    writeln!(m, "k\t{k}").unwrap();
    match source {
        Source::Files { paths, skip, max } => {
            for p in paths {
                let abs = std::fs::canonicalize(p).unwrap_or_else(|_| p.clone());
                writeln!(m, "insert\t{}\t{skip}\t{max}", abs.display()).unwrap();
            }
        }
        Source::Memory(v) => {
            let mut fa = Vec::new();
            for (i, r) in v.iter().enumerate() {
                fa.extend_from_slice(format!(">r{i}\n").as_bytes());
                fa.extend_from_slice(r);
                fa.push(b'\n');
            }
            std::fs::write(dir.join("insert.fa"), fa).unwrap_or_else(|e| die(e));
            let abs = std::fs::canonicalize(dir.join("insert.fa")).unwrap();
            writeln!(m, "insert\t{}\t0\t{}", abs.display(), u64::MAX).unwrap();
        }
    }
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

fn main() {
    let args = parse_args();
    let k = args.cfg.k;
    let mut rng = Rng::new(args.seed);

    let source = if let Some((len, repeats)) = args.synthetic {
        let genome = if repeats {
            repeat_rich_genome(&mut rng, len)
        } else {
            rng.dna(len)
        };
        let reads = match args.reads {
            Some(n) => sample_reads(
                &mut rng,
                &genome,
                n,
                args.read_len.0,
                args.read_len.1,
                args.error_rate,
            ),
            None => vec![genome],
        };
        let mut all = reads;
        if !args.inputs.is_empty() {
            let mut cur = Vec::new();
            Source::Files {
                paths: args.inputs.clone(),
                skip: 0,
                max: args.max_records,
            }
            .for_each_chunk(|c| match c {
                Some(c) => cur.extend_from_slice(c),
                None => all.push(std::mem::take(&mut cur)),
            });
        }
        Source::Memory(all)
    } else {
        Source::Files {
            paths: args.inputs.clone(),
            skip: 0,
            max: args.max_records,
        }
    };

    println!("== unihash-bench ==");
    println!("config: {:?}", args.cfg);
    println!(
        "input:  {}",
        match (&source, args.synthetic) {
            (Source::Memory(v), Some((len, rep))) => format!(
                "synthetic {} genome, {len} bp, {} records",
                if rep { "repeat-rich" } else { "uniform" },
                v.len()
            ),
            _ => format!("{:?}", args.inputs),
        }
    );

    // ------------------------------------------------- insertion, table only
    // Input is decoded into ~4 MiB batches, each inserted under one timer.
    let insert_pass = |sink: &mut dyn FnMut(Option<&[u8]>)| -> (Duration, u64, u64) {
        let mut t = Duration::ZERO;
        let (mut records, mut bases) = (0u64, 0u64);
        let mut batch = Batch::default();
        let mut run = |batch: &mut Batch, t: &mut Duration| {
            let t0 = Instant::now();
            batch.replay(&mut *sink);
            *t += t0.elapsed();
            batch.clear();
        };
        source.for_each_chunk(|c| {
            match c {
                Some(c) => {
                    bases += c.len() as u64;
                    batch.push_chunk(c);
                }
                None => {
                    records += 1;
                    batch.push_end();
                }
            }
            if batch.bytes() >= BATCH_BYTES {
                run(&mut batch, &mut t);
            }
        });
        run(&mut batch, &mut t);
        (t, records, bases)
    };

    let (rss_base, _) = rss();
    let heap_base = counting::current();
    counting::reset_peak();
    let mut table = KmerTable::new(args.cfg.clone()).unwrap_or_else(|e| die(e));
    let (t_table, records, bases) = insert_pass(&mut |c| match c {
        Some(c) => table.insert_chunk(c),
        None => table.end_read(),
    });
    let heap_peak_insert = counting::peak() - heap_base;
    let heap_live_before_shrink = counting::current() - heap_base;
    let insert_buffers = table.insert_buffer_bytes();
    table.shrink_to_fit();
    let defrag = args.defrag.then(|| {
        let before = counting::current() - heap_base;
        counting::reset_peak();
        let t0 = Instant::now();
        let rep = table.defrag();
        let d = t0.elapsed();
        table.shrink_to_fit();
        let peak = counting::peak() - heap_base;
        (rep, d, before, peak)
    });
    let heap_table = counting::current() - heap_base;
    let (rss_table, hwm_table) = rss();
    let kmers_in = table.insert_stats().kmers_processed;
    let st = table.stats();

    println!("\n-- insertion --");
    println!("records {records}, bases {bases}, k-mers {kmers_in}");
    println!(
        "table   {:>9.3} s  {}",
        t_table.as_secs_f64(),
        rate(kmers_in, t_table)
    );

    // ------------------------------------------------- insertion, oracle
    let mut oracle = args.oracle.then(|| Oracle::new(k));
    let heap_before_oracle = counting::current();
    if let Some(o) = oracle.as_mut() {
        let (t_oracle, _, _) = insert_pass(&mut |c| match c {
            Some(c) => o.insert_chunk(c),
            None => o.end_read(),
        });
        o.shrink_to_fit();
        println!(
            "oracle  {:>9.3} s  {}",
            t_oracle.as_secs_f64(),
            rate(kmers_in, t_oracle)
        );
    }
    let heap_oracle = counting::current() - heap_before_oracle;

    if let Some((rep, d, before, peak)) = &defrag {
        println!(
            "defrag  {:>9.3} s  heap {before} → {heap_table} B (peak {peak} B); {rep}",
            d.as_secs_f64()
        );
    }

    println!("\n-- table stats --\n{st}");

    // ------------------------------------------------- memory
    // Without the oracle: eager run indexing never stores a k-mer twice;
    // otherwise count the distinct stored k-mers (a temporary hash set,
    // after the table measurements above).
    let distinct = match &oracle {
        Some(o) => o.len() as u64,
        None if args.cfg.run_indexing == RunIndexing::Eager => st.stored_kmers,
        None => table.distinct_kmers(),
    };
    let per = |b: usize| b as f64 / distinct.max(1) as f64;
    let blocks = st.buckets as usize + st.overflow_blocks as usize;
    println!("\n-- memory: table alone (measured before the oracle is built) --");
    println!(
        "distinct k-mers              {distinct}{}",
        if oracle.is_none() {
            " (counted from the pieces)"
        } else {
            ""
        }
    );
    println!(
        "redundant stored k-mers      {}",
        st.stored_kmers.saturating_sub(distinct)
    );
    println!(
        "pieces (2-bit + offsets)     {:>12} B  {:6.3} B/k-mer  ({} bases, {} pieces)",
        st.piece_bytes_reserved,
        per(st.piece_bytes_reserved),
        st.total_bases,
        st.pieces
    );
    println!(
        "index (flat blocks)          {:>12} B  {:6.3} B/k-mer  ({blocks} × {BLOCK_BYTES} B blocks: {} heads + {} overflow; {:.1} B/entry)",
        st.index_bytes_reserved,
        per(st.index_bytes_reserved),
        st.buckets,
        st.overflow_blocks,
        st.index_bytes_reserved as f64 / st.entries.max(1) as f64
    );
    println!(
        "layer 2 index (ℓ2 = {:>2})      {:>12} B  {:6.3} B/k-mer  ({} entries, {} buckets; {} hot ℓ-minimizers)",
        st.l2,
        st.index2_bytes_reserved,
        per(st.index2_bytes_reserved),
        st.entries2,
        st.buckets2,
        st.insert.hot1
    );
    println!(
        "layer 3 k-mer set            {:>12} B  {:6.3} B/k-mer  ({} k-mers; {} hot ℓ2-minimizers)",
        st.layer3_bytes_reserved,
        per(st.layer3_bytes_reserved),
        st.layer3_kmers,
        st.insert.hot2
    );
    println!(
        "changed-piece bitset         {:>12} B  {:6.3} B/k-mer  (pieces created or cut since the last defrag)",
        st.dirty_bytes_reserved,
        per(st.dirty_bytes_reserved)
    );
    let accounted = st.piece_bytes_reserved
        + st.index_bytes_reserved
        + st.index2_bytes_reserved
        + st.layer3_bytes_reserved
        + st.dirty_bytes_reserved;
    println!(
        "table total, accounted       {:>12} B  {:6.3} B/k-mer",
        accounted,
        per(accounted)
    );
    println!(
        "table heap, allocator        {:>12} B  {:6.3} B/k-mer  (live after shrink_to_fit; accounted − allocator = {} B)",
        heap_table,
        per(heap_table),
        accounted as i64 - heap_table as i64
    );
    println!(
        "heap peak during insertion   {:>12} B  {:6.3} B/k-mer  (index growth, Vec slack, I/O and insertion buffers; {} B live before shrink)",
        heap_peak_insert,
        per(heap_peak_insert),
        heap_live_before_shrink
    );
    println!(
        "insertion buffers            {:>12} B  (bounded by ~14 × max_piece_len, not by input length)",
        insert_buffers
    );
    println!(
        "index growth spill           {} rebuilds; {} B via disk, {} B via memory{}",
        st.spill.rebuilds,
        st.spill.bytes_to_disk,
        st.spill.bytes_in_memory,
        match &args.cfg.spill_dir {
            Some(d) => format!(" (spill dir {})", d.display()),
            None => String::new(),
        }
    );
    println!(
        "RSS: process before table {}, after table {} (+{}), peak {}",
        mib(rss_base),
        mib(rss_table),
        mib(rss_table.saturating_sub(rss_base)),
        mib(hwm_table)
    );
    if let Some(o) = &oracle {
        println!(
            "oracle HashSet<u64>, allocator {:>10} B  {:6.3} B/k-mer  (estimate from capacity {} B)",
            heap_oracle,
            per(heap_oracle),
            o.heap_bytes()
        );
        println!(
            "ideal u64 table @7/8 load    {:>12} B  {:6.3} B/k-mer",
            o.ideal_bytes(),
            per(o.ideal_bytes())
        );
    }
    println!(
        "mean piece length            {:.1} bases (index entries / piece {:.2})",
        st.mean_piece_len(),
        st.entries as f64 / st.pieces.max(1) as f64
    );

    // ------------------------------------------------- query workloads
    let mut mismatches = 0u64;
    let mut sets: Vec<(String, Source)> = Vec::new();
    sets.push((
        "inserted input (positives)".into(),
        match &source {
            Source::Memory(v) => Source::Memory(v.clone()),
            Source::Files { paths, skip, max } => Source::Files {
                paths: paths.clone(),
                skip: *skip,
                max: *max,
            },
        },
    ));
    if args.holdout {
        sets.push((
            format!(
                "held-out records {}..{} of each input file",
                args.max_records + 1,
                2 * args.max_records
            ),
            Source::Files {
                paths: args.inputs.clone(),
                skip: args.max_records,
                max: args.max_records,
            },
        ));
    }
    for q in &args.queries {
        sets.push((
            format!("query set {}", q.display()),
            Source::Files {
                paths: vec![q.clone()],
                skip: 0,
                max: args.max_records,
            },
        ));
    }
    if let Some(dir) = &args.dump {
        dump_inputs(dir, k, &source);
    }
    for (i, (label, src)) in sets.iter().enumerate() {
        let w = build_workload(src, k, args.query_kmers, args.window);
        let dump = args.dump.as_deref().map(|d| (d, i));
        mismatches += measure_queries(label, &w, &table, oracle.as_ref(), &mut rng, dump);
    }

    // ------------------------------------------------- negative queries
    let mask = kmer_mask(k);
    let mut negs = Vec::with_capacity(args.neg as usize);
    while (negs.len() as u64) < args.neg {
        let x = rng.next_u64() & mask;
        let present = match &oracle {
            Some(o) => o.contains(x),
            None => false,
        };
        if !present {
            negs.push(x);
        }
    }
    if let Some(dir) = &args.dump {
        use std::io::Write;
        write_u64s(&dir.join("neg.u64"), &negs);
        let mut m = std::fs::OpenOptions::new()
            .append(true)
            .open(dir.join("manifest.tsv"))
            .unwrap_or_else(|e| die(e));
        writeln!(m, "neg\t{}\t0", negs.len()).unwrap_or_else(|e| die(e));
    }
    let mut neg = QueryStats::default();
    let ((), t_neg) = timed(|| {
        for &x in &negs {
            let mut tr = QueryTrace::default();
            let hit = table.contains_traced(x, &mut tr);
            neg.record(hit, &tr);
        }
    });
    println!("\n-- negative queries (uniform random k-mers absent from the set) --");
    println!(
        "table   {:>9.3} s  {}",
        t_neg.as_secs_f64(),
        rate(neg.queries, t_neg)
    );
    if let Some(o) = &oracle {
        let (c, t) = timed(|| negs.iter().filter(|&&x| o.contains(x)).count());
        println!(
            "oracle  {:>9.3} s  {}",
            t.as_secs_f64(),
            rate(negs.len() as u64, t)
        );
        assert_eq!(c, 0);
    }
    println!("{neg}");
    if oracle.is_some() {
        println!("  FALSE POSITIVES              {}", neg.hits);
    }

    // ------------------------------------------------- minimizer skew
    let mut occ: HashMap<u64, u32> = HashMap::new();
    for m in table.entry_minimizers() {
        *occ.entry(m).or_insert(0) += 1;
    }
    let mut hist = [0u64; 33];
    for &c in occ.values() {
        hist[(32 - c.leading_zeros()) as usize] += 1;
    }
    println!("\n-- minimizer skew --");
    println!("distinct minimizers {}, entries {}", occ.len(), st.entries);
    println!("entries/minimizer   #minimizers");
    for (b, &n) in hist.iter().enumerate().skip(1) {
        if n > 0 {
            let lo = 1u64 << (b - 1);
            let hi = (1u64 << b) - 1;
            println!("  {:>7}-{:<9} {n}", lo, hi);
        }
    }
    let mut hot: Vec<(u64, u32)> = occ.into_iter().collect();
    hot.sort_unstable_by_key(|h| std::cmp::Reverse(h.1));
    println!("hottest minimizers:");
    for (m, c) in hot.iter().take(args.top) {
        println!(
            "  {}  {c}",
            decode_kmer(lmer_of_order_value(*m), args.cfg.l)
        );
    }
    let mut bmax = 0usize;
    let mut bempty = 0usize;
    for b in 0..table.index().num_buckets() {
        let n = table.index().bucket_len(b);
        bmax = bmax.max(n);
        bempty += (n == 0) as usize;
    }
    println!(
        "buckets: {} total ({} slots each), {} empty, max chain occupancy {}, overflow blocks {}, grew {}×",
        table.index().num_buckets(),
        BLOCK_SLOTS,
        bempty,
        bmax,
        table.index().overflow_blocks(),
        table.index().grows()
    );

    let (_, hwm) = rss();
    println!(
        "\npeak RSS {} (whole run, includes the oracle and query workloads)",
        mib(hwm)
    );
    if mismatches > 0 || (oracle.is_some() && neg.hits > 0) {
        eprintln!("ERROR: table disagrees with the oracle");
        exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(v: &[&str]) -> Result<Args, String> {
        parse_args_from(v.iter().map(|s| s.to_string()))
    }

    #[test]
    fn paths_keep_underscores_and_numbers_accept_separators() {
        let dir = std::env::temp_dir().join(format!("unihash_args_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let reads = dir.join("SRR5833294.fastq");
        let query = dir.join("SRR5901135_1.fastq");
        std::fs::write(&reads, "@r\nACGT\n+\nIIII\n").unwrap();
        std::fs::write(&query, "@r\nACGT\n+\nIIII\n").unwrap();
        let a = args(&[
            reads.to_str().unwrap(),
            "--query",
            query.to_str().unwrap(),
            "--neg",
            "1_000",
            "--max-records",
            "2_000_000",
            "--holdout",
        ])
        .unwrap();
        assert_eq!(a.queries, vec![query.clone()]);
        assert_eq!(a.inputs, vec![reads.clone()]);
        assert_eq!(a.neg, 1000);
        assert_eq!(a.max_records, 2_000_000);
        assert!(a.holdout);
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn missing_files_are_reported_before_running() {
        let e = args(&[
            "--synthetic",
            "1000",
            "--query",
            "/nonexistent/SRR5901135_1.fastq.gz",
        ])
        .err()
        .unwrap();
        assert!(e.contains("/nonexistent/SRR5901135_1.fastq.gz"), "{e}");
        assert!(args(&["--k", "40", "--synthetic", "100"]).is_err());
        assert!(args(&["--bogus"]).is_err());
        assert!(args(&["--synthetic", "100", "--holdout"]).is_err());
    }

    fn workload_of(records: Vec<Vec<u8>>, target: u64, window: usize) -> Workload {
        build_workload(&Source::Memory(records), 31, target, window)
    }

    #[test]
    fn genome_workload_spreads_windows_over_the_whole_input() {
        let mut rng = Rng::new(5);
        let genome = rng.dna(1_000_000);
        let w = workload_of(vec![genome.clone()], 50_000, 1000);
        // ~50 k k-mers in ~52 windows of 1000 bases (970 k-mers each)
        assert!(
            w.kmers.len() >= 45_000 && w.kmers.len() <= 55_000,
            "{}",
            w.kmers.len()
        );
        assert_eq!(
            w.windows as usize,
            w.seq.events.iter().filter(|e| e.is_none()).count()
        );
        // Window i is genome[i·stride .. i·stride + 1000] (the last one may
        // be cut by the end of the record); together they reach the end.
        let stride = w.stride as usize;
        let mut expect = Vec::new();
        for i in 0..w.windows as usize {
            let s = i * stride;
            expect.extend_from_slice(&genome[s..(s + 1000).min(genome.len())]);
        }
        assert_eq!(w.seq.data, expect);
        assert!((w.windows as usize - 1) * stride > 950_000);
    }

    #[test]
    fn read_workload_takes_whole_records_at_a_stride() {
        let mut rng = Rng::new(6);
        let reads: Vec<Vec<u8>> = (0..10_000).map(|_| rng.dna(100)).collect();
        let w = workload_of(reads.clone(), 70_000, 1000);
        // 70 k-mers per read → every 10th read
        assert_eq!(w.windows, 1000);
        assert_eq!(w.kmers.len(), 70_000);
        assert_eq!(&w.seq.data[..100], &reads[0][..]);
        assert_eq!(&w.seq.data[100..200], &reads[10][..]);
    }
}
