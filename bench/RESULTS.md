# Benchmark results

k = 31, ℓ = 17, release build, single thread, 1 M negative queries per run, and at
most 20 M positive queries. Raw outputs are in `results/`. Every run reported
**0 mismatches** against the `HashSet<u64>` oracle and **0 false positives**.

## Optimization round 2: insertion shortcuts, 64 B blocks, in-place defrag, concurrent streaming (2026-10-29, sandbox, 2 cores)

Against the previous commit on the same machine; `unihash-bench` (one
thread) and `unihash-scale --capacity --reps 3 --check` (1 and 2 threads).
Every run: 0 mismatches against the oracle, every concurrent table equal to
the single-threaded one with an exact index. Synthetic reads:
`--synthetic 4000000 --reads 400000` (38 M k-mers, 8.1 M distinct, 12×).

| | before | after |
|---|---|---|
| `KmerTable` insertion, chr21 / reads (min of 5) | 5.68 / 3.91 s | 5.03 / 2.04 s |
| `KmerTable` after `shrink_to_fit`, chr21 / E. coli / reads | 2.020 / 1.819 / 2.715 B/k-mer | 1.927 / 1.756 / 2.600 |
| concurrent insertion, chr21, 1 / 2 threads | 6.73 / 3.46 s | 4.93 / 2.62 s |
| concurrent insertion, reads, 1 / 2 threads | 5.03 / 2.90 s | 2.30 / 1.35 s |
| concurrent table (`with_capacity`), chr21 / reads | 2.658 / 2.908 B/k-mer | 2.450 / 2.786 |
| after `into_table`, chr21 / reads | 2.007 / 2.697 B/k-mer | 1.866 / 2.503 |
| concurrent deletion (every 4th item), chr21, 1 / 2 threads | 2.75 / 1.36 s | 2.58 / 1.32 s |
| concurrent `defrag`, chr21 after deletion (every 4th item) | 2.97 s (conversion) | 1.42 s, same pieces |
| concurrent `defrag_incremental`, same | 1.53 s (conversion) | 0.29 s |
| concurrent streaming queries, chr21 / reads, 1 thread | (only on `KmerTable`) | 34.3 / 25.2 M k-mers/s |

* **Insertion lookups** (`KmerTable` and `Writer`): hit extension (one base
  read when the previous k-mer was present) and a candidate cache per
  minimizer (the bucket is scanned once while successive k-mers share their
  minimizer; the run's own new entries are added to it). On chr21 4.9 M of
  the 5.5 M present k-mers are extensions and 30 M of 35 M lookups reuse
  candidates; on the reads 28.6 M of 38 M k-mers are extensions. The writer
  re-reads the piece record before extending (deletions shrink pieces) and
  uses cached candidates only under the shard lock it collected them under.
  Tests: `tests/m13_insert_lookups.rs`. Of 9 injected faults, 7 fail them
  (no cache reset between chunks or after a deferred run is indexed, own
  entries not added, no base check on extension; in the writer: no
  piece-record check, cache used across lock acquisitions, own entries not
  added). Two guards are not exercised: the cache reset after a promotion
  (unreachable today: a hot minimizer is looked up before the cache), and the
  writer's check that a cursor's bases are the hit's, which matters only if
  the piece is emptied and its id reused between the lookup and the cursor.
* **64 B index blocks** (7 entries, 3-byte offsets, one aligned cache line):
  per slot 9.1 B instead of 10 B (−8.6 %). On one input, though, the size is
  set by where its entry count falls between two ×1.5 growth steps (heads
  57–85 % full), and 7 slots per block move those steps: over 40 genome sizes
  (`--synthetic`, 2–8 Mbp, 3.5 % apart, both builds with the insertion
  shortcuts) the 64 B index is 20 % smaller for 26 sizes and 20 % larger for
  14, −6.0 % on average (whole table −5.4 %). chr21 and the read set fall on
  the larger side (index +19 %). So `shrink_to_fit` now re-sizes an index to
  ≈ 80 % when that saves 1/16 of its heads (0.81 s on chr21, 2.30 →
  1.93 B/k-mer), and `into_table` sizes the merged index for its exact entry
  count; the table above is after these. Query speed showed no consistent
  change (within this machine's ±15 % noise between runs). A ×1.25 growth
  factor was tried instead, to narrow the 57–85 % range: chr21 insertion
  5.0 → 6.8 s (twice as many re-placements), not kept.
* **In-place defragmentation of `ConcurrentKmerTable`**: `defrag` plans the
  same chains as `KmerTable::defrag` (shared code), writes them into a new
  segmented base array, and re-indexes every piece from all cores into shards
  sized for 80 % load; hot minimizers stay hot. `defrag_incremental` rewrites
  the base array without dead bases (ids and index entries unchanged), drops
  index orphans, and re-sizes shards that would save 1/16 of their heads at
  80 % load (all cores). Tests (`tests/m14_concurrent_defrag.rs`): the same
  pieces as through a conversion, base by base, then concurrent writers and
  deleters on the defragmented table; 7 injected faults (overlap copied twice,
  reverse complement, word-copy offset, hot marks, frontier, stale chunk tails,
  shared bases) all caught. `defrag_incremental` leaves the table 3 % larger
  than the conversion did (2.32 vs 2.26 B/k-mer on chr21) for 5× less time;
  after `defrag` 2.24 vs 2.22.
* **`ConcurrentStreamingQuery`** (`ConcurrentKmerTable::query_stream`): hit
  extension with the piece record re-read, candidates cached while it holds the
  minimizer's shard read lock without interruption (within one `push`), the
  answers of a chunk handed to the callback after the lock is released (so it
  may insert into the table). 61–102 % of `KmerTable`'s streaming speed
  (chr21, 1 / 2 threads: 34 / 42 vs 46 / 69 M k-mers/s; reads: 25 / 37 vs
  25 / 45), 13–28× the rate of independent concurrent lookups of positives.
  Tests (`tests/m15_concurrent_stream.rs`):
  equal to `contains` in every configuration, exact across inserts, deletions
  and defragmentation between pushes, readers during insertion and deletion,
  a writing callback; dropping the cache reset on release or the
  piece-record check fails them, calling back under the lock deadlocks
  (caught by the 300 s timeout).
* Found on the way: adding a generic call of the index rebuild (in
  `shrink_to_fit`) changed how LLVM compiled the insertion path and cost 15 %
  of chr21 insertion, although the new code ran after it. The rebuild is now
  one non-generic, never-inlined function taking `&dyn Fn`.

## Index reads without locks: tried, not kept (2026-10-29)

A version of the concurrent table whose queries and deletion lookups read the
index shards without locks was built and measured, then dropped. Shards held
the single-threaded block layout in atomic words; additions were published by
a release store of the block length, other changes ran under a per-shard
seqlock, and grown shards were swapped in with `arc-swap`.

* **Cost per lookup.** In safe Rust the reader reaches the blocks through
  `arc-swap` → a storage header → the blocks: one more dependent memory access
  than the read lock's word → blocks (an uncontended `arc-swap` load costs the
  same 15 ns as an uncontended read lock). On 2 cores, random positive queries
  were 20–23 % slower with 1 thread (chr21 1.38 → 1.06 M/s, synthetic reads
  2.21 → 1.74) and 17–18 % slower with 2; negative queries were equal or faster
  (reads, 2 threads: 6.78 → 8.49 M/s); insertion and deletion within noise.
* **What it could save.** Read locks write the shard's lock word on every
  lookup; on the laptop that held concurrent queries to 0.89 (chr21) and 0.92
  (SRR5833294) of `KmerTable` at 8 threads (`results/ab/*_prev.txt`): at most
  about 10 % to gain, less than the lock-free path costs. The laptop A/B
  (`results/ab/`) is not conclusive on its own: its baseline predates round 1
  below, and the machine's speed changed by up to 2× between the runs (same
  `KmerTable` code: 0.94 vs 1.83 M positive queries/s). In the lock-free run,
  concurrent queries reached 0.53–0.57 of `KmerTable` at 8 threads.
* Writers looking k-mers up without the lock first (after a present k-mer):
  4 % slower on chr21, no gain on reads.
* Kept: `tests/m12_reader_visibility.rs`. The lock-free version had a race the
  locks hide (a cut adds its right part's items, then shrinks the old piece; a
  lookup between the two could miss a surviving k-mer); these tests, which cut
  a piece from the front one k-mer at a time and make minimizers hot while
  readers check, caught it and each of three injected faults 10 times out of
  10. They pass unchanged on the locked table.
* Two lessons for atomic index blocks, should they come back: store piece ids
  and offsets as `u32` atomics (writing half of a 64-bit word needs a read of
  it, a cache miss plain stores never wait for), and use `#[repr(C)]` (the
  default layout moved the length word to the block's end, onto a second cache
  line: +8 % insertion time).

## Optimization round 1: parallel conversions, `with_capacity` for reads (2026-10-29, sandbox)

* **`into_table`** recomputed every index entry's minimizer from the pieces and
  inserted it, one thread, about 320 ns per entry (two random reads in the
  pieces, two in the new index). Minimizers are now computed by all cores, chunk
  by chunk, and each chunk's head blocks are read in batches of 32 before
  inserting (about 20 % faster insertion). **`from_table`** grew every shard
  from its initial size while filling it, one thread; shards are now sized for
  their share first and filled by all cores.

  | 2 cores | `into_table` before → after | `from_table` before → after |
  |---|---|---|
  | chr21 (4.3 M entries) | 1.05 → 0.69 s | 3.02 → 1.15 s |
  | synthetic reads (4.8 M entries) | 1.76 → 1.07 s | 4.40 → 1.38 s |

  On chr21 the remaining merge time is 0.27 s of minimizers (all cores) and
  0.31 s of insertion (one thread, 72 ns per entry).
* **`with_capacity`** assumed 0.125 index entries per k-mer (2 / (w + 1)); read
  sets need 0.151 (about one more entry per piece, and pieces are short), so
  every shard grew once near the end, under its lock, and ended 1.5× too large.
  It now sizes for 2.5 / (w + 1) at 80 % head load. Synthetic reads (1.6 M):
  6521 → 2892 shard growths (the rest are tiny layer-2 shards), growth time 1.1
  → 0.08 s, concurrent table 3.05 → 3.13 B/k-mer. chr21: 2.41 → 2.66 B/k-mer
  while concurrent (heads at 66 % load); `into_table` gives the same 2.007.
  One writer is now as fast as `KmerTable` on chr21 (1.01×).
* New test: conversions large enough to use all cores, both ways, exact.

## Concurrent insertion and deletion on the laptop (2026-10-28, second run)

`CAPACITY=1 tools/scale.sh` (deletion of every 4th work item after each build).
Every table, after insertion and after deletion, equal to the single-threaded
reference with an exact index. Speed-up over 1 thread:

| data | insert 2 / 4 / 8 threads | delete 2 / 4 / 8 threads |
|---|---|---|
| E. coli | 2.00 / 3.86 / **5.44** | 1.95 / 3.65 / **4.31** |
| SRR5833294 (2 M reads) | 1.79 / 4.07 / **5.19** | 1.50 / 3.65 / **4.07** |
| chr21 | 1.44 / 2.86 / 1.41 | 1.95 / 3.31 / 0.99 |

* **chr21 rerun alone, laptop plugged in**: insertion 1.89 / 3.42 / **4.53×**
  (8 threads: 1.70 s, 23.6 M k-mers/s, 4.43× `KmerTable`), deletion 1.91 / 3.34
  / **4.40×** (8 threads: 0.68 s, 4.70× `KmerTable`); lock waits 0.37 s of 13.6
  thread-seconds (insertion) and 0.005 s (deletion); all exact. The 8-thread
  drop above came from the machine, as the analysis below suggested.
* The machine was much slower during this run than during the previous one, in
  code that has not changed since: single-threaded `KmerTable` insertion took
  11.4 s on chr21 (7.25 s before) and 37.8 s on SRR (24.9 s), one-thread
  lock-free queries ran at half the rate, `into_table` took twice as long.
  Absolute rates are not comparable with the previous run.
* chr21 at 8 threads: neither insertion nor deletion was held up by the table's
  locks (lock waits 0.6 s of 50 thread-seconds for insertion, 0.05 s for
  deletion), yet deletion took as long as with 1 thread; insertion saw 1993
  involuntary context switches (55 at 4 threads): the threads were preempted.
  `unihash-scale` now also reports, per row, the share of thread-time spent on a
  CPU, a single-thread CPU speed check relative to the start of the run, piece
  lock waits, waits on open pieces and stale retries of deletions, to separate
  table contention from a busy or throttled machine.

## Concurrent deletion (2026-10-28, sandbox)

`ConcurrentKmerTable::delete_read` from several threads (`unihash-scale
--capacity --check --delete-frac 0.25`: after each build, every 4th work item
is deleted). Rates in M k-mers looked up per second; every result checked
against `KmerTable::delete_read` on one thread (same k-mers, exact index).

| data | `KmerTable` | 1 thread | 2 threads | 2 vs 1 | 2 vs `KmerTable` |
|---|---|---|---|---|---|
| chr21 (179 windows of 64 kbp, 9.1 M k-mers removed) | 3.14 | 3.61 | 5.43 | 1.50× | 1.73× |
| reads (synthetic, 100 k reads, 4.3 M removed) | 3.34 | 3.47 | 5.46 | 1.58× | 1.63× |

* One concurrent deleter is as fast as `KmerTable` or faster: the right part of a
  cut is always a new piece, so items are added and removed rather than re-keyed
  in place, which turns out no slower.
* Lock waits are negligible (0.03 % of acquisitions contended).
* Tested in `tests/m11_concurrent_delete.rs` with deleters cutting overlapping
  stretches, writers inserting new reads and a reader checking that a k-mer never
  deleted is always found and that nothing never inserted is ever found. Breaking
  the cut order (shrinking the piece before publishing the right part) or
  cutting open runs makes these tests fail. Under load (two test processes on
  2 cores) a debug assertion caught a second cut reaching a right part before
  all its index items were filed; right parts are now flagged like open runs
  until filed (80 stressed runs clean). Skipping the gate that keeps
  promotions out of cuts did not make them fail on 2 cores: that window is
  narrow, and the gate is kept by construction.

## Concurrent insertion on the laptop with `with_capacity` (2026-10-28)

`CAPACITY=1 tools/scale.sh` on the user's PC (4 cores / 8 SMT threads): the
concurrent tables are sized with `with_capacity` (E. coli 1024 shards, chr21
8192, SRR5833294 16384). All tables exact. Insertion speed-up over 1 writer
(and over `KmerTable`):

| data | 2 threads | 4 threads | 8 threads | 8 threads before (1024 shards, growing) |
|---|---|---|---|---|
| E. coli | 1.92× | 3.47× | 3.90× (8.8× vs `KmerTable`¹) | 4.82× |
| chr21 | 1.82× | 1.25ײ | **3.38×** (3.27×) | 1.71× |
| SRR5833294, 2 M reads | 2.04× | 3.84× | **5.05×** (5.65×) | 2.22× |

* SRR5833294: 22.2 s with 1 writer (12 % faster than `KmerTable`, which grows
  its index) → 4.4 s with 8 (20.9 M k-mers/s); near-linear up to the 4 physical
  cores, +31 % from SMT. Lock waits: 0.13 s of 35 thread-seconds at 8 threads
  (0.04 % of acquisitions contended), against seconds before. Every input
  fraction now scales 5–6× at 8 threads (1.9× before from half of the reads on).
* chr21 at 8 threads: 3.38× (2.2 s) instead of 1.71×.
* ² chr21 at 4 threads (6.0 s, slower than 2 threads) is a measurement outlier,
  not contention: lock waits 0.04 s, 80 involuntary context switches. The same
  run shows other signs of outside load: `into_table` (single-threaded,
  independent of the thread count) took 5.9 s after the 8-thread build instead
  of 2 s, and ¹ the single-threaded E. coli reference took 1.2 s instead of
  0.5 s (hence the inflated "vs `KmerTable`" ratios on E. coli).
* Memory while inserting: SRR 3.23 B/k-mer, below `KmerTable`'s peak (4.03);
  chr21 2.41 (`KmerTable` 2.02, peak 2.58). `into_table` gives the
  `KmerTable` size (3.05 / 2.01).
* Next bottleneck: `into_table` is single-threaded (5.7–7 s on SRR, more than
  the 8-thread insertion itself); on reads `with_capacity` under-estimates the
  index slightly (about 0.15 entries per k-mer instead of the 0.125 assumed:
  20 k shard growths, 2 thread-seconds).

## Concurrent insertion on the laptop: SMT and lock convoys (2026-10-28)

User's PC: 4 physical cores, 8 hardware threads (SMT). `tools/scale.sh`, 1024
shards, every table checked against the single-threaded one (all exact).
Insertion speed-up over 1 writer:

| data | 2 threads | 4 threads | 8 threads | lock-free queries 4 → 8 threads |
|---|---|---|---|---|
| E. coli (8 MB table) | 1.98× | 3.58× | 4.82× | +40 % |
| chr21 (70 MB) | 1.98× | 2.38× | **1.71×** | +29 % |
| SRR5833294, 2 M reads (245 MB) | 1.91× | 3.33× | **2.22×** | +34 % |

**SMT explains E. coli only.** A second hardware thread on a core typically
adds 20–40 %; E. coli insertion gains 35 % from 4 to 8 threads, like the
lock-free queries on every data set. SMT cannot make 8 threads *slower* than 4,
nor explain chr21 at 4 threads (2.38× while its queries scale 4.1×). On SRR, 8
threads run at 4.3–5.0× on the first quarter of the reads and drop to 1.9× from
half of them on: the loss grows with the table.

**Diagnosis (sandbox, instrumented: lock waits, shard growths, context
switches).** Writers hold a shard's write lock for each lookup-and-append and
across consecutive k-mers of the same shard; with 1024 shards and about 660 k
lock acquisitions per second per writer, every writer needs every shard about
once per millisecond. So whenever a holder is slow, the others pile up on its
shard and sleep:

* Shard growth runs under the lock: 16 k growths take 15 % of chr21's insertion
  time, about 75 µs each on chr21 and growing with the shard. Growth is what
  makes the loss depend on table size.
* A preempted holder is off the CPU for a scheduler time slice (milliseconds).
  With as many writers as hardware threads (8 on 8, plus desktop processes)
  writers get preempted, and every other writer stalls on the preempted
  holder's shard within a millisecond. Reproduced on the 2-core sandbox with 4
  writers (oversubscribed): involuntary context switches 140 → 2952, lock waits
  0.24 → 5.8 thread-seconds, and 4 writers slower than 2 (5.30 vs 5.09 s).
* Plain collisions are rare (0.1 % of acquisitions) and cheap once spinning.

| chr21, sandbox | 1 writer | 2 writers | 4 writers (2 cores) | lock waits at 2 / 4 |
|---|---|---|---|---|
| 1024 shards, growing | 7.96 s | 5.09 s | 5.30 s | 0.24 / 5.79 s |
| `with_capacity` (8192 shards, no growth) | **6.77 s** | **3.97 s** | **3.82 s** | 0.04 / 2.29 s |

Changes: `ConcurrentKmerTable::with_capacity(cfg, kmers)` (shards scale with
the table and start large enough not to grow; one writer becomes as fast as
`KmerTable`), shard locks on separate cache lines, spinning before sleeping on
a contended lock, and contention counters (`ConcurrentStats::contention`,
printed by `unihash-scale` with context switches and page faults). The
remaining limit is that lookups hold the lock; lock-free lookups (bucket
version counters) would leave only appends under it.

## Concurrent insertion: thread and input-size scaling (2026-10-28, sandbox)

`ConcurrentKmerTable` (1024 layer-1 shards) built by `unihash-scale --check`
with 1 and 2 writer threads on the 2-core sandbox (the laptop run with up to 8
threads is `tools/scale.sh`). Every multi-threaded table was checked against
the single-threaded `KmerTable`: same k-mers, none stored twice, exact index.
Rates in M k-mers/s (k-mers without N); memory per distinct k-mer.

| data | `KmerTable` | 1 writer | 2 writers | 2 vs 1 | 2 vs `KmerTable` | concurrent table / peak | after `into_table` (`KmerTable`) |
|---|---|---|---|---|---|---|---|
| E. coli | 10.64 | 8.48 | 15.33 | 1.81× | 1.44× | 2.02 / 2.02 B | 1.82 B (1.82 B, peak 2.45) |
| chr21 | 6.42 | 5.24 | 8.91 | 1.70× | 1.39× | 2.33 / 2.33 B | 2.01 B (2.02 B, peak 2.58) |
| reads (synthetic, 1.6 M × 100–150 b) | 6.32 | 5.15 | 9.51 | 1.85× | 1.50× | 3.05 / 3.05 B | 3.44 B (3.45 B, peak 4.66) |

Queries with 1 / 2 threads (M/s): `KmerTable` (no lock) and the concurrent table
(shard read locks), uniform positives and random negatives:

| data | `KmerTable` positive | random | streaming | concurrent positive | random |
|---|---|---|---|---|---|
| E. coli | 3.84 / 7.72 | 5.01 / 10.01 | 45.9 / 90.4 | 3.40 / 6.43 | 4.47 / 7.61 |
| chr21 | 1.95 / 4.25 | 2.55 / 6.42 | 39.0 / 45.8 | 1.74 / 2.81 | 2.62 / 4.86 |
| reads | 1.85 / 3.62 | 3.25 / 6.57 | 16.7 / 27.0 | 1.55 / 3.25 | 3.12 / 5.95 |

* **Exactness.** Also tested in `tests/m10_concurrent.rs`: up to 8 writers
  inserting the *same* reads in different orders and chunkings (1 to 64
  shards, hot thresholds down to 1, pieces of 9 bases), readers during
  insertion (a completed read is always found, nothing outside the final set
  is ever reported), and round trips with deletion and both defrags. A
  deliberately injected race (lock released between lookup and append) makes
  these tests fail immediately.
* **One writer** is 15–20 % slower than `KmerTable`: a lock per run of k-mers
  sharing a shard, atomic piece records, the segmented base array (candidate
  checks in the writer's own piece and chunk bypass it) and index growth that
  reads pieces through it.
* **Memory.** The concurrent table carries the fixed costs and slack of 1152
  shards (+11–15 % on E. coli and chr21), but its peak equals its size: shards
  grow one at a time, so there is no whole-index growth spike (the
  `KmerTable` peak is 1.3–1.35× its size). On reads it is even 12 % smaller
  than `KmerTable`: 1024 shards at independent points of their ×1.5 growth
  average out the load factor, while the single index happened to have just
  grown. `into_table` (1 s on chr21) gives the same size as `KmerTable`.
* **Input size.** On prefixes of chr21 (¼, ½, all) `KmerTable` inserts 6.4–7.1
  M k-mers/s and on prefixes of the reads (⅛ … all) 6.4–6.9: the cost per k-mer
  hardly grows with the table. Bytes per k-mer follow the index's growth steps
  (2.02–2.75 B on chr21, 2.72–3.45 B on reads).
* The sandbox's two vCPUs scale lock-free queries by only 1.5–2.2×, so the
  insertion efficiencies above (85–92 %) are close to what the hardware gives.

## Eager run indexing (2026-10-28, sandbox)

`run_indexing: Eager` files each new k-mer in the index as soon as it is
appended to its piece, instead of indexing the whole run when it ends. Same
binary, both modes, alternated, 3 runs each; *baseline* is the build before the
change (deferred only). Insertion rate is the median in M k-mers/s (2-core
sandbox, ±10 % run to run). Memory per distinct k-mer, after `shrink_to_fit`.

| data | baseline | deferred | **eager** | redundant k-mers (deferred → eager) | table (deferred → eager) | pieces (deferred → eager) |
|---|---|---|---|---|---|---|
| E. coli | 9.57 | 9.91 | **11.29** | 2 760 → 0 | 1.819 → 1.819 B/k-mer | 1 643 → 1 724 |
| chr21 | 6.09 | 6.32 | **6.75** | 81 008 → 0 | 2.021 → 2.020 B/k-mer | 119 490 → 122 335 |
| reads (synthetic repeats, 400 k × 100–150 b, 0.5 % errors) | 7.93 | 7.78 | 7.78 | 14 215 → 0 | 2.913 → 2.910 B/k-mer | 258 822 → 259 088 |

Laptop (`tools/eager_compare.sh`, medians of 3 runs, deferred → eager):
E. coli 9.82 → 10.33, chr21 5.81 → 5.89, SRR5833294 (2 M reads) 3.71 → 3.76
M k-mers/s; identical memory (3.056 B/k-mer on SRR5833294); k-mers stored twice
2 760 / 81 008 / 5 810 → 0. Eager run indexing is now the default.

* Same k-mer sets (34 625 350 distinct on chr21 in both modes). Eager stores
  each exactly once; deferred stored k-mers repeated inside one run twice
  (0.23 % of chr21). Eager has slightly more pieces because a run now ends at
  an intra-run repeat instead of storing it again; the table is the same size.
* Insertion is as fast or faster: there is no pending ASCII buffer to re-encode
  and no second pass computing the piece's ℓ-mer values and minimizer
  positions. The gain is largest on long records (E. coli, chr21), where runs
  are long.
* Peak heap during insertion and query rates are unchanged (within noise).
* The deferred mode of the new build matches the baseline: the change does not
  slow the existing path.
* Correctness: `tests/m9_eager.rs` (both modes against the oracle and each
  other, identical pieces and index entries on repeat-free input, chunking
  independence, open-read visibility, deletion and both defrags) and eager
  variants in `tests/m8_delete.rs`, the bench and the churn smoke runs.

## Incremental defragmentation on the laptop (2026-10-25)

User's PC, `CHECK=1 tools/churn_real.sh`. Every policy on every data set ended
with **exact set equality** with the HashSet (E. coli 4 299 170, chr21
33 132 009, SRR5833294 78 022 310 k-mers in both, 0 only in either) and an
exact index, also after the final full defrag.

| data | policy | full + inc | defrag share | mean table vs start | end table |
|---|---|---|---|---|---|
| E. coli | never | 0 + 0 | 0 % | +18 % | 11.38 MB |
| | `overhead:0.1` | 5 + 0 | 16.0 % | +1.0 % | 8.00 MB |
| | **`inc:0.02+full:0.1`** | 0 + 15 | **1.0 %** | +4.0 % | 8.71 MB |
| chr21 | never | 0 + 0 | 0 % | +5.2 % | 76.28 MB |
| | `overhead:0.1` | 1 + 0 | 11.9 % | +0.9 % | 70.46 MB |
| | **`inc:0.02+full:0.1`** | 0 + 5 | **0.6 %** | **+1.2 %** | 69.71 MB |
| SRR5833294 (500 k reads) | never / `overhead:0.1` | 0 + 0 | 0 % | +2.4 % | 211.56 MB |
| | **`inc:0.02+full:0.1`** | 0 + 2 | **3.0 %** | **+1.2 %** | 203.13 MB |

* chr21: the same memory as full defrags for 1/20 of the defrag time
  (0.09 s vs 2.18 s).
* SRR5833294: the full-only policy never triggered (estimated overhead 7.2 %,
  8.1 % reclaimed by a final `defrag`). Two incremental calls (0.2 s each)
  removed the dead bases: the end table is 203.1 MB instead of 211.6 MB.
* Queries were unchanged under every policy.

## Incremental defragmentation (2026-10-25, sandbox)

`defrag_incremental()` reclaims what deletions free without rebuilding the
index. It slides the base array over dead bases (pieces keep their ids, so only
their starts change; 32-base word moves), drops unlinked overflow blocks, and
re-sizes an index layer once its slack is 1/16 of it. Policy `inc:X+full:Y`
calls it when `incremental_overhead() > X` (dead bases + re-sizable slack) and
calls a full `defrag()` when `piece_overhead(k) > Y` (extra pieces, shared
bases, emptied records). Raw outputs: `results/churn/*_incremental.txt`. Every
run's answers and index were checked in the tests and in the CI churn run.

Mean table relative to the starting (defragmented) size. Defrag share is the
fraction of edit + defrag time.

| | policy | full + incremental calls | defrag time | share | mean table |
|---|---|---|---|---|---|
| E. coli (40 k × 300 b) | never | 0 + 0 | 0 | 0 % | +18 % |
| | `overhead:0.1` (full only) | 5 + 0 | 0.80 s | 15.4 % | +1.0 % |
| | `inc:0.02+full:0.05` | 3 + 13 | 0.51 s | 10.3 % | −0.1 % |
| | **`inc:0.02+full:0.1`** | 0 + 15 | **0.04 s** | **1.2 %** | +4.0 % |
| chr21 (30 k × 1000 b) | never | 0 + 0 | 0 | 0 % | +5.2 % |
| | `overhead:0.1` (full only) | 1 + 0 | 1.86 s | 10.5 % | +0.9 % |
| | **`inc:0.02+full:0.1`** | 0 + 5 | **0.10 s** | **0.7 %** | **+1.2 %** |
| | `inc:0.05+full:0.1` | 0 + 1 | 0.02 s | 0.1 % | +2.7 % |
| Reads (synthetic, 200 k × read) | never | 0 + 0 | 0 | 0 % | +12 % |
| | `overhead:0.1` (full only) | 6 + 0 | 2.76 s | 20.9 % | −9.0 % |
| | `inc:0.02+full:0.05` | 3 + 13 | 1.83 s | 14.8 % | −8.8 % |
| | **`inc:0.02+full:0.1`** | 0 + 16 | **0.61 s** | **5.8 %** | −4.5 % |

* An incremental call costs 3–40 ms here (E. coli 3 ms, chr21 20 ms, reads
  38 ms when it re-sizes the index), against 0.17 s, 2 s and 0.5 s for a full
  `defrag`.
* **Dead bases are most of the reclaimable memory** (chr21 never-defragmented:
  6.4 MB of 10.1 MB reclaimed); they need no index work. Extra pieces from cuts
  and re-insertions grow slowly: `piece_overhead` never reached 0.1 in these
  runs, so `inc:0.02+full:0.1` made no full `defrag` at all.
* Full-only policies reach lower memory on reads (−9 %), because merging
  re-inserted reads with their neighbours shrinks the table below its start.
  `inc:0.02+full:0.05` gets the same memory for 30 % less defrag time.
* **Merging incrementally does not pay.** An earlier version also joined the
  pieces changed since the last defrag with their neighbours (8 index lookups
  per changed piece, then rewriting the chain and re-filing its index items).
  On chr21 that cost ≈ 230 ms per call for 8.5 k joins (≈ 0.24 MB): ≈ 30 µs per
  28 B saved, about 5× more per reclaimed byte than a full `defrag`. It stays
  available as `defrag_incremental_with(max_mean_len)` and is off by default.
* Queries are unaffected by any policy (random 2.1–2.3 M/s, streaming
  19–21 M/s on chr21).

**Recommended policy:**

1. Call `defrag()` once after bulk-loading reads.
2. Then call `defrag_incremental()` whenever `incremental_overhead() > 0.02`.
3. Call `defrag()` whenever `piece_overhead(k) > 0.10`.

## Deletion and defragmentation on the laptop (2026-10-25)

User's PC, `tools/churn_real.sh` and `unihash-bench --defrag` (outputs pasted in
the session; the churn runs used the overhead estimate of 2026-10-24, before the
index-slack term below was added).

### Defragmenting after loading SRR5833294 (2 M reads)

| | before | after `defrag` |
|---|---|---|
| table | 244.8 MB (3.05 B/k-mer) | **201.4 MB (2.51 B/k-mer), −18 %** |
| pieces | 2 069 193 | 1 756 874 (312 319 joins) |
| bases | 142.3 M | 132.9 M |
| vs kache-hash (20.29 B/k-mer) | 6.6× less | **8.1× less** |

`defrag` took 7.1 s (insertion 23.6 s; ≈ 29 ms per MB). Every query pattern
still had 0 mismatches and 0 false positives. Query speeds were close to the
earlier run; the HashSet baseline was also slower in this run (random 3.2 vs
5.3 M/s), so the machine was busier. Its peak heap was 369 MB, above the
insertion peak (336 MB), because the end-pairing temporaries (≈ 124 MB) were
built while the indexes were still allocated. **Fixed**: the indexes are now
freed first and the new heads allocated only after the old base array is
dropped. On the synthetic reads, the defrag peak fell from 48.0 MB to 30.9 MB,
i.e. the size of the table before defrag.

### Churn on the laptop

| data | policy | defrags | defrag share | mean table vs start | queries |
|---|---|---|---|---|---|
| E. coli (40 k × 300 b) | never | 0 | 0 % | +19 % (end +39 %) | unchanged |
| | **0.10** | 4 | **14 %** | **+1.3 %** | |
| | 0.05 | 10 | 26 % | −1.9 % | |
| | 0.20 | 1 | 4.7 % | +8.5 % | |
| chr21 (30 k × 1000 b) | never | 0 | 0 % | +5.7 % (end +12 %) | unchanged |
| | **0.10** | 1 | **12 %** | **+3.3 %** | |
| | 0.05 | 2 | 21 % | −0.3 % | |
| SRR5833294 (500 k reads) | never / 0.10–0.40 | 0 | 0 % | +3.1 % (end +6.7 %) | unchanged |
| | 0.05 | 1 | 36 % | +0.7 % | |

Edit time: 3.6 s (E. coli), 15.2 s (chr21), 11.7 s (500 k reads deleted and
re-inserted, 44 M k-mer operations). The read set fragments slowly: a quarter
of the reads deleted and re-inserted raised the reclaimable fraction to 9.5 %
(214.9 → 194.4 MB at the final defrag, 6.6 s).

### Overhead estimate recalibrated

The 2026-10-24 estimate ignored index slack: slots freed by deletions and partly
empty overflow blocks, which `defrag` returns by re-sizing the index. On reads
that slack is ~ 3 % of the table, so the estimate ran low (7.1 % estimated vs
9.5 % reclaimed on SRR5833294). Each index now records its bytes per entry after
its last rebuild (growth or `defrag`), and `Fragmentation::index_slack_bytes` is
the excess over that. Estimate vs reclaimed by one final `defrag` (`never`
policy, sandbox):

| | before | now | reclaimed |
|---|---|---|---|
| E. coli | 31.3 % | 32.5 % | 29.1 % |
| chr21 | 12.0 % | 15.8 % | 14.1 % |
| synthetic reads | 29.4 % | 35.4 % | 34.0 % |

The estimate now lands within 1.5–3.5 points above the actual reclaim on all
three. `unihash-churn` prints both values.

### Decision

* `defrag()` once after bulk-loading reads (−18 % on SRR5833294, −20 % on the
  synthetic reads).
* During deletions and re-insertions, `defrag()` when
  `fragmentation().overhead(k) > 0.10`: ≈ 12–15 % of the edit time, mean memory
  within +1–3 % of the compact layout on every data set. 0.05 or less doubles or
  triples that time for ≈ 2 % less memory.
* Fragmentation never slowed queries measurably.

## Deletion and defragmentation (2026-10-24, sandbox)

Same 2-core sandbox as the next section (noisy: ±20 % between identical runs);
compare within a table. Raw outputs: `results/churn/`. All runs with `--check`
in CI and the deletion tests report 0 mismatches against the oracle and an
exact index (`verify_index`).

### Deletion throughput

Deleting random units from a freshly built and defragmented table, no
re-insertion (`unihash-churn --lag ∞ --policy never`):

| | units | input k-mers | stored k-mers removed | time | input k-mers/s |
|---|---|---|---|---|---|
| chr21 | 20 000 × 1000 bases | 19.4 M | 12.3 M | 6.23 s | 3.1 M |
| Reads (synthetic, 400 k reads) | 100 000 reads | 12.0 M | 5.4 M | 3.48 s | 3.4 M |

About half the insertion rate (5.9 and 7.8 M/s). A cut re-files the index items
of the right part of the piece, so its cost grows with the piece length; this
is why `max_piece_len` now defaults to 4096 (memory cost on chr21: +0.02 %).
Memory is not returned by deletion itself: dead bases and freed index slots
stay allocated until `defrag`.

### Defragmenting a freshly built table

| | table before | after | pieces | `defrag` time |
|---|---|---|---|---|
| E. coli | 8.28 MB | 8.32 MB | 1643 → 1629 | 0.2 s |
| chr21 | 69.96 MB | 68.86 MB (−1.6 %) | 119 490 → 93 105 | 2.0 s |
| Reads (synthetic, 400 k) | 30.92 MB | **24.79 MB (−20 %)** | 284 303 → 197 224 | 0.55 s |

Overlapping reads leave pieces that extend each other by exactly k − 1 bases.
Merging them saves the piece records, the overlaps and the junction index
entries: **3.01 → 2.42 B/k-mer on reads**, with unchanged query speed
(`results/churn/reads_bench_{defrag,nodefrag}.txt`). The defrag peak heap (48 MB)
stays below the insertion peak (50 MB). On E. coli the index is re-sized to a
lower load, so the table grows by 0.5 %.

### Churn: when to call `defrag`

Each operation deletes one unit and re-inserts the unit deleted `lag`
operations earlier, so the k-mer set stays roughly constant while pieces
fragment. Policy `overhead:X` calls `defrag` when
`fragmentation().overhead(k) > X`. "Mean table" is sampled 200 times over the
run. Throughputs are for the final state.

E. coli, 40 000 × 300-base units (2.6 × the genome churned), lag 1000, start
8.32 MB:

| policy | defrags | edit s | defrag s | defrag share | mean table | max table | random | streaming |
|---|---|---|---|---|---|---|---|---|
| never | 0 | 4.61 | 0 | 0 % | 9.91 MB (+19 %) | 11.53 MB | 3.98 M/s | 36.5 M/s |
| 0.02 | 26 | 6.13 | 4.39 | 42 % | 7.99 MB | 8.38 MB | 3.89 | 28.3 |
| 0.05 | 10 | 5.39 | 1.74 | 24 % | 8.16 MB | 8.60 MB | 4.90 | 38.3 |
| **0.10** | 4 | 4.48 | 0.64 | **13 %** | **8.43 MB (+1 %)** | 8.98 MB | 4.06 | 35.5 |
| 0.20 | 1 | 4.00 | 0.17 | 4 % | 9.03 MB (+9 %) | 9.97 MB | 3.68 | 29.4 |

chr21, 30 000 × 1000-base units (0.64 × the genome), lag 2000, start 68.86 MB:

| policy | defrags | edit s | defrag s | defrag share | mean table | max table | random | streaming |
|---|---|---|---|---|---|---|---|---|
| never | 0 | 15.59 | 0 | 0 % | 72.80 MB (+6 %) | 77.12 MB | 2.06 M/s | 18.8 M/s |
| 0.05 | 2 | 15.94 | 3.63 | 19 % | 68.62 MB | 71.68 MB | 2.01 | 19.5 |
| **0.10** | 1 | 16.14 | 2.03 | **11 %** | 71.10 MB | 75.47 MB | 2.38 | 20.8 |
| 0.20 | 0 | 14.66 | 0 | 0 % | 72.80 MB | 77.12 MB | 1.85 | 17.5 |

Reads (synthetic, 400 k reads), 200 000 reads deleted and re-inserted, lag
20 000, start 24.79 MB:

| policy | defrags | edit s | defrag s | defrag share | mean table | max table | random | streaming |
|---|---|---|---|---|---|---|---|---|
| never | 0 | 9.88 | 0 | 0 % | 28.53 MB (+15 %) | 32.79 MB | 3.12 M/s | 20.8 M/s |
| 0.05 | 10 | 10.46 | 4.52 | 30 % | 22.26 MB | 25.34 MB | 3.33 | 21.9 |
| **0.10** | 4 | 10.54 | 1.81 | **15 %** | 23.38 MB | 26.41 MB | 3.21 | 22.1 |
| 0.20 | 1 | 10.33 | 0.47 | 4 % | 25.67 MB | 29.22 MB | 3.23 | 20.8 |

Findings:

* **Fragmentation does not slow queries** measurably (random and streaming
  columns, within noise). It costs memory: dead bases (2 bits each) and, per
  cut or re-inserted unit, an extra piece record, k − 1 overlap bases and
  junction entries. Without `defrag`, the table grew by a fraction of its size per
  genome's worth of churned bases (E. coli +39 % after 2.6 genomes, chr21 +12 %
  after 0.64: 0.15–0.2 per genome).
* **`defrag` costs ≈ 25 ms per MB of table** (0.2 s for E. coli, 2 s for chr21,
  0.5 s for 25 MB of reads), i.e. as much as churning ≈ 10 % of the table.
* **`overhead(k)` tracks the reclaimable fraction**: 31.3 % estimated vs 29 %
  reclaimed on E. coli; 12.0 % vs 14 % on chr21 (`never` rows, one defrag at the
  end).
* **Recommended policy: call `defrag` when `fragmentation().overhead(k) > 0.10`**.
  Memory then stays within ≈ 10 % of the compact layout (mean +1 to +3 %; −6 % on
  reads, whose re-inserted reads merge better than the original ones), for
  11–15 % of the edit time. Thresholds of 0.05 or less double or triple that
  time for ≈ 2 % less memory. Use 0.2 when edit throughput matters more: 4 % of
  the time, mean memory +4 to +9 %. And defragment once after bulk-loading
  reads.
* On E. coli, defragmenting often also slows the edits themselves (4.0 → 6.1 s
  from 0.2 to 0.02): merged pieces are long, and a cut re-files the right part.
  With `--max-piece-len 1024` the edit time is flat (3.3–3.9 s) at +1 % memory
  (`results/churn/ecoli_maxpiece1024.txt`).

## Hot-minimizer layers, spilled growth, kache-hash comparison (2026-10-24, sandbox)

2-core cloud sandbox, so absolute rates are lower than on the laptop; compare
within a table. Raw outputs: `results/kache/` and `results/sandbox-2026-10-24b/`.
"Reads" is a synthetic set of 800 k reads (2 × 400 k, `--holdout`) with
errors from a repeat-rich genome, plus a low-overlap query set
(`reads_lo_1`). SRR5833294 was not available in the sandbox.

### unihash vs kache-hash vs `HashSet<u64>`

Memory in heap bytes per distinct k-mer. All other rows in M k-mers/s.
kache-hash is single-threaded, k = 31, ℓ = 17, initial capacity 2^20, and every
answer is checked against the oracle (0 mismatches, 0 false positives).

| E. coli | unihash | kache-hash | HashSet |
|---|---|---|---|
| Memory, B/k-mer | **1.82** | 22.34 | 16.58 |
| Insert | 8.84 | 9.42 | 11.02 |
| k-mers in order | 6.02 | 3.35 | 11.34 |
| k-mers random order | 4.99 | 1.95 | 10.86 |
| sequence, streaming | **46.10** | 19.99 | 10.18 ¹ |
| random negatives | 6.18 | 2.32 | 24.09 |

| chr21 | unihash | kache-hash | HashSet |
|---|---|---|---|
| Memory, B/k-mer | **2.01** | 23.50 | 17.44 |
| Insert | 6.20 | 5.26 | 8.22 |
| k-mers in order | 4.95 | 2.90 | 7.22 |
| k-mers random order | 2.30 | 1.71 | 7.22 |
| sequence, streaming | **31.49** | 15.32 | 5.68 ¹ |
| random negatives | 3.39 | 2.04 | 10.95 |

| Reads (synthetic) | unihash | kache-hash | HashSet |
|---|---|---|---|
| Memory, B/k-mer | **3.01** | 19.82 | 14.71 |
| Insert | 9.43 | 11.34 | 11.37 |
| inserted: in order / random / streaming | 5.39 / 2.89 / **25.41** | 2.92 / 1.67 / 14.46 | 8.74 / 7.33 / 7.27 ¹ |
| held-out: in order / random / streaming | 5.22 / 3.11 / **24.48** | 3.10 / 1.73 / 16.03 | 9.52 / 9.16 / 7.46 ¹ |
| low-overlap set: in order / random / streaming | 5.64 / 3.87 / 19.63 | 2.94 / 2.01 / **21.29** | 15.09 / 14.80 / 15.26 ¹ |
| random negatives | 4.52 | 2.21 | 17.22 |

¹ HashSet "streaming" is independent lookups of the k-mers parsed from the same
sequence.

* **Memory.** unihash uses **7–12× less memory than kache-hash on genomes and 6.6× on reads** and 5–9× less
  than a `HashSet<u64>`. kache-hash's ~20–24 B/k-mer includes its power-of-two
  capacity: 64 M slots for 34.6 M k-mers on chr21. Its peak RSS on chr21 was
  1280 MiB, against 83–102 MiB peak RSS for unihash.
* **Streaming queries**, kache-hash's native mode, are **1.5–2.3× faster** in
  unihash on genomes and on reads with hits. kache-hash wins (21.3 vs 19.6) only on
  the low-overlap query set, where most k-mers are absent and unihash's
  hit-extension shortcut rarely applies.
* **Pre-parsed k-mers.** kache-hash only takes k-mers through a `Kmer_Window`
  built from characters. Decoding to ASCII is done outside the timer, but `init`
  (hash and minimizer) is timed, like unihash recomputing the minimizer of a
  packed k-mer. unihash is 1.3–2.6× faster here.
* **Insertion.** kache-hash inserts reads faster (11.3 vs 9.4) and chr21 slower
  (5.3 vs 6.2).
* **Build caveat.** kache-hash (commit `f7cfa5b`) must be compiled with
  `-fno-strict-aliasing`. At `-O3` without it, GCC 13 produced wrong answers
  (streaming found 36 k of 4.5 M present k-mers). `tools/kache_compare.sh` sets
  the flag and fails on any mismatch.

### Hot-minimizer routing (layers 2 and 3)

chr21, `--no-oracle`, sweep of `--hot-threshold` (0 = routing off; default 64):

| T | insert | heap B/k-mer | hot ℓ / ℓ2 | L2 entries | L3 k-mers | in order | random | streaming | negatives |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 4.49 | 1.924 | 0 / 0 | 0 | 0 | 2.80 | 1.51 | 18.20 | 3.49 |
| 16 | **6.82** | 2.157 | 4623 / 2112 | 423 k | 153 k | 4.69 | **2.32** | 28.64 | 3.33 |
| **64** | 6.60 | 2.007 | 692 / 98 | 228 k | 22.8 k | **4.82** | 2.30 | **33.05** | 3.43 |
| 256 | 6.20 | 1.959 | 92 / 1 | 97 k | 709 | 4.57 | 2.30 | 29.88 | 3.70 |

* With T = 64, chr21 gains **+47 % insertion, +72 % in-order, +52 % random and
  +82 % streaming** for **+4 % memory** (0.08 B/k-mer: layer 2 0.105, layer 3
  0.009, minus the entries drained from layer 1).
* Only 692 ℓ-minimizers (of millions) are hot, and 98 of them are hot again at
  ℓ2 = 25. The layer-3 `HashSet` holds 22.8 k k-mers (0.07 %).
* E. coli has a single hot minimizer, and the synthetic reads have none, so
  routing changes nothing there (insert 11.42 vs 12.38 M/s on E. coli, which is
  within noise; identical table bytes on reads).

### Growth without coexisting indexes

Peak heap during insertion (counting allocator) with the in-memory spill
(default) and with `--spill-dir`. Before this round, the old and new block
arrays coexisted, and chr21 peaked at 137 MB.

| | final table | peak, memory spill | peak, disk spill | peak RSS, memory → disk |
|---|---|---|---|---|
| E. coli | 8.27 MB | 28.5 MB | 26.5 MB | — |
| chr21 | 69.7 MB | 98.3 MB (1.41×) | **85.3 MB (1.22×)** | 102.1 → 83.2 MiB |
| Reads | 30.9 MB | 49.7 MB (1.61×) | **39.9 MB (1.29×)** | 49.0 → 39.8 MiB |

The remaining overhead is the insertion buffers, `Vec` growth slack of the piece
store, and the overflow area copied by compaction. On E. coli the
fixed-size insertion buffers dominate. The disk spill writes 8 B per entry per
rebuild (73 MB in total over 14 rebuilds on chr21) and removes the file after
each rebuild.

## Flat buckets + packed pieces, uniform workloads (2026-10-24)

**Changes.**

* **Index.** Flat fixed-capacity two-choice buckets (80 B blocks of 8 entries,
  9 B per entry, no per-bucket `Vec` header), with automatic growth by 1.5×.
* **Pieces.** Packed at 2 bits per base.
* **Benchmark.**
  * Query workloads now cover each source uniformly.
  * In-order and random-order use the same k-mers, measured the same way.
  * Table memory is measured alone, before the oracle exists, and cross-checked
    with a counting allocator.
  * A held-out query set was added (`--holdout`).

**Before/after, same machine** (2-core sandbox, k = 31, ℓ = 17). "Before" is the
previous commit with `bucket_bits` tuned per dataset (16 / 20 / 19); "after" starts
from 2^10 buckets and grows. Raw outputs are in `results/sandbox-2026-10-24/`. The
synthetic reads are 400 k × 150 bp from a 5 Mb random genome, 0.5 % errors, random
strand. 0 mismatches and 0 false positives everywhere.

| | E. coli | chr21 | Synthetic reads |
|---|---|---|---|
| **Table, B/distinct k-mer** | 2.86 → **1.82** (−37 %) | 3.40 → **1.93** (−43 %) | 5.07 → **3.01** (−41 %) |
| Pieces, B/k-mer | 1.01 → 0.25 | 1.13 → 0.30 | 2.05 → 0.68 |
| Index, B/k-mer | 1.85 → 1.56 | 2.27 → 1.63 | 3.02 → 2.33 |
| Insert, k-mers/s | 12.5 M → 11.0 M | 3.00 M → **4.05 M** | 7.42 M → **8.85 M** |
| Random order (pre-parsed), k-mers/s | 3.32 M → **4.37 M** | 0.91 M → **1.40 M** | 1.63 M → **2.48 M** |
| In order, from sequence, k-mers/s | 5.02 M → 6.18 M | 2.71 M → 2.71 M | 4.72 M → 5.10 M |
| In order, streaming, k-mers/s | 43.9 M → 46.7 M | 20.8 M → 16.4 M ¹ | 26.6 M → 24.8 M |
| Random negatives, k-mers/s | 4.86 M → 4.38 M | 2.52 M → **3.27 M** | 2.75 M → **4.14 M** |

¹ Not like-for-like. "Before" streamed the first 20 Mb in one piece; "after"
streams 20,623 separate 1 kb windows, each starting cold. On the same contiguous
workload (`--window 20000000`) the two are equal: 20.1 vs 19.7 M/s.

**Table alone, exact (after).**

| | E. coli | chr21 | Synthetic reads |
|---|---|---|---|
| Heap, accounted = allocator ² | 8.27 MB | 66.8 MB | 30.9 MB |
| Peak heap during insertion | 32.6 MB ³ | 137 MB | 64.7 MB |
| RSS added by the table | 8.2 MiB | 64.9 MiB | 32.8 MiB |
| Oracle `HashSet<u64>` heap | 75.5 MB (16.6 B/k-mer) | 604 MB (17.4) | 151 MB (14.7) |

² The two differ by the same 256 B in every run, so the accounting covers the table exactly.
³ Includes the bounded insertion buffers: 12 MB on E. coli, where one run of absent
k-mers reaches `max_piece_len` = 1 Mb. The peak also covers index growth, during
which the old and new block arrays coexist, and `Vec` slack.

**Same measuring for all patterns: the earlier anomaly is gone.** With the same
k-mers, pre-parsed in both cases, the oracle runs at the same speed in order and
shuffled (chr21: 5.55 vs 5.52 M/s; reads: 6.87 vs 7.55 M/s). The earlier gap was a
measurement artifact. For the table, random order costs about 2× in-order on chr21
and the reads (cache misses on pieces and blocks), and 1.4× on E. coli.

**Held-out reads (records 400,001–800,000), synthetic.** Hit rate 0.891. In-order
lookups run at 5.1–5.2 M/s, random order at 2.79 M/s, and streaming at 24.1 M/s
(86 % answered by extension). The oracle runs at 7.2–7.5 M/s.

**What got worse or did not move.**

* **More entries scanned on human data.** Positive chr21 queries scan more entries
  (20.3 vs 7.5 per query) because a lookup always visits two chains, and hot
  (Alu, poly-A) chains are shared by more minimizers. Tags are compact, so this
  still runs faster in random order, but it caps in-order independent lookups at
  the previous speed. Routing hot minimizers is the next step.
* **E. coli negatives are 10 % slower.** Its previous tuned configuration had one
  bucket per query; now two blocks are read.
* **Insertion peak is about 2× the final table**, because of growth.

## Laptop run, all query patterns (`tools/bench_real.sh`, 2026-10-23)

Raw outputs: `results/{ecoli,chr21,SRR5833294}.txt`. k-mers/s, table / `HashSet<u64>`
oracle, with the table-to-oracle speed ratio in brackets. 0 mismatches and 0 false
positives everywhere.

| | E. coli | chr21 | SRR5833294 (2 M reads) |
|---|---|---|---|
| Memory, B per distinct k-mer | 2.86 / 16.58 | 3.40 / 17.44 | 5.01 / 15.06 |
| Insert | 11.5 M / 8.0 M (1.44×) | 3.18 M / 5.66 M (0.56×) | 2.11 M / 5.15 M (0.41×) |
| In order, independent | 5.39 M / 4.35 M (1.24×) | 3.07 M / 3.65 M (0.84×) | 3.89 M / 2.92 M (1.33×) |
| **In order, streaming** | **46.5 M** / 4.35 M (10.7×) | **21.8 M** / 3.65 M (6.0×) | **17.2 M** / 2.92 M (5.9×) |
| Random order (shuffled sample) | 2.08 M / 8.19 M (0.25×) | 0.89 M / 6.51 M (0.14×) | 0.57 M / 5.77 M (0.10×) |
| Random negatives | 3.23 M / 14.96 M (0.22×) | 2.50 M / 11.97 M (0.21×) | 2.38 M / 11.24 M (0.21×) |
| Streaming answered by extension | 99.96 % | 98.05 % | 97.08 % |

Low-hit query sets SRR5901135, queried in order against the SRR5833294 table.
These are effectively negatives: 6,652 hits in 341 M k-mers for `_1` and 12,240 in
343 M for `_2`.

| | Independent | **Streaming** | Oracle |
|---|---|---|---|
| SRR5901135_1 | 6.25 M | **15.4 M** | 6.65 M |
| SRR5901135_2 | 6.23 M | **15.4 M** | 6.61 M |

Streaming on misses comes entirely from the candidate cache. 87 % of lookups
reuse the previous minimizer's candidates, and entries scanned per k-mer fall
from 2.85 to 0.37.

**Comparison points that are still missing or not like-for-like:**

1. **In order and random order are not timed the same way.** In-order passes
   include parsing k-mers from the sequence; random-order passes time pre-parsed
   `u64`s. On this laptop the oracle is also *slower* in order than in random
   order (for example 2.92 vs 5.77 M/s on SRR5833294). In the sandbox the two were
   equal, and hash-set access is random in both cases. Until in-order is also
   measured on pre-parsed k-mers, compare only within one pattern.
2. **The in-order and random-order samples cover different k-mers.** In-order
   stops after `--pos-limit` bases, so chr21 covers its first 20 Mb (13.4 M
   k-mers). Random order samples the whole input. Entries scanned per query differ
   (5.3 vs 7.5 on chr21).
3. **No realistic high-hit or mixed query set.** The inserted input is queried
   against itself (100 % hits), and SRR5901135 is about 0 % hits. Held-out reads
   of SRR5833294 (records after the first 2 M) would give a high-hit set with
   sequencing errors. A mixed set would exercise extension and the cache together.
4. **No streaming baseline.** The oracle has no in-order fast path. The natural
   comparisons are the reference structures of this project, kache-hash (the
   project paper) and an SSHash-style streaming query. Neither is benchmarked.
5. **No compact flat hash table baseline.** "Ideal u64 table" is a memory
   estimate only (9.14 B/k-mer). There is no speed measurement for an
   open-addressing `u64` table, which would be the fair fast baseline.
6. **No table-only memory on these runs.** Peak RSS includes the oracle (2.18 GB
   on SRR5833294, of which 1.2 GB is the oracle). Table-only RSS was measured once,
   for chr21 in the sandbox.
7. **Partial data and single runs.** Only the first 2 M records of each read set
   were used, there was one run per configuration (no variance), and
   `bucket_bits` was not tuned per dataset (2^22 for reads, 25 % of the table in
   bucket headers).

## Query access patterns: in order, streaming, random order

The bench tool now measures queries in three patterns (2-core sandbox, k-mers/s,
0 mismatches everywhere):

* **in order, independent**: `contains` on each successive k-mer of the input.
* **in order, streaming**: `query_stream`, which reuses the rolling minimizer,
  extends the previous hit by one base, and caches the candidates while the
  minimizer is unchanged.
* **random order**: a uniform sample of the input's k-mers, shuffled.

| Dataset | In order, independent | **In order, streaming** | Oracle, in order | Random order | Oracle, random |
|---|---|---|---|---|---|
| E. coli (4.6 M k-mers, 2^16 buckets) | 4.50 M | **42.1 M** | 6.82 M | 2.48 M | 6.15 M |
| chr21 (first 20 Mb, 13.4 M k-mers; 5 M random) | 2.57 M | **18.0 M** | 5.44 M | 0.85 M | 5.65 M |
| Synthetic reads, 100 bp (1.4 M k-mers) | 4.99 M | **20.6 M** | 32.0 M ¹ | 4.12 M | 35.6 M ¹ |
| Synthetic query set, 50 % hits (0.35 M k-mers) | 5.09 M | **27.7 M** | 29.7 M ¹ | — | — |

¹ These tables are only 2 MB, so the oracle fits in cache.

How each streaming k-mer was answered:

| Dataset | By extension | Index lookup, cached minimizer | Index lookup, bucket scan |
|---|---|---|---|
| E. coli | 99.96 % | 0.01 % | 0.03 % |
| chr21 | 98.05 % | 0.70 % | 1.25 % |
| Synthetic reads (inserted) | 96.2 % | 0.4 % | 3.4 % |
| Synthetic query set, 50 % hits | 47.7 % | 43.5 % | 8.8 % |

* **Streaming is 7× (chr21) to 9× (E. coli) faster than independent lookups**, and 3.3–6×
  faster than the hash set queried in the same order. Almost every positive k-mer
  is answered by one base comparison.
* **On mixed hit/miss data, the candidate cache does the work.** Misses cannot be
  extended, but the minimizer usually stays the same, so 43 % of lookups reuse the
  cached candidates without scanning the bucket.
* **Random order costs about 3× more than in-order independent lookups on chr21**
  (0.85 vs 2.57 M/s). The two structures (pieces 39 MB, index 78 MB) are then hit
  at random, so each query takes cache misses in the bucket array, the entries and
  the piece bytes. The hash set's random-order rate barely changes. This is the
  case a flat bucket layout and packed pieces should improve.

## Update: hashed order only, streaming insertion

Following the chr21 results below, the lexicographic order was removed: ℓ-mers are
always ranked by the hashed order. Insertion is now streaming (`insert_chunk` /
`end_read`), and the bench tool streams FASTA records in chunks instead of loading
whole chromosomes.

chr21, 2^20 buckets, `--no-oracle`, measured in the 2-core sandbox:

| | Before | After |
|---|---|---|
| Peak RSS at the end of insertion | 492.6 MiB | **162.2 MiB** |
| Table heap (pieces + index) | 112.2 MiB | 112.2 MiB (identical pieces and index) |
| Insertion buffers | 16 B per base of the longest segment | 2.6 MB (bounded by about 14 B × `max_piece_len`) |
| Insert throughput | 2.98 M/s | 3.05 M/s |

With the oracle, chr21 peak RSS fell from 1.37 GB to 1.01 GB, which is now mostly
the 604 MB oracle. The remaining overhead above the table is `Vec` growth slack
before `shrink_to_fit` plus I/O buffers. Rows for the hashed order in the tables
below still apply, because the stored pieces and index entries are unchanged.

## Read set SRR5833294 (laptop, first 2 M reads)

`results/SRR5833294.txt`. 2 M reads of 76 bp, 91.9 M k-mers, 80.2 M distinct,
2^22 buckets.

| | Table | Oracle `HashSet<u64>` |
|---|---|---|
| Memory per distinct k-mer | **5.01 B** (pieces 1.98, bucket headers 1.26, entries 1.77) | 15.06 B |
| Insert | 1.60 M/s ¹ | 3.44 M/s ¹ |
| Positive query | 3.28 M/s | 5.42 M/s |
| Negative query | 2.39 M/s | 11.22 M/s |

* Pieces are short (mean 69 bp, 0.007 % redundancy), so the index is 60 % of the
  table, and a quarter of the table is bucket headers (2.8 entries per bucket). With
  2^21 buckets the table would save 0.63 B/k-mer but scan about twice as many
  entries per negative query. A flat bucket layout would avoid having to choose.
* 13 % of k-mers were already present on insertion. This set is high-hit, so
  overlapping reads are mostly skipped.
* ¹ Insert rates from this run, and the E. coli/chr21 insert rates from the same
  run, are understated (see the next point). Re-run `tools/bench_real.sh` after the
  fix.
* **This run aborted before querying SRR5901135**, because the bench tool stripped
  `_` from every option value, so `SRR5901135_1.fastq.gz` became
  `SRR59011351.fastq.gz`. This is fixed. Numbers still accept `1_000_000`, file
  paths are kept verbatim, and all input and query files are checked before
  insertion starts. `bench_real.sh` now queries both SRR5901135 mates.
* **Insert timing artifact, fixed.** After the streaming change, the bench timed
  every chunk, and FASTA chunks were single 50–80 bp lines. Table and oracle slowed
  by the same 2.4 s on chr21 (oracle 5.0 → 3.9 M/s), which points to clock-read
  overhead, not to the data structures. The reader now coalesces FASTA sequence
  into 64 KiB chunks, and the bench decodes 4 MiB batches and times each structure
  once per batch. Same sandbox machine, chr21 insert:

  | Version | Table | Oracle |
  |---|---|---|
  | Before streaming (whole record) | 2.55 M/s | 7.29 M/s |
  | Streaming, per-line chunks | 2.51–2.62 M/s | 6.16 M/s |
  | Streaming, 64 KiB chunks, batch timing | **2.79 M/s** | 7.10 M/s |

## Real data (`tools/bench_real.sh`, laptop, before the update)

`results/{ecoli,chr21}_{lex,hash}.txt`.

### Throughput (k-mers/s, table / oracle)

| Dataset | Order | Insert | Positive query | Negative query |
|---|---|---|---|---|
| E. coli K-12 (4.6 Mb), 2^16 buckets | lex | 9.0 M / 8.0 M | 6.1 M / 8.0 M | 3.5 M / 15.2 M |
| E. coli K-12 | hash | 9.5 M / 8.2 M | 5.5 M / 8.0 M | 3.3 M / 15.2 M |
| chr21 GRCh38 (40 M k-mers), 2^20 buckets | lex | **0.58 M** / 5.1 M | **1.2 M** / 5.9 M | 2.7 M / 12.0 M |
| chr21 GRCh38 | hash | 2.8 M / 5.0 M | 3.0 M / 6.3 M | 2.5 M / 11.3 M |

### Memory (heap bytes per distinct k-mer)

| Dataset | Order | Pieces | Bucket headers | Entries | **Table total** | Oracle | Ideal u64 table at load 7/8 |
|---|---|---|---|---|---|---|---|
| E. coli (mean piece 8.6 kb) | lex | 1.01 | 0.35 | 1.78 | **3.14** | 16.58 | 9.14 |
| E. coli | hash | 1.01 | 0.35 | 1.51 | **2.86** | 16.58 | 9.14 |
| chr21 (mean piece 328 bp) | lex | 1.13 | 0.73 | 1.89 | **3.75** | 17.44 | 9.14 |
| chr21 | hash | 1.13 | 0.73 | 1.54 | **3.40** | 17.44 | 9.14 |

On chr21 with the hashed order, the table is 5.1× smaller than `HashSet<u64>` and
2.7× smaller than an ideal u64 table. Intra-run redundancy is 0.26 % of stored
k-mers on chr21 and 0.29 % on E. coli.

### Query cost and skew on chr21

| | lex | hash |
|---|---|---|
| Entries scanned per positive query | 10.5 | 5.2 |
| Candidate k-mers per positive query | 22.7 | 5.5 |
| Entries scanned per negative query | 5.2 | 4.2 |
| Negative queries reaching a string read | 2.1 % | 1.7 % |
| Largest minimizer (index entries) | 21 511 (poly-A) | 2 115 |
| Minimizers with ≥ 1024 entries | 9 | 7 |

* **The hashed order is clearly better on human sequence.** Insertion is 4.8×
  faster and positive queries are 2.5× faster. The index is 19 % smaller, because
  pieces have fewer tied minimizer positions.
* **Repeats, not the order, cause the remaining hot minimizers.** The hottest
  lexicographic minimizers include Alu fragments (`AAAGTGCTGGGATTACA`,
  `AATCCCAGCACTTTGGG`). The hashed order also has minimizers with thousands of
  entries. These buckets grow roughly in proportion to genome size, and chr21 is
  about 1.5 % of GRCh38.
* **Negative-query cost equals bucket occupancy.** Entries scanned per negative
  query is almost exactly entries per bucket: 4.2 on chr21 and 8.7 on E. coli,
  where 2^16 buckets is undersized for 571 k entries. Tags reject more than 99 %
  of them. Larger tables cost 24 B per extra bucket header, so a flat
  fixed-capacity layout is needed.
* **Peak RSS is 1.36 GB, but the table is only 118 MB.** The rest is the oracle
  (604 MB) and per-segment scratch space. Insertion currently keeps two `Vec<u64>`
  of length equal to the segment (ℓ-mer values and window minima), which is
  16 B per base of the longest ACGT segment. On chr21 that is hundreds of MB of
  temporary memory. Segments must be processed in chunks before running on whole
  genomes.

## Synthetic data (2-core cloud sandbox)

### Throughput (k-mers/s)

| Dataset | Order | Insert (table / oracle) | Positive query (table / oracle) | Negative query (table / oracle) |
|---|---|---|---|---|
| Uniform genome, 4.6 Mb, inserted whole | lex | 11.2 M / 14.8 M | 5.9 M / 6.9 M | 5.8 M / 21.8 M |
| 300k × 150 bp reads from it, 0.5 % errors, random strand | lex | 6.9 M / 8.8 M | 4.9 M / 6.5 M | 3.4 M / 14 M |
| same reads | hash | 7.2 M / 8.8 M | 4.8 M / 6.5 M | 3.0 M / 14 M |
| Repeat-rich genome, 4.6 Mb | lex | **0.045 M** / 19.9 M | **0.040 M** / 7.1 M | 3.8 M / 16 M |
| Repeat-rich genome, 4.6 Mb | hash | 3.7 M / 19.1 M | 2.5 M / 7.7 M | 3.7 M / 17 M |

### Memory (heap bytes per distinct k-mer, after `shrink_to_fit`)

| Dataset | Pieces | Index | **Table total** | Oracle `HashSet<u64>` | Ideal u64 table at load 7/8 |
|---|---|---|---|---|---|
| Uniform genome (mean piece 920 kb, 2^16 buckets) | 1.00 | 2.05 | **3.05** | 16.41 | 9.14 |
| Reads, lex (mean piece 68 bp, 2^18 buckets) | 2.01 | 2.72 | **4.73** | 17.69 | 9.14 |
| Reads, hash | 2.01 | 2.52 | **4.52** | 17.69 | 9.14 |
| Repeat-rich, lex (mean piece 541 bp, 2^20 buckets) | — | — | 10.24 | 10.49 | 9.14 |
| Repeat-rich, hash | — | — | 9.87 | 10.49 | 9.14 |

The repeat-rich rows use 2^20 buckets, so 25 MB of the total is `Vec` bucket headers
(24 B per bucket). The bucket count has to match the expected number of entries.
On the uniform genome with 2^20 buckets the total was 8.18 B/k-mer instead of 3.05.

### Observations

1. **The lexicographic minimizer order fails on low-complexity sequence.** On the
   repeat-rich genome the poly-A ℓ-mer alone has 27 k index entries, and queries scan
   about 100 entries and 750 candidate k-mers each, which makes insertion and queries
   60–80× slower than with the hashed order. The
   hashed order (a bijective mix of the canonical ℓ-mer, with the same tie structure)
   removes most of the skew. Homopolymer windows remain a hot spot under both orders
   (max bucket 9.7 k entries with hash), which hot-minimizer routing later
   addressed.
2. **The tie-offset candidate scan** (`--scan ties`, the default) checks only
   `s = q − p` for `p ∈ Q(x) ∪ mirror(Q(x))`, not the full `[q−(k−ℓ), q]` range
   range. Positive queries do 1.5 candidate comparisons instead of 8.3, and
   insertion is 27 % faster. Both modes are exact and both are tested against the
   oracle.
3. **Negative queries are 3–5× slower than the hash set.** They compute `m(x)`
   (O(k)) and touch a random bucket. Only about 1.5 % of negative queries reach a string
   comparison (0.08 candidate k-mers per query), so the 8-bit tag is sufficient.
4. **Memory.** When pieces are long, the table needs about 3 B/k-mer, against about
   16 B for `HashSet<u64>` and 9.1 B for an ideal u64 table. With 68 bp pieces (reads
   with errors) it needs 4.5–4.7 B/k-mer, which is still smaller. Index entries take
   12 B each (unpacked struct, 3 B of padding), and there are about 6 per piece.
