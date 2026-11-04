# Paper measurements on compiu-pro (2026-11-02T19:30:26+01:00)

Machine: see machine.txt. Commit: 8993d8a.

| paper table / text | files |
|---|---|
| memory, single-thread throughput, per-query counters | kache/<data>_{summary.md,unihash.txt,kache.txt} |
| concurrent insertion, concurrent queries | kache/<data>_scale_{summary.md,unihash.txt,kache.txt} |
| ConcurrentKmerTable scaling and deletion, KmerTable deletion | scale/<data>_capacity.txt |
| eager vs deferred run indexing | eager/summary.txt |
| hot threshold, candidate scan, max_piece_len | ablation/summary.txt |
| defrag after bulk loading | defrag/*.txt |
| churn | churn/<data>.txt |

Not re-run by this script: the lexicographic minimizer order (no longer an
option of the code), the 64-byte block size sweep and the A/B of the
insertion shortcuts (comparisons against earlier commits), and the
concurrent in-place defrag timings (sandbox).

## Step timings

```
tests         57 s  finished 2026-11-02T17:30:26+01:00
kache       1648 s  finished 2026-11-02T18:43:49+01:00
scale       1121 s  finished 2026-11-02T19:02:30+01:00
eager        320 s  finished 2026-11-02T19:07:50+01:00
ablation     251 s  finished 2026-11-02T19:12:01+01:00
defrag       193 s  finished 2026-11-02T19:15:14+01:00
churn        912 s  finished 2026-11-02T19:30:26+01:00
```

## chr21_scale_summary.md

Insertion, M k-mers/s (seconds); every structure holds the same 34625350 distinct k-mers.

| threads | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| 1, single-threaded structure | 7.51 (5.34 s) | 7.56 (5.30 s) | 8.31 (4.82 s) |
| 1, concurrent structure | 7.88 (5.09 s) | 11.65 (3.44 s) | 12.86 (3.12 s) |
| 2, concurrent structure | 14.39 (2.79 s) | 22.70 (1.77 s) | 23.36 (1.72 s) |
| 4, concurrent structure | 27.93 (1.44 s) | 41.24 (0.97 s) | 37.59 (1.07 s) |
| 8, concurrent structure | 44.70 (0.90 s) | 63.75 (0.63 s) | 44.16 (0.91 s) |
| 12, concurrent structure | 48.75 (0.82 s) | 75.01 (0.53 s) | 44.75 (0.90 s) |

Memory, B per distinct k-mer (heap after the build):

| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| single-threaded structure | 1.927 | 23.500 | 17.443 |
| concurrent structure (12 threads) | 2.450 | 23.500 | 17.445 |
| unihash concurrent table after `into_table` | 1.866 | | |

Queries, positive k-mers (shuffled), M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 2.13 | 1.99 | 2.02 | 10.66 |
| 2 | 4.25 | 3.58 | 4.11 | 21.17 |
| 4 | 8.45 | 6.77 | 7.94 | 37.29 |
| 8 | 13.45 | 11.51 | 12.77 | 46.80 |
| 12 | 15.85 | 14.89 | 17.54 | 57.71 |

Queries, uniform random k-mers, M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 3.50 | 3.60 | 2.30 | 18.43 |
| 2 | 6.95 | 6.37 | 4.58 | 36.37 |
| 4 | 13.50 | 12.03 | 8.29 | 69.28 |
| 8 | 19.02 | 17.88 | 14.31 | 80.30 |
| 12 | 22.04 | 22.30 | 19.14 | 102.22 |

Queries, streaming (sequences), M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 36.65 | 27.52 | 22.78 | 12.93 |
| 2 | 64.60 | 44.78 | 37.66 | 21.88 |
| 4 | 127.36 | 82.29 | 71.25 | 40.61 |
| 8 | 163.79 | 127.26 | 106.88 | 52.82 |
| 12 | 221.37 | 144.75 | 89.58 | 61.39 |

Single-threaded structures grow from empty (kache-hash from its initial capacity); the concurrent ones are created for the distinct k-mer count (unihash `with_capacity`, kache-hash's `max_sz`, `HashSet::with_capacity` per shard). KmerTable, kache-hash and HashSet are queried after the single-threaded build, shared read-only (no locks); ConcurrentKmerTable takes its shard read locks. The concurrent HashSet is 1024 shards of Mutex<HashSet<u64>> filled through per-thread buffers. HashSet streaming: independent lookups of the k-mers parsed from the same sequences.

## chr21_summary.md

| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| Memory, B per distinct k-mer | 1.93 | 23.50 | 17.44 |
| Insert, M k-mers/s | 7.10 | 9.03 | 7.86 |
| inserted input (positives): k-mers in order | 4.49 | 2.99 | 9.79 |
| inserted input (positives): k-mers random order | 2.28 | 2.04 | 9.77 |
| inserted input (positives): sequence, streaming | 29.10 | 19.11 | 5.09 |
| Random negatives, M k-mers/s | 3.72 | 2.31 | 16.32 |

Streaming row, HashSet column: independent lookups parsed from the same sequence.

## ecoli_scale_summary.md

Insertion, M k-mers/s (seconds); every structure holds the same 4554269 distinct k-mers.

| threads | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| 1, single-threaded structure | 10.94 (0.42 s) | 11.12 (0.42 s) | 11.81 (0.39 s) |
| 1, concurrent structure | 11.38 (0.41 s) | 14.26 (0.33 s) | 18.96 (0.24 s) |
| 2, concurrent structure | 23.96 (0.19 s) | 27.73 (0.17 s) | 33.50 (0.14 s) |
| 4, concurrent structure | 47.28 (0.10 s) | 54.06 (0.09 s) | 48.85 (0.10 s) |
| 8, concurrent structure | 69.94 (0.07 s) | 84.73 (0.06 s) | 54.33 (0.09 s) |
| 12, concurrent structure | 71.13 (0.07 s) | 81.69 (0.06 s) | 53.50 (0.09 s) |

Memory, B per distinct k-mer (heap after the build):

| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| single-threaded structure | 1.754 | 22.335 | 16.577 |
| concurrent structure (12 threads) | 2.169 | 22.335 | 16.593 |
| unihash concurrent table after `into_table` | 1.696 | | |

Queries, positive k-mers (shuffled), M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 4.38 | 3.41 | 2.21 | 11.88 |
| 2 | 9.44 | 7.12 | 4.47 | 23.28 |
| 4 | 19.08 | 14.93 | 8.98 | 46.38 |
| 8 | 26.24 | 21.09 | 14.04 | 50.48 |
| 12 | 35.49 | 28.35 | 19.40 | 70.40 |

Queries, uniform random k-mers, M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 5.16 | 4.51 | 2.46 | 28.83 |
| 2 | 10.74 | 9.33 | 4.90 | 64.57 |
| 4 | 22.23 | 18.51 | 9.83 | 131.76 |
| 8 | 30.90 | 26.01 | 15.12 | 172.58 |
| 12 | 39.83 | 32.72 | 20.33 | 255.25 |

Queries, streaming (sequences), M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 44.60 | 31.68 | 24.21 | 12.09 |
| 2 | 88.20 | 62.90 | 47.29 | 24.58 |
| 4 | 174.24 | 122.69 | 90.48 | 47.34 |
| 8 | 217.58 | 147.13 | 127.89 | 49.52 |
| 12 | 305.38 | 218.89 | 174.19 | 61.08 |

Single-threaded structures grow from empty (kache-hash from its initial capacity); the concurrent ones are created for the distinct k-mer count (unihash `with_capacity`, kache-hash's `max_sz`, `HashSet::with_capacity` per shard). KmerTable, kache-hash and HashSet are queried after the single-threaded build, shared read-only (no locks); ConcurrentKmerTable takes its shard read locks. The concurrent HashSet is 1024 shards of Mutex<HashSet<u64>> filled through per-thread buffers. HashSet streaming: independent lookups of the k-mers parsed from the same sequences.

## ecoli_summary.md

| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| Memory, B per distinct k-mer | 1.76 | 22.34 | 16.58 |
| Insert, M k-mers/s | 9.57 | 10.92 | 11.90 |
| inserted input (positives): k-mers in order | 5.58 | 3.06 | 12.05 |
| inserted input (positives): k-mers random order | 4.03 | 2.14 | 11.97 |
| inserted input (positives): sequence, streaming | 40.16 | 24.08 | 6.46 |
| Random negatives, M k-mers/s | 4.93 | 2.41 | 25.92 |

Streaming row, HashSet column: independent lookups parsed from the same sequence.

## SRR5833294_scale_summary.md

Insertion, M k-mers/s (seconds); every structure holds the same 80209364 distinct k-mers.

| threads | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| 1, single-threaded structure | 5.28 (17.40 s) | 7.08 (12.99 s) | 6.04 (15.21 s) |
| 1, concurrent structure | 6.45 (14.26 s) | 9.39 (9.79 s) | 9.63 (9.55 s) |
| 2, concurrent structure | 11.95 (7.69 s) | 18.22 (5.05 s) | 18.04 (5.09 s) |
| 4, concurrent structure | 24.07 (3.82 s) | 32.53 (2.83 s) | 29.45 (3.12 s) |
| 8, concurrent structure | 36.85 (2.50 s) | 42.37 (2.17 s) | 33.69 (2.73 s) |
| 12, concurrent structure | 43.99 (2.09 s) | 46.95 (1.96 s) | 34.08 (2.70 s) |

Memory, B per distinct k-mer (heap after the build):

| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| single-threaded structure | 2.506 | 20.289 | 15.060 |
| concurrent structure (12 threads) | 2.797 | 20.289 | 15.061 |
| unihash concurrent table after `into_table` | 2.428 | | |

Queries, positive k-mers (shuffled), M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 1.73 | 1.50 | 1.75 | 7.98 |
| 2 | 3.44 | 2.80 | 3.57 | 15.79 |
| 4 | 6.64 | 5.52 | 7.09 | 29.53 |
| 8 | 10.65 | 8.64 | 11.45 | 36.27 |
| 12 | 14.09 | 11.44 | 15.03 | 41.35 |

Queries, uniform random k-mers, M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 3.40 | 2.83 | 2.03 | 14.36 |
| 2 | 6.78 | 5.67 | 3.97 | 27.54 |
| 4 | 13.07 | 11.36 | 8.10 | 50.46 |
| 8 | 18.18 | 16.11 | 12.75 | 59.02 |
| 12 | 21.85 | 19.14 | 16.80 | 79.77 |

Queries, streaming (sequences), M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 19.85 | 14.91 | 13.54 | 7.83 |
| 2 | 38.94 | 28.44 | 25.40 | 15.54 |
| 4 | 76.36 | 56.75 | 49.05 | 29.41 |
| 8 | 108.55 | 76.89 | 75.53 | 35.44 |
| 12 | 151.96 | 99.99 | 102.95 | 44.66 |

Single-threaded structures grow from empty (kache-hash from its initial capacity); the concurrent ones are created for the distinct k-mer count (unihash `with_capacity`, kache-hash's `max_sz`, `HashSet::with_capacity` per shard). KmerTable, kache-hash and HashSet are queried after the single-threaded build, shared read-only (no locks); ConcurrentKmerTable takes its shard read locks. The concurrent HashSet is 1024 shards of Mutex<HashSet<u64>> filled through per-thread buffers. HashSet streaming: independent lookups of the k-mers parsed from the same sequences.

## SRR5833294_summary.md

| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| Memory, B per distinct k-mer | 2.51 | 20.29 | 15.06 |
| Insert, M k-mers/s | 5.48 | 7.40 | 7.07 |
| inserted input (positives): k-mers in order | 4.29 | 2.79 | 8.59 |
| inserted input (positives): k-mers random order | 1.83 | 1.76 | 8.65 |
| inserted input (positives): sequence, streaming | 19.93 | 14.88 | 4.46 |
| held-out records 2000001..4000000 of each input file: k-mers in order | 4.41 | 2.97 | 12.68 |
| held-out records 2000001..4000000 of each input file: k-mers random order | 2.33 | 2.03 | 11.44 |
| held-out records 2000001..4000000 of each input file: sequence, streaming | 11.15 | 16.15 | 7.24 |
| query set data/SRR5901135_1.fastq.gz: k-mers in order | 5.84 | 3.04 | 15.98 |
| query set data/SRR5901135_1.fastq.gz: k-mers random order | 3.76 | 2.21 | 15.52 |
| query set data/SRR5901135_1.fastq.gz: sequence, streaming | 18.13 | 22.95 | 9.38 |
| query set data/SRR5901135_2.fastq.gz: k-mers in order | 5.81 | 3.07 | 15.46 |
| query set data/SRR5901135_2.fastq.gz: k-mers random order | 3.67 | 2.10 | 15.55 |
| query set data/SRR5901135_2.fastq.gz: sequence, streaming | 17.77 | 21.95 | 9.32 |
| Random negatives, M k-mers/s | 3.64 | 2.06 | 15.39 |

Streaming row, HashSet column: independent lookups parsed from the same sequence.

## eager/summary.txt

```
ecoli       deferred rep 1  insert 0.446 s 10.40 M/s  1.756 B/k-mer  pieces 1643  stored 4557029  redundant 2760
ecoli       eager    rep 1  insert 0.382 s 12.14 M/s  1.756 B/k-mer  pieces 1724  stored 4554269  redundant 0
ecoli       deferred rep 2  insert 0.426 s 10.90 M/s  1.756 B/k-mer  pieces 1643  stored 4557029  redundant 2760
ecoli       eager    rep 2  insert 0.382 s 12.15 M/s  1.756 B/k-mer  pieces 1724  stored 4554269  redundant 0
ecoli       deferred rep 3  insert 0.423 s 10.99 M/s  1.756 B/k-mer  pieces 1643  stored 4557029  redundant 2760
ecoli       eager    rep 3  insert 0.383 s 12.11 M/s  1.756 B/k-mer  pieces 1724  stored 4554269  redundant 0
chr21       deferred rep 1  insert 5.520 s 7.26 M/s  1.929 B/k-mer  pieces 119490  stored 34706358  redundant 81008
chr21       eager    rep 1  insert 5.052 s 7.94 M/s  1.927 B/k-mer  pieces 122335  stored 34625350  redundant 0
chr21       deferred rep 2  insert 5.543 s 7.23 M/s  1.929 B/k-mer  pieces 119490  stored 34706358  redundant 81008
chr21       eager    rep 2  insert 5.050 s 7.94 M/s  1.927 B/k-mer  pieces 122335  stored 34625350  redundant 0
chr21       deferred rep 3  insert 5.498 s 7.29 M/s  1.929 B/k-mer  pieces 119490  stored 34706358  redundant 81008
chr21       eager    rep 3  insert 5.054 s 7.93 M/s  1.927 B/k-mer  pieces 122335  stored 34625350  redundant 0
SRR5833294  deferred rep 1  insert 17.427 s 5.27 M/s  2.505 B/k-mer  pieces 2069193  stored 80215174  redundant 5810
SRR5833294  eager    rep 1  insert 15.861 s 5.80 M/s  2.506 B/k-mer  pieces 2069503  stored 80209364  redundant 0
SRR5833294  deferred rep 2  insert 17.408 s 5.28 M/s  2.505 B/k-mer  pieces 2069193  stored 80215174  redundant 5810
SRR5833294  eager    rep 2  insert 15.868 s 5.79 M/s  2.506 B/k-mer  pieces 2069503  stored 80209364  redundant 0
SRR5833294  deferred rep 3  insert 17.431 s 5.27 M/s  2.505 B/k-mer  pieces 2069193  stored 80215174  redundant 5810
SRR5833294  eager    rep 3  insert 15.906 s 5.78 M/s  2.506 B/k-mer  pieces 2069503  stored 80209364  redundant 0
```

## ablation/summary.txt

```
chr21 T=0  insert 6.836 s 5.86 M/s  1.877 B/k-mer  pieces 122335  entries 4436747  hot 0/0  in-order 2.75  random 1.73  streaming 15.01  negatives 4.07  candidates 11.976
chr21 T=16  insert 5.019 s 7.99 M/s  2.061 B/k-mer  pieces 122335  entries 4214995  hot 4564/2055  in-order 4.95  random 2.77  streaming 33.85  negatives 4.07  candidates 1.839
chr21 T=64  insert 5.109 s 7.85 M/s  1.927 B/k-mer  pieces 122335  entries 4323807  hot 691/97  in-order 4.70  random 2.62  streaming 31.53  negatives 4.03  candidates 2.969
chr21 T=256  insert 5.620 s 7.13 M/s  1.894 B/k-mer  pieces 122335  entries 4388212  hot 92/1  in-order 4.20  random 2.37  streaming 27.19  negatives 4.06  candidates 4.958
ecoli T=0  insert 0.358 s 12.96 M/s  1.754 B/k-mer  pieces 1724  entries 570887  hot 0/0  in-order 5.99  random 4.91  streaming 40.81  negatives 5.86  candidates 1.549
ecoli T=64  insert 0.383 s 12.12 M/s  1.756 B/k-mer  pieces 1724  entries 570803  hot 1/0  in-order 5.86  random 4.89  streaming 40.74  negatives 5.86  candidates 1.545
reads scan=ties  insert 2.215 s 17.15 M/s  2.600 B/k-mer  pieces 251124  entries 1237363  hot 0/0  in-order 5.19  random 3.29  streaming 28.25  negatives 4.71  candidates 1.838
reads scan=range  insert 2.503 s 15.18 M/s  2.600 B/k-mer  pieces 251124  entries 1237363  hot 0/0  in-order 4.54  random 2.95  streaming 27.34  negatives 4.64  candidates 10.848
chr21 max_piece_len=4096  insert 5.073 s 7.90 M/s  1.927 B/k-mer  pieces 122335  entries 4323807  hot 691/97  in-order 4.71  random 2.62  streaming 31.57  negatives 4.01  candidates 2.969
chr21 max_piece_len=1048576  insert 5.069 s 7.91 M/s  1.926 B/k-mer  pieces 120106  entries 4321877  hot 691/97  in-order 4.73  random 2.64  streaming 31.66  negatives 4.06  candidates 2.927
```
