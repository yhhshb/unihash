# Paper measurements on compiu-pro (2026-11-04T20:25:38+01:00)

Machine: see machine.txt. Commit: fac7269.

| paper table / text | files |
|---|---|
| memory, single-thread throughput, per-query counters | kache/<data>_{summary.md,unihash.txt,kache.txt} |
| concurrent insertion, concurrent queries | kache/<data>_scale_{summary.md,unihash.txt,kache.txt} |
| ConcurrentKmerTable scaling and deletion, KmerTable deletion | scale/<data>_capacity.txt |
| cost of the read locks (locked queries vs read snapshot, same table) | scale/<data>_capacity.txt, section "read locks" |
| concurrent defrag in place vs round trip through KmerTable | scale/<data>_capacity.txt, section "defragmentation after deletion" |
| eager vs deferred run indexing | eager/summary.txt |
| hot threshold, candidate scan, max_piece_len | ablation/summary.txt |
| defrag after bulk loading | defrag/*.txt |
| churn | churn/<data>.txt |

Not re-run by this script: the lexicographic minimizer order (no longer an
option of the code), the 64-byte block size sweep and the A/B of the
insertion shortcuts (comparisons against earlier commits), and the
lock-free read path that was tried and dropped (not in the code; the
read-lock section gives the most any lock-free read path could gain).

## Step timings

```
tests         59 s  finished 2026-11-04T18:24:16+01:00
kache       2733 s  finished 2026-11-04T19:09:49+01:00
scale       2866 s  finished 2026-11-04T19:57:35+01:00
eager        319 s  finished 2026-11-04T20:02:54+01:00
ablation     250 s  finished 2026-11-04T20:07:04+01:00
defrag       193 s  finished 2026-11-04T20:10:17+01:00
churn        921 s  finished 2026-11-04T20:25:38+01:00
```

## chr21_scale_summary.md

Insertion, M k-mers/s (seconds); every structure holds the same 34625350 distinct k-mers.

| threads | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| 1, single-threaded structure | 7.84 (5.12 s) | 9.02 (4.45 s) | 8.85 (4.53 s) |
| 1, concurrent structure | 8.23 (4.87 s) | 12.74 (3.15 s) | 13.96 (2.87 s) |
| 2, concurrent structure | 16.08 (2.49 s) | 24.20 (1.66 s) | 25.27 (1.59 s) |
| 3, concurrent structure | 24.22 (1.66 s) | 34.76 (1.15 s) | 34.57 (1.16 s) |
| 4, concurrent structure | 31.83 (1.26 s) | 43.16 (0.93 s) | 41.90 (0.96 s) |
| 5, concurrent structure | 38.96 (1.03 s) | 51.57 (0.78 s) | 46.17 (0.87 s) |
| 6, concurrent structure | 47.16 (0.85 s) | 59.05 (0.68 s) | 48.68 (0.82 s) |
| 7, concurrent structure | 47.30 (0.85 s) | 62.93 (0.64 s) | 49.55 (0.81 s) |
| 8, concurrent structure | 49.58 (0.81 s) | 67.36 (0.59 s) | 49.90 (0.80 s) |
| 9, concurrent structure | 51.28 (0.78 s) | 71.36 (0.56 s) | 50.36 (0.80 s) |
| 10, concurrent structure | 52.49 (0.76 s) | 73.61 (0.55 s) | 50.82 (0.79 s) |
| 11, concurrent structure | 53.48 (0.75 s) | 77.60 (0.52 s) | 50.69 (0.79 s) |
| 12, concurrent structure | 54.95 (0.73 s) | 79.78 (0.50 s) | 51.05 (0.79 s) |

Memory, B per distinct k-mer (heap after the build):

| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| single-threaded structure | 1.927 | 23.500 | 17.443 |
| concurrent structure (12 threads) | 2.449 | 23.500 | 17.445 |
| unihash concurrent table after `into_table` | 1.866 | | |

Queries, positive k-mers (shuffled), M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 2.50 | 2.14 | 2.17 | 11.44 |
| 2 | 5.05 | 4.24 | 4.33 | 22.51 |
| 3 | 7.62 | 6.44 | 6.53 | 32.70 |
| 4 | 10.19 | 8.56 | 8.70 | 42.55 |
| 5 | 12.74 | 10.57 | 10.85 | 51.86 |
| 6 | 15.23 | 12.57 | 12.98 | 60.78 |
| 7 | 13.74 | 11.75 | 12.05 | 46.10 |
| 8 | 15.63 | 13.40 | 13.71 | 50.56 |
| 9 | 17.32 | 15.02 | 15.25 | 54.59 |
| 10 | 18.39 | 15.92 | 16.22 | 53.59 |
| 11 | 20.08 | 17.51 | 17.77 | 58.22 |
| 12 | 21.91 | 19.17 | 19.34 | 62.76 |

Queries, uniform random k-mers, M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 3.75 | 3.40 | 2.42 | 19.66 |
| 2 | 7.56 | 6.72 | 4.85 | 38.59 |
| 3 | 11.35 | 10.10 | 7.27 | 57.07 |
| 4 | 15.18 | 13.45 | 9.71 | 73.89 |
| 5 | 18.85 | 16.62 | 12.09 | 90.24 |
| 6 | 22.44 | 19.85 | 14.42 | 105.50 |
| 7 | 19.09 | 17.28 | 13.20 | 76.91 |
| 8 | 21.60 | 19.62 | 14.95 | 85.84 |
| 9 | 23.91 | 21.84 | 16.72 | 89.58 |
| 10 | 24.28 | 22.60 | 17.65 | 98.96 |
| 11 | 26.62 | 24.72 | 19.33 | 106.67 |
| 12 | 28.90 | 26.86 | 21.02 | 115.05 |

Queries, streaming (sequences), M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 40.83 | 28.96 | 24.08 | 13.70 |
| 2 | 72.14 | 50.64 | 40.00 | 22.74 |
| 3 | 106.86 | 74.57 | 59.67 | 33.03 |
| 4 | 140.16 | 99.17 | 80.69 | 43.20 |
| 5 | 164.01 | 118.76 | 96.93 | 52.55 |
| 6 | 204.26 | 139.82 | 119.26 | 61.78 |
| 7 | 170.87 | 128.91 | 118.29 | 59.16 |
| 8 | 173.43 | 139.48 | 123.45 | 54.92 |
| 9 | 198.13 | 132.52 | 134.83 | 57.38 |
| 10 | 194.53 | 135.66 | 144.37 | 60.31 |
| 11 | 208.11 | 145.09 | 154.48 | 65.52 |
| 12 | 227.69 | 156.53 | 163.57 | 65.44 |

Single-threaded structures grow from empty (kache-hash from its initial capacity); the concurrent ones are created for the distinct k-mer count (unihash `with_capacity`, kache-hash's `max_sz`, `HashSet::with_capacity` per shard). KmerTable, kache-hash and HashSet are queried after the single-threaded build, shared read-only (no locks); ConcurrentKmerTable takes its shard read locks. The concurrent HashSet is 1024 shards of Mutex<HashSet<u64>> filled through per-thread buffers. HashSet streaming: independent lookups of the k-mers parsed from the same sequences.

## chr21_summary.md

| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| Memory, B per distinct k-mer | 1.93 | 23.50 | 17.44 |
| Insert, M k-mers/s | 7.80 | 9.17 | 8.61 |
| inserted input (positives): k-mers in order | 4.69 | 3.06 | 11.03 |
| inserted input (positives): k-mers random order | 2.57 | 2.07 | 10.76 |
| inserted input (positives): sequence, streaming | 31.52 | 20.69 | 5.82 |
| Random negatives, M k-mers/s | 4.02 | 2.39 | 17.99 |

Streaming row, HashSet column: independent lookups parsed from the same sequence.

## ecoli_scale_summary.md

Insertion, M k-mers/s (seconds); every structure holds the same 4554269 distinct k-mers.

| threads | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| 1, single-threaded structure | 12.24 (0.38 s) | 11.39 (0.41 s) | 14.06 (0.33 s) |
| 1, concurrent structure | 12.87 (0.36 s) | 13.94 (0.33 s) | 22.00 (0.21 s) |
| 2, concurrent structure | 24.76 (0.19 s) | 25.29 (0.18 s) | 35.68 (0.13 s) |
| 3, concurrent structure | 36.58 (0.13 s) | 36.77 (0.13 s) | 49.07 (0.10 s) |
| 4, concurrent structure | 48.22 (0.10 s) | 46.22 (0.10 s) | 57.77 (0.08 s) |
| 5, concurrent structure | 58.61 (0.08 s) | 52.37 (0.09 s) | 61.32 (0.08 s) |
| 6, concurrent structure | 71.80 (0.07 s) | 57.83 (0.08 s) | 65.38 (0.07 s) |
| 7, concurrent structure | 71.74 (0.07 s) | 58.50 (0.08 s) | 67.03 (0.07 s) |
| 8, concurrent structure | 71.63 (0.07 s) | 61.19 (0.08 s) | 67.74 (0.07 s) |
| 9, concurrent structure | 71.39 (0.07 s) | 63.56 (0.07 s) | 69.48 (0.07 s) |
| 10, concurrent structure | 71.36 (0.07 s) | 63.81 (0.07 s) | 69.02 (0.07 s) |
| 11, concurrent structure | 75.83 (0.06 s) | 61.13 (0.08 s) | 69.07 (0.07 s) |
| 12, concurrent structure | 75.80 (0.06 s) | 58.53 (0.08 s) | 69.57 (0.07 s) |

Memory, B per distinct k-mer (heap after the build):

| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| single-threaded structure | 1.754 | 22.335 | 16.577 |
| concurrent structure (12 threads) | 2.169 | 22.335 | 16.594 |
| unihash concurrent table after `into_table` | 1.696 | | |

Queries, positive k-mers (shuffled), M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 4.63 | 4.01 | 2.36 | 14.83 |
| 2 | 9.23 | 7.69 | 4.75 | 29.66 |
| 3 | 13.79 | 11.24 | 7.17 | 44.41 |
| 4 | 18.39 | 15.01 | 9.55 | 59.27 |
| 5 | 22.96 | 18.46 | 11.98 | 73.37 |
| 6 | 27.58 | 22.07 | 14.35 | 86.94 |
| 7 | 21.63 | 18.17 | 12.86 | 55.93 |
| 8 | 24.41 | 20.44 | 14.60 | 62.92 |
| 9 | 25.29 | 21.92 | 16.09 | 69.89 |
| 10 | 27.88 | 23.59 | 17.48 | 76.48 |
| 11 | 30.68 | 25.72 | 19.21 | 83.46 |
| 12 | 33.47 | 27.95 | 20.91 | 90.41 |

Queries, uniform random k-mers, M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 5.33 | 4.77 | 2.57 | 41.29 |
| 2 | 10.57 | 9.16 | 5.15 | 81.52 |
| 3 | 15.86 | 13.30 | 7.71 | 122.38 |
| 4 | 21.14 | 17.63 | 10.27 | 163.43 |
| 5 | 26.42 | 21.75 | 12.83 | 205.28 |
| 6 | 31.79 | 26.10 | 15.31 | 245.17 |
| 7 | 25.09 | 21.83 | 13.76 | 188.37 |
| 8 | 28.47 | 24.46 | 15.64 | 203.41 |
| 9 | 30.94 | 26.80 | 17.40 | 227.28 |
| 10 | 30.65 | 27.44 | 18.19 | 250.90 |
| 11 | 33.64 | 29.95 | 19.96 | 275.70 |
| 12 | 36.69 | 32.59 | 21.70 | 300.60 |

Queries, streaming (sequences), M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 46.15 | 32.32 | 24.59 | 15.00 |
| 2 | 90.86 | 63.51 | 48.34 | 29.80 |
| 3 | 135.89 | 95.33 | 71.69 | 44.40 |
| 4 | 180.82 | 127.02 | 94.91 | 59.04 |
| 5 | 216.17 | 149.80 | 113.69 | 70.38 |
| 6 | 269.35 | 188.24 | 142.20 | 86.24 |
| 7 | 229.36 | 155.75 | 138.38 | 68.64 |
| 8 | 220.05 | 146.18 | 134.73 | 64.12 |
| 9 | 246.86 | 164.02 | 150.00 | 70.20 |
| 10 | 246.89 | 163.74 | 152.04 | 69.78 |
| 11 | 281.23 | 186.93 | 166.47 | 78.20 |
| 12 | 327.71 | 215.03 | 194.60 | 89.92 |

Single-threaded structures grow from empty (kache-hash from its initial capacity); the concurrent ones are created for the distinct k-mer count (unihash `with_capacity`, kache-hash's `max_sz`, `HashSet::with_capacity` per shard). KmerTable, kache-hash and HashSet are queried after the single-threaded build, shared read-only (no locks); ConcurrentKmerTable takes its shard read locks. The concurrent HashSet is 1024 shards of Mutex<HashSet<u64>> filled through per-thread buffers. HashSet streaming: independent lookups of the k-mers parsed from the same sequences.

## ecoli_summary.md

| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| Memory, B per distinct k-mer | 1.76 | 22.34 | 16.58 |
| Insert, M k-mers/s | 12.13 | 11.85 | 13.61 |
| inserted input (positives): k-mers in order | 5.95 | 3.18 | 13.50 |
| inserted input (positives): k-mers random order | 4.90 | 2.34 | 13.39 |
| inserted input (positives): sequence, streaming | 40.78 | 24.97 | 8.02 |
| Random negatives, M k-mers/s | 5.71 | 2.51 | 26.96 |

Streaming row, HashSet column: independent lookups parsed from the same sequence.

## SRR5833294_scale_summary.md

Insertion, M k-mers/s (seconds); every structure holds the same 80209364 distinct k-mers.

| threads | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| 1, single-threaded structure | 5.65 (16.27 s) | 8.63 (10.65 s) | 7.62 (12.07 s) |
| 1, concurrent structure | 6.95 (13.22 s) | 11.31 (8.13 s) | 11.87 (7.75 s) |
| 2, concurrent structure | 13.86 (6.63 s) | 21.83 (4.21 s) | 22.04 (4.17 s) |
| 3, concurrent structure | 20.74 (4.43 s) | 31.38 (2.93 s) | 30.74 (2.99 s) |
| 4, concurrent structure | 27.58 (3.33 s) | 40.12 (2.29 s) | 37.92 (2.42 s) |
| 5, concurrent structure | 34.28 (2.68 s) | 46.78 (1.97 s) | 42.72 (2.15 s) |
| 6, concurrent structure | 41.18 (2.23 s) | 53.00 (1.73 s) | 45.40 (2.02 s) |
| 7, concurrent structure | 42.94 (2.14 s) | 54.30 (1.69 s) | 45.84 (2.00 s) |
| 8, concurrent structure | 44.77 (2.05 s) | 54.50 (1.69 s) | 46.25 (1.99 s) |
| 9, concurrent structure | 46.52 (1.98 s) | 64.97 (1.42 s) | 45.94 (2.00 s) |
| 10, concurrent structure | 48.11 (1.91 s) | 64.77 (1.42 s) | 45.47 (2.02 s) |
| 11, concurrent structure | 49.57 (1.85 s) | 65.24 (1.41 s) | 47.23 (1.95 s) |
| 12, concurrent structure | 51.43 (1.79 s) | 61.86 (1.49 s) | 47.44 (1.94 s) |

Memory, B per distinct k-mer (heap after the build):

| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| single-threaded structure | 2.506 | 20.289 | 15.060 |
| concurrent structure (12 threads) | 2.798 | 20.289 | 15.061 |
| unihash concurrent table after `into_table` | 2.427 | | |

Queries, positive k-mers (shuffled), M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 2.02 | 1.69 | 2.13 | 10.91 |
| 2 | 4.07 | 3.43 | 4.28 | 21.72 |
| 3 | 6.14 | 5.15 | 6.41 | 31.74 |
| 4 | 8.06 | 6.84 | 8.54 | 40.89 |
| 5 | 10.01 | 8.40 | 10.61 | 49.45 |
| 6 | 11.83 | 9.98 | 12.64 | 57.92 |
| 7 | 11.25 | 9.80 | 11.85 | 44.02 |
| 8 | 12.75 | 11.15 | 13.48 | 49.24 |
| 9 | 14.23 | 12.31 | 15.01 | 52.61 |
| 10 | 15.26 | 13.44 | 15.98 | 51.50 |
| 11 | 16.53 | 14.56 | 17.51 | 55.61 |
| 12 | 17.94 | 15.96 | 19.03 | 59.93 |

Queries, uniform random k-mers, M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 3.60 | 3.25 | 2.39 | 18.09 |
| 2 | 7.28 | 6.45 | 4.79 | 35.97 |
| 3 | 10.85 | 9.73 | 7.20 | 52.55 |
| 4 | 14.51 | 12.86 | 9.53 | 63.48 |
| 5 | 18.03 | 16.02 | 11.91 | 82.61 |
| 6 | 21.35 | 19.07 | 14.16 | 96.57 |
| 7 | 18.25 | 16.75 | 13.00 | 71.61 |
| 8 | 20.57 | 18.90 | 14.75 | 78.22 |
| 9 | 22.73 | 20.95 | 16.42 | 83.30 |
| 10 | 23.30 | 21.78 | 17.31 | 90.73 |
| 11 | 25.40 | 23.76 | 18.98 | 97.93 |
| 12 | 27.58 | 26.01 | 20.62 | 105.54 |

Queries, streaming (sequences), M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 23.16 | 16.91 | 16.96 | 10.77 |
| 2 | 44.96 | 32.54 | 32.65 | 21.40 |
| 3 | 67.13 | 48.65 | 47.32 | 31.32 |
| 4 | 88.63 | 64.90 | 62.61 | 40.28 |
| 5 | 110.71 | 80.52 | 76.03 | 48.77 |
| 6 | 133.49 | 95.67 | 90.98 | 57.35 |
| 7 | 104.35 | 77.75 | 81.64 | 43.84 |
| 8 | 120.66 | 86.57 | 94.11 | 48.63 |
| 9 | 134.85 | 95.08 | 103.23 | 50.87 |
| 10 | 149.37 | 105.71 | 112.74 | 50.85 |
| 11 | 163.34 | 116.27 | 122.38 | 55.38 |
| 12 | 177.14 | 125.94 | 131.99 | 59.64 |

Single-threaded structures grow from empty (kache-hash from its initial capacity); the concurrent ones are created for the distinct k-mer count (unihash `with_capacity`, kache-hash's `max_sz`, `HashSet::with_capacity` per shard). KmerTable, kache-hash and HashSet are queried after the single-threaded build, shared read-only (no locks); ConcurrentKmerTable takes its shard read locks. The concurrent HashSet is 1024 shards of Mutex<HashSet<u64>> filled through per-thread buffers. HashSet streaming: independent lookups of the k-mers parsed from the same sequences.

## SRR5833294_summary.md

| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| Memory, B per distinct k-mer | 2.51 | 20.29 | 15.06 |
| Insert, M k-mers/s | 5.61 | 8.77 | 7.77 |
| inserted input (positives): k-mers in order | 4.65 | 3.05 | 10.60 |
| inserted input (positives): k-mers random order | 2.08 | 2.11 | 10.60 |
| inserted input (positives): sequence, streaming | 22.05 | 17.05 | 5.57 |
| held-out records 2000001..4000000 of each input file: k-mers in order | 4.70 | 3.13 | 15.04 |
| held-out records 2000001..4000000 of each input file: k-mers random order | 2.59 | 2.31 | 13.76 |
| held-out records 2000001..4000000 of each input file: sequence, streaming | 12.01 | 18.42 | 8.45 |
| query set data/SRR5901135_1.fastq.gz: k-mers in order | 6.03 | 3.26 | 15.43 |
| query set data/SRR5901135_1.fastq.gz: k-mers random order | 4.07 | 2.49 | 17.09 |
| query set data/SRR5901135_1.fastq.gz: sequence, streaming | 19.26 | 25.82 | 10.72 |
| query set data/SRR5901135_2.fastq.gz: k-mers in order | 6.02 | 3.25 | 17.85 |
| query set data/SRR5901135_2.fastq.gz: k-mers random order | 4.04 | 2.48 | 17.51 |
| query set data/SRR5901135_2.fastq.gz: sequence, streaming | 19.26 | 25.72 | 10.70 |
| Random negatives, M k-mers/s | 3.89 | 2.36 | 16.91 |

Streaming row, HashSet column: independent lookups parsed from the same sequence.

## eager/summary.txt

```
ecoli       deferred rep 1  insert 0.423 s 10.97 M/s  1.756 B/k-mer  pieces 1643  stored 4557029  redundant 2760
ecoli       eager    rep 1  insert 0.381 s 12.17 M/s  1.756 B/k-mer  pieces 1724  stored 4554269  redundant 0
ecoli       deferred rep 2  insert 0.425 s 10.93 M/s  1.756 B/k-mer  pieces 1643  stored 4557029  redundant 2760
ecoli       eager    rep 2  insert 0.381 s 12.17 M/s  1.756 B/k-mer  pieces 1724  stored 4554269  redundant 0
ecoli       deferred rep 3  insert 0.423 s 10.96 M/s  1.756 B/k-mer  pieces 1643  stored 4557029  redundant 2760
ecoli       eager    rep 3  insert 0.385 s 12.05 M/s  1.756 B/k-mer  pieces 1724  stored 4554269  redundant 0
chr21       deferred rep 1  insert 5.500 s 7.29 M/s  1.929 B/k-mer  pieces 119490  stored 34706358  redundant 81008
chr21       eager    rep 1  insert 5.024 s 7.98 M/s  1.927 B/k-mer  pieces 122335  stored 34625350  redundant 0
chr21       deferred rep 2  insert 5.522 s 7.26 M/s  1.929 B/k-mer  pieces 119490  stored 34706358  redundant 81008
chr21       eager    rep 2  insert 5.078 s 7.89 M/s  1.927 B/k-mer  pieces 122335  stored 34625350  redundant 0
chr21       deferred rep 3  insert 5.491 s 7.30 M/s  1.929 B/k-mer  pieces 119490  stored 34706358  redundant 81008
chr21       eager    rep 3  insert 5.072 s 7.90 M/s  1.927 B/k-mer  pieces 122335  stored 34625350  redundant 0
SRR5833294  deferred rep 1  insert 17.410 s 5.28 M/s  2.505 B/k-mer  pieces 2069193  stored 80215174  redundant 5810
SRR5833294  eager    rep 1  insert 15.839 s 5.80 M/s  2.506 B/k-mer  pieces 2069503  stored 80209364  redundant 0
SRR5833294  deferred rep 2  insert 17.407 s 5.28 M/s  2.505 B/k-mer  pieces 2069193  stored 80215174  redundant 5810
SRR5833294  eager    rep 2  insert 15.796 s 5.82 M/s  2.506 B/k-mer  pieces 2069503  stored 80209364  redundant 0
SRR5833294  deferred rep 3  insert 17.361 s 5.30 M/s  2.505 B/k-mer  pieces 2069193  stored 80215174  redundant 5810
SRR5833294  eager    rep 3  insert 15.769 s 5.83 M/s  2.506 B/k-mer  pieces 2069503  stored 80209364  redundant 0
```

## ablation/summary.txt

```
chr21 T=0  insert 6.821 s 5.88 M/s  1.877 B/k-mer  pieces 122335  entries 4436747  hot 0/0  in-order 2.76  random 1.73  streaming 15.05  negatives 4.08  candidates 11.976
chr21 T=16  insert 5.013 s 8.00 M/s  2.061 B/k-mer  pieces 122335  entries 4214995  hot 4564/2055  in-order 4.97  random 2.77  streaming 33.80  negatives 4.06  candidates 1.839
chr21 T=64  insert 5.057 s 7.93 M/s  1.927 B/k-mer  pieces 122335  entries 4323807  hot 691/97  in-order 4.71  random 2.62  streaming 31.59  negatives 4.05  candidates 2.969
chr21 T=256  insert 5.602 s 7.16 M/s  1.894 B/k-mer  pieces 122335  entries 4388212  hot 92/1  in-order 4.20  random 2.38  streaming 27.16  negatives 4.05  candidates 4.958
ecoli T=0  insert 0.355 s 13.09 M/s  1.754 B/k-mer  pieces 1724  entries 570887  hot 0/0  in-order 6.02  random 4.94  streaming 40.83  negatives 5.88  candidates 1.549
ecoli T=64  insert 0.382 s 12.14 M/s  1.756 B/k-mer  pieces 1724  entries 570803  hot 1/0  in-order 5.95  random 5.00  streaming 39.48  negatives 5.89  candidates 1.545
reads scan=ties  insert 2.204 s 17.24 M/s  2.600 B/k-mer  pieces 251124  entries 1237363  hot 0/0  in-order 5.21  random 3.30  streaming 28.33  negatives 4.80  candidates 1.838
reads scan=range  insert 2.479 s 15.33 M/s  2.600 B/k-mer  pieces 251124  entries 1237363  hot 0/0  in-order 4.56  random 2.99  streaming 27.49  negatives 4.73  candidates 10.848
chr21 max_piece_len=4096  insert 5.060 s 7.92 M/s  1.927 B/k-mer  pieces 122335  entries 4323807  hot 691/97  in-order 4.72  random 2.63  streaming 31.62  negatives 4.03  candidates 2.969
chr21 max_piece_len=1048576  insert 5.055 s 7.93 M/s  1.926 B/k-mer  pieces 120106  entries 4321877  hot 691/97  in-order 4.74  random 2.64  streaming 31.58  negatives 4.04  candidates 2.927
```
