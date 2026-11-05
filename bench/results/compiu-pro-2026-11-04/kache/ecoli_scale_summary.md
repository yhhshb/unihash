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
