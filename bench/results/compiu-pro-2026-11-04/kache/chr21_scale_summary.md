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
