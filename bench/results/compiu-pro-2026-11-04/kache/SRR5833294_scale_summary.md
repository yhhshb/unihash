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
