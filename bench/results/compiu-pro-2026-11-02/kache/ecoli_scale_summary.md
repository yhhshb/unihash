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
