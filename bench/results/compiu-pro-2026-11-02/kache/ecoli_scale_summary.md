Insertion, M k-mers/s (seconds); every structure holds the same 4554269 distinct k-mers.

| threads | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| 1, single-threaded structure | 12.23 (0.38 s) | 11.61 (0.40 s) | 14.08 (0.33 s) |
| 1, concurrent structure | 13.14 (0.35 s) | 15.11 (0.31 s) | 23.30 (0.20 s) |
| 2, concurrent structure | 25.34 (0.18 s) | 29.59 (0.16 s) | 41.30 (0.11 s) |
| 4, concurrent structure | 49.42 (0.09 s) | 57.73 (0.08 s) | 60.65 (0.08 s) |
| 8, concurrent structure | 73.07 (0.06 s) | 94.25 (0.05 s) | 70.65 (0.07 s) |
| 12, concurrent structure | 75.72 (0.06 s) | 99.04 (0.05 s) | 72.52 (0.06 s) |

Memory, B per distinct k-mer (heap after the build):

| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| single-threaded structure | 1.754 | 22.335 | 16.577 |
| concurrent structure (12 threads) | 2.170 | 22.335 | 16.593 |
| unihash concurrent table after `into_table` | 1.695 | | |

Queries, positive k-mers (shuffled), M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 5.09 | 4.29 | 2.37 | 14.79 |
| 2 | 10.13 | 8.20 | 4.74 | 29.94 |
| 4 | 20.30 | 15.97 | 9.52 | 59.17 |
| 8 | 26.71 | 21.96 | 14.70 | 62.48 |
| 12 | 37.73 | 30.52 | 20.95 | 89.50 |

Queries, uniform random k-mers, M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 5.92 | 5.36 | 2.56 | 41.15 |
| 2 | 11.86 | 10.25 | 5.10 | 82.59 |
| 4 | 23.69 | 19.82 | 10.18 | 166.64 |
| 8 | 31.73 | 27.29 | 15.68 | 202.58 |
| 12 | 42.00 | 36.35 | 21.85 | 299.76 |

Queries, streaming (sequences), M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 44.94 | 32.37 | 25.08 | 14.93 |
| 2 | 88.49 | 63.76 | 48.48 | 29.67 |
| 4 | 176.29 | 126.56 | 94.66 | 58.87 |
| 8 | 220.28 | 152.04 | 134.59 | 71.65 |
| 12 | 327.57 | 223.26 | 190.06 | 89.76 |

Single-threaded structures grow from empty (kache-hash from its initial capacity); the concurrent ones are created for the distinct k-mer count (unihash `with_capacity`, kache-hash's `max_sz`, `HashSet::with_capacity` per shard). KmerTable, kache-hash and HashSet are queried after the single-threaded build, shared read-only (no locks); ConcurrentKmerTable takes its shard read locks. The concurrent HashSet is 1024 shards of Mutex<HashSet<u64>> filled through per-thread buffers. HashSet streaming: independent lookups of the k-mers parsed from the same sequences.
