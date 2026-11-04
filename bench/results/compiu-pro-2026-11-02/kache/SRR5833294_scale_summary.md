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
