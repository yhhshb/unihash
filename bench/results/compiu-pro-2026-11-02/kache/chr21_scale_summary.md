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
