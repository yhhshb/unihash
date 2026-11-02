Insertion, M k-mers/s (seconds); every structure holds the same 4554269 distinct k-mers.

| threads | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| 1, single-threaded structure | 9.57 (0.48 s) | 7.30 (0.64 s) | 7.91 (0.59 s) |
| 1, concurrent structure | 11.70 (0.40 s) | 13.81 (0.34 s) | 9.06 (0.51 s) |
| 2, concurrent structure | 22.43 (0.21 s) | 25.65 (0.18 s) | 16.88 (0.28 s) |
| 4, concurrent structure | 40.34 (0.12 s) | 44.48 (0.10 s) | 26.54 (0.17 s) |
| 8, concurrent structure | 51.18 (0.09 s) | 54.49 (0.09 s) | 29.66 (0.16 s) |

Memory, B per distinct k-mer (heap after the build):

| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| single-threaded structure | 1.754 | 22.335 | 16.577 |
| concurrent structure (8 threads) | 2.170 | 22.335 | 16.593 |
| unihash concurrent table after `into_table` | 1.695 | | |

Queries, positive k-mers (shuffled), M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 3.75 † | 3.21 † | 1.70 † | 8.34 † |
| 2 | 7.69 † | 6.45 † | 3.47 † | 17.16 † |
| 4 | 14.97 | 12.09 | 6.57 | 30.02 |
| 8 | 19.92 | 17.42 | 9.95 | 32.43 |

Queries, uniform random k-mers, M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 4.16 † | 3.79 † | 2.01 † | 16.38 † |
| 2 | 8.80 † | 7.63 † | 4.04 † | 34.86 † |
| 4 | 17.76 | 14.14 | 7.47 | 68.13 |
| 8 | 22.32 | 19.65 | 11.34 | 75.74 |

Queries, streaming (sequences), M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 42.84 † | 31.98 † | 17.27 † | 8.51 † |
| 2 | 78.94 † | 60.17 † | 33.44 † | 17.42 † |
| 4 | 150.73 | 101.82 | 67.03 | 32.56 |
| 8 | 229.01 | 138.93 | 80.43 | 34.82 |

† measured while the one-thread CPU speed was 56–71 % of the run's median (CPU check): the machine was throttled or busy; repeat before drawing conclusions from these cells.

Single-threaded structures grow from empty (kache-hash from its initial capacity); the concurrent ones are created for the distinct k-mer count (unihash `with_capacity`, kache-hash's `max_sz`, `HashSet::with_capacity` per shard). KmerTable, kache-hash and HashSet are queried after the single-threaded build, shared read-only (no locks); ConcurrentKmerTable takes its shard read locks. The concurrent HashSet is 1024 shards of Mutex<HashSet<u64>> filled through per-thread buffers. HashSet streaming: independent lookups of the k-mers parsed from the same sequences.
