Insertion, M k-mers/s (seconds); every structure holds the same 80209364 distinct k-mers.

| threads | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| 1, single-threaded structure | 4.11 (22.38 s) | 5.47 (16.82 s) | 4.83 (19.04 s) |
| 1, concurrent structure | 5.68 (16.17 s) | 8.95 (10.27 s) | 5.95 (15.46 s) |
| 2, concurrent structure | 10.62 (8.65 s) | 16.73 (5.50 s) | 11.24 (8.18 s) |
| 4, concurrent structure | 18.70 (4.92 s) | 27.28 (3.37 s) | 17.92 (5.13 s) |
| 8, concurrent structure | 22.31 (4.12 s) | 31.28 (2.94 s) | 19.08 (4.82 s) |

Memory, B per distinct k-mer (heap after the build):

| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| single-threaded structure | 2.506 | 20.289 | 15.060 |
| concurrent structure (8 threads) | 2.799 | 20.289 | 15.061 |
| unihash concurrent table after `into_table` | 2.428 | | |

Queries, positive k-mers (shuffled), M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 1.48 † | 1.13 † | 1.41 † | 5.60 † |
| 2 | 3.15 † | 2.29 † | 2.83 | 11.60 † |
| 4 | 6.06 | 4.43 | 5.04 | 18.78 † |
| 8 | 7.80 | 5.73 | 6.70 | 19.70 |

Queries, uniform random k-mers, M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 3.05 † | 2.32 † | 1.69 † | 10.94 † |
| 2 | 6.33 † | 4.69 † | 3.32 | 22.39 † |
| 4 | 12.54 | 8.54 | 5.76 | 40.66 † |
| 8 | 13.83 | 10.65 | 7.40 | 41.89 |

Queries, streaming (sequences), M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 19.56 † | 13.79 † | 10.06 † | 5.68 † |
| 2 | 37.34 † | 23.51 † | 18.38 | 10.78 † |
| 4 | 71.72 | 37.13 | 35.17 | 19.94 † |
| 8 | 71.18 | 47.27 | 40.98 | 20.36 |

† measured while the one-thread CPU speed was 54–82 % of the run's median (CPU check): the machine was throttled or busy; repeat before drawing conclusions from these cells.

Single-threaded structures grow from empty (kache-hash from its initial capacity); the concurrent ones are created for the distinct k-mer count (unihash `with_capacity`, kache-hash's `max_sz`, `HashSet::with_capacity` per shard). KmerTable, kache-hash and HashSet are queried after the single-threaded build, shared read-only (no locks); ConcurrentKmerTable takes its shard read locks. The concurrent HashSet is 1024 shards of Mutex<HashSet<u64>> filled through per-thread buffers. HashSet streaming: independent lookups of the k-mers parsed from the same sequences.
