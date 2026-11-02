Insertion, M k-mers/s (seconds); every structure holds the same 34625350 distinct k-mers.

| threads | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| 1, single-threaded structure | 5.71 (7.01 s) | 5.66 (7.09 s) | 5.56 (7.21 s) |
| 1, concurrent structure | 7.28 (5.51 s) | 10.91 (3.67 s) | 6.87 (5.84 s) |
| 2, concurrent structure | 13.77 (2.91 s) | 19.64 (2.04 s) | 12.67 (3.16 s) |
| 4, concurrent structure | 23.38 (1.71 s) | 30.87 (1.30 s) | 21.17 (1.89 s) |
| 8, concurrent structure | 29.20 (1.37 s) | 35.70 (1.12 s) | 22.68 (1.77 s) |

Memory, B per distinct k-mer (heap after the build):

| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| single-threaded structure | 1.927 | 23.500 | 17.443 |
| concurrent structure (8 threads) | 2.451 | 23.500 | 17.445 |
| unihash concurrent table after `into_table` | 1.867 | | |

Queries, positive k-mers (shuffled), M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 1.81 | 1.41 | 1.50 † | 6.95 † |
| 2 | 3.86 † | 2.22 † | 3.00 † | 14.03 † |
| 4 | 6.92 | 5.58 | 5.37 | 25.92 † |
| 8 | 8.56 † | 7.05 † | 7.25 † | 25.14 |

Queries, uniform random k-mers, M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 3.14 | 2.95 | 1.87 † | 12.46 † |
| 2 | 6.44 † | 5.68 † | 3.54 † | 25.41 † |
| 4 | 11.49 | 11.50 | 6.21 | 48.26 † |
| 8 | 12.98 † | 13.08 † | 8.13 † | 51.22 |

Queries, streaming (sequences), M queries/s:

| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |
|---|---|---|---|---|
| 1 | 37.11 | 25.81 | 15.58 † | 7.96 † |
| 2 | 62.06 † | 35.98 † | 24.54 † | 13.10 † |
| 4 | 87.30 | 71.74 | 46.31 | 27.59 † |
| 8 | 114.03 † | 76.74 † | 56.42 † | 29.88 |

† measured while the one-thread CPU speed was 49–84 % of the run's median (CPU check): the machine was throttled or busy; repeat before drawing conclusions from these cells.

Single-threaded structures grow from empty (kache-hash from its initial capacity); the concurrent ones are created for the distinct k-mer count (unihash `with_capacity`, kache-hash's `max_sz`, `HashSet::with_capacity` per shard). KmerTable, kache-hash and HashSet are queried after the single-threaded build, shared read-only (no locks); ConcurrentKmerTable takes its shard read locks. The concurrent HashSet is 1024 shards of Mutex<HashSet<u64>> filled through per-thread buffers. HashSet streaming: independent lookups of the k-mers parsed from the same sequences.
