| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| Memory, B per distinct k-mer | 1.93 | 23.50 | 17.44 |
| Insert, M k-mers/s | 6.25 | 6.06 | 5.96 |
| inserted input (positives): k-mers in order | 4.87 | 2.67 | 6.79 |
| inserted input (positives): k-mers random order | 1.89 | 1.56 | 6.57 |
| inserted input (positives): sequence, streaming | 31.66 | 13.41 | 3.46 |
| Random negatives, M k-mers/s | 3.17 | 1.95 | 12.32 |

Streaming row, HashSet column: independent lookups parsed from the same sequence.
