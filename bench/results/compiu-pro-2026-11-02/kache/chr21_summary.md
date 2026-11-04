| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| Memory, B per distinct k-mer | 1.93 | 23.50 | 17.44 |
| Insert, M k-mers/s | 7.10 | 9.03 | 7.86 |
| inserted input (positives): k-mers in order | 4.49 | 2.99 | 9.79 |
| inserted input (positives): k-mers random order | 2.28 | 2.04 | 9.77 |
| inserted input (positives): sequence, streaming | 29.10 | 19.11 | 5.09 |
| Random negatives, M k-mers/s | 3.72 | 2.31 | 16.32 |

Streaming row, HashSet column: independent lookups parsed from the same sequence.
