| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| Memory, B per distinct k-mer | 2.02 | 23.50 | 17.44 |
| Insert, M k-mers/s | 5.46 | 5.89 | 5.88 |
| inserted input (positives): k-mers in order | 4.39 | 2.50 | 6.42 |
| inserted input (positives): k-mers random order | 1.53 | 1.53 | 5.84 |
| inserted input (positives): sequence, streaming | 29.93 | 13.42 | 3.28 |
| Random negatives, M k-mers/s | 3.17 | 1.92 | 11.97 |

Streaming row, HashSet column: independent lookups parsed from the same sequence.
