| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| Memory, B per distinct k-mer | 1.76 | 22.34 | 16.58 |
| Insert, M k-mers/s | 12.13 | 11.85 | 13.61 |
| inserted input (positives): k-mers in order | 5.95 | 3.18 | 13.50 |
| inserted input (positives): k-mers random order | 4.90 | 2.34 | 13.39 |
| inserted input (positives): sequence, streaming | 40.78 | 24.97 | 8.02 |
| Random negatives, M k-mers/s | 5.71 | 2.51 | 26.96 |

Streaming row, HashSet column: independent lookups parsed from the same sequence.
