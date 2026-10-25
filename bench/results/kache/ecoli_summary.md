| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| Memory, B per distinct k-mer | 1.82 | 22.34 | 16.58 |
| Insert, M k-mers/s | 9.28 | 7.38 | 8.37 |
| inserted input (positives): k-mers in order | 6.46 | 2.66 | 8.03 |
| inserted input (positives): k-mers random order | 2.85 | 1.69 | 8.36 |
| inserted input (positives): sequence, streaming | 43.71 | 16.49 | 4.33 |
| Random negatives, M k-mers/s | 4.14 | 2.04 | 15.07 |

Streaming row, HashSet column: independent lookups parsed from the same sequence.
