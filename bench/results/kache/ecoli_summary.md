| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| Memory, B per distinct k-mer | 1.76 | 22.34 | 16.58 |
| Insert, M k-mers/s | 8.07 | 6.91 | 7.47 |
| inserted input (positives): k-mers in order | 6.03 | 2.71 | 8.05 |
| inserted input (positives): k-mers random order | 3.37 | 1.67 | 6.66 |
| inserted input (positives): sequence, streaming | 44.76 | 14.37 | 4.19 |
| Random negatives, M k-mers/s | 3.92 | 1.80 | 10.17 |

Streaming row, HashSet column: independent lookups parsed from the same sequence.
