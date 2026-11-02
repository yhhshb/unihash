| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| Memory, B per distinct k-mer | 1.76 | 22.34 | 16.58 |
| Insert, M k-mers/s | 12.00 | 11.80 | 13.39 |
| inserted input (positives): k-mers in order | 5.88 | 3.13 | 13.29 |
| inserted input (positives): k-mers random order | 4.71 | 2.32 | 13.72 |
| inserted input (positives): sequence, streaming | 40.90 | 25.04 | 7.65 |
| Random negatives, M k-mers/s | 5.77 | 2.52 | 29.90 |

Streaming row, HashSet column: independent lookups parsed from the same sequence.
