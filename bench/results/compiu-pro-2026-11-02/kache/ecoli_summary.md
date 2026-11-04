| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| Memory, B per distinct k-mer | 1.76 | 22.34 | 16.58 |
| Insert, M k-mers/s | 9.57 | 10.92 | 11.90 |
| inserted input (positives): k-mers in order | 5.58 | 3.06 | 12.05 |
| inserted input (positives): k-mers random order | 4.03 | 2.14 | 11.97 |
| inserted input (positives): sequence, streaming | 40.16 | 24.08 | 6.46 |
| Random negatives, M k-mers/s | 4.93 | 2.41 | 25.92 |

Streaming row, HashSet column: independent lookups parsed from the same sequence.
