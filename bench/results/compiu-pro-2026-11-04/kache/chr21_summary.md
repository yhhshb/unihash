| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| Memory, B per distinct k-mer | 1.93 | 23.50 | 17.44 |
| Insert, M k-mers/s | 7.80 | 9.17 | 8.61 |
| inserted input (positives): k-mers in order | 4.69 | 3.06 | 11.03 |
| inserted input (positives): k-mers random order | 2.57 | 2.07 | 10.76 |
| inserted input (positives): sequence, streaming | 31.52 | 20.69 | 5.82 |
| Random negatives, M k-mers/s | 4.02 | 2.39 | 17.99 |

Streaming row, HashSet column: independent lookups parsed from the same sequence.
