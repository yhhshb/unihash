| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| Memory, B per distinct k-mer | 1.93 | 23.50 | 17.44 |
| Insert, M k-mers/s | 7.72 | 9.25 | 8.63 |
| inserted input (positives): k-mers in order | 4.67 | 3.02 | 10.67 |
| inserted input (positives): k-mers random order | 2.45 | 2.14 | 10.22 |
| inserted input (positives): sequence, streaming | 31.28 | 20.65 | 5.84 |
| Random negatives, M k-mers/s | 3.99 | 2.40 | 18.07 |

Streaming row, HashSet column: independent lookups parsed from the same sequence.
