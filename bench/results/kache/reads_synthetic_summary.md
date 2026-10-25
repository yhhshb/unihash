| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| Memory, B per distinct k-mer | 3.01 | 19.82 | 14.71 |
| Insert, M k-mers/s | 9.43 | 11.34 | 11.37 |
| inserted input (positives): k-mers in order | 5.39 | 2.92 | 8.74 |
| inserted input (positives): k-mers random order | 2.89 | 1.67 | 7.33 |
| inserted input (positives): sequence, streaming | 25.41 | 14.46 | 7.27 |
| held-out records 400001..800000 of each input file: k-mers in order | 5.22 | 3.10 | 9.52 |
| held-out records 400001..800000 of each input file: k-mers random order | 3.11 | 1.73 | 9.16 |
| held-out records 400001..800000 of each input file: sequence, streaming | 24.48 | 16.03 | 7.46 |
| query set synthetic/reads_lo_1.fastq.gz: k-mers in order | 5.64 | 2.94 | 15.09 |
| query set synthetic/reads_lo_1.fastq.gz: k-mers random order | 3.87 | 2.01 | 14.80 |
| query set synthetic/reads_lo_1.fastq.gz: sequence, streaming | 19.63 | 21.29 | 15.26 |
| Random negatives, M k-mers/s | 4.52 | 2.21 | 17.22 |

Streaming row, HashSet column: independent lookups parsed from the same sequence.
