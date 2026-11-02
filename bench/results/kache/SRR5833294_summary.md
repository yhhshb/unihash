| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| Memory, B per distinct k-mer | 2.51 | 20.29 | 15.06 |
| Insert, M k-mers/s | 4.23 | 5.57 | 5.41 |
| inserted input (positives): k-mers in order | 4.64 | 2.54 | 5.92 |
| inserted input (positives): k-mers random order | 1.50 | 1.43 | 5.91 |
| inserted input (positives): sequence, streaming | 20.09 | 10.29 | 2.97 |
| held-out records 2000001..4000000 of each input file: k-mers in order | 4.59 | 2.64 | 9.45 |
| held-out records 2000001..4000000 of each input file: k-mers random order | 1.92 | 1.71 | 8.47 |
| held-out records 2000001..4000000 of each input file: sequence, streaming | 10.20 | 11.64 | 5.03 |
| query set data/SRR5901135_1.fastq.gz: k-mers in order | 6.41 | 2.84 | 11.81 |
| query set data/SRR5901135_1.fastq.gz: k-mers random order | 3.25 | 1.87 | 11.69 |
| query set data/SRR5901135_1.fastq.gz: sequence, streaming | 17.46 | 17.21 | 6.68 |
| query set data/SRR5901135_2.fastq.gz: k-mers in order | 6.31 | 2.84 | 11.34 |
| query set data/SRR5901135_2.fastq.gz: k-mers random order | 3.04 | 1.86 | 11.62 |
| query set data/SRR5901135_2.fastq.gz: sequence, streaming | 16.80 | 17.08 | 6.62 |
| Random negatives, M k-mers/s | 3.07 | 1.79 | 11.59 |

Streaming row, HashSet column: independent lookups parsed from the same sequence.
