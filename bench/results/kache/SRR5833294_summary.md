| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| Memory, B per distinct k-mer | 3.05 | 20.29 | 15.06 |
| Insert, M k-mers/s | 3.73 | 5.73 | 5.47 |
| inserted input (positives): k-mers in order | 4.77 | 2.55 | 5.88 |
| inserted input (positives): k-mers random order | 1.48 | 1.44 | 5.88 |
| inserted input (positives): sequence, streaming | 19.79 | 10.26 | 2.95 |
| held-out records 2000001..4000000 of each input file: k-mers in order | 4.92 | 2.64 | 9.40 |
| held-out records 2000001..4000000 of each input file: k-mers random order | 1.96 | 1.71 | 8.53 |
| held-out records 2000001..4000000 of each input file: sequence, streaming | 10.44 | 11.67 | 5.07 |
| query set data/SRR5901135_1.fastq.gz: k-mers in order | 7.19 | 2.84 | 11.76 |
| query set data/SRR5901135_1.fastq.gz: k-mers random order | 3.43 | 1.87 | 11.63 |
| query set data/SRR5901135_1.fastq.gz: sequence, streaming | 17.79 | 17.30 | 6.66 |
| query set data/SRR5901135_2.fastq.gz: k-mers in order | 7.12 | 2.83 | 11.69 |
| query set data/SRR5901135_2.fastq.gz: k-mers random order | 3.38 | 1.86 | 11.63 |
| query set data/SRR5901135_2.fastq.gz: sequence, streaming | 17.73 | 17.18 | 6.62 |
| Random negatives, M k-mers/s | 3.19 | 1.80 | 11.50 |

Streaming row, HashSet column: independent lookups parsed from the same sequence.
