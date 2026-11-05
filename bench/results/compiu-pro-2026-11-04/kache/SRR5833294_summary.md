| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| Memory, B per distinct k-mer | 2.51 | 20.29 | 15.06 |
| Insert, M k-mers/s | 5.61 | 8.77 | 7.77 |
| inserted input (positives): k-mers in order | 4.65 | 3.05 | 10.60 |
| inserted input (positives): k-mers random order | 2.08 | 2.11 | 10.60 |
| inserted input (positives): sequence, streaming | 22.05 | 17.05 | 5.57 |
| held-out records 2000001..4000000 of each input file: k-mers in order | 4.70 | 3.13 | 15.04 |
| held-out records 2000001..4000000 of each input file: k-mers random order | 2.59 | 2.31 | 13.76 |
| held-out records 2000001..4000000 of each input file: sequence, streaming | 12.01 | 18.42 | 8.45 |
| query set data/SRR5901135_1.fastq.gz: k-mers in order | 6.03 | 3.26 | 15.43 |
| query set data/SRR5901135_1.fastq.gz: k-mers random order | 4.07 | 2.49 | 17.09 |
| query set data/SRR5901135_1.fastq.gz: sequence, streaming | 19.26 | 25.82 | 10.72 |
| query set data/SRR5901135_2.fastq.gz: k-mers in order | 6.02 | 3.25 | 17.85 |
| query set data/SRR5901135_2.fastq.gz: k-mers random order | 4.04 | 2.48 | 17.51 |
| query set data/SRR5901135_2.fastq.gz: sequence, streaming | 19.26 | 25.72 | 10.70 |
| Random negatives, M k-mers/s | 3.89 | 2.36 | 16.91 |

Streaming row, HashSet column: independent lookups parsed from the same sequence.
