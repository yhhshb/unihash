| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| Memory, B per distinct k-mer | 2.51 | 20.29 | 15.06 |
| Insert, M k-mers/s | 5.81 | 8.85 | 7.73 |
| inserted input (positives): k-mers in order | 4.64 | 3.04 | 10.61 |
| inserted input (positives): k-mers random order | 2.07 | 2.11 | 10.61 |
| inserted input (positives): sequence, streaming | 22.03 | 16.99 | 5.57 |
| held-out records 2000001..4000000 of each input file: k-mers in order | 4.69 | 3.12 | 15.01 |
| held-out records 2000001..4000000 of each input file: k-mers random order | 2.59 | 2.31 | 13.72 |
| held-out records 2000001..4000000 of each input file: sequence, streaming | 12.04 | 18.31 | 8.44 |
| query set data/SRR5901135_1.fastq.gz: k-mers in order | 6.00 | 3.25 | 17.59 |
| query set data/SRR5901135_1.fastq.gz: k-mers random order | 3.97 | 2.48 | 17.41 |
| query set data/SRR5901135_1.fastq.gz: sequence, streaming | 19.00 | 25.68 | 10.31 |
| query set data/SRR5901135_2.fastq.gz: k-mers in order | 5.97 | 3.25 | 17.54 |
| query set data/SRR5901135_2.fastq.gz: k-mers random order | 3.98 | 2.47 | 17.35 |
| query set data/SRR5901135_2.fastq.gz: sequence, streaming | 19.26 | 25.57 | 10.57 |
| Random negatives, M k-mers/s | 3.85 | 2.36 | 16.74 |

Streaming row, HashSet column: independent lookups parsed from the same sequence.
