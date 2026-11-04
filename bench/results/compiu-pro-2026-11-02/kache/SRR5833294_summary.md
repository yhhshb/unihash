| | unihash | kache-hash | HashSet<u64> |
|---|---|---|---|
| Memory, B per distinct k-mer | 2.51 | 20.29 | 15.06 |
| Insert, M k-mers/s | 5.48 | 7.40 | 7.07 |
| inserted input (positives): k-mers in order | 4.29 | 2.79 | 8.59 |
| inserted input (positives): k-mers random order | 1.83 | 1.76 | 8.65 |
| inserted input (positives): sequence, streaming | 19.93 | 14.88 | 4.46 |
| held-out records 2000001..4000000 of each input file: k-mers in order | 4.41 | 2.97 | 12.68 |
| held-out records 2000001..4000000 of each input file: k-mers random order | 2.33 | 2.03 | 11.44 |
| held-out records 2000001..4000000 of each input file: sequence, streaming | 11.15 | 16.15 | 7.24 |
| query set data/SRR5901135_1.fastq.gz: k-mers in order | 5.84 | 3.04 | 15.98 |
| query set data/SRR5901135_1.fastq.gz: k-mers random order | 3.76 | 2.21 | 15.52 |
| query set data/SRR5901135_1.fastq.gz: sequence, streaming | 18.13 | 22.95 | 9.38 |
| query set data/SRR5901135_2.fastq.gz: k-mers in order | 5.81 | 3.07 | 15.46 |
| query set data/SRR5901135_2.fastq.gz: k-mers random order | 3.67 | 2.10 | 15.55 |
| query set data/SRR5901135_2.fastq.gz: sequence, streaming | 17.77 | 21.95 | 9.32 |
| Random negatives, M k-mers/s | 3.64 | 2.06 | 15.39 |

Streaming row, HashSet column: independent lookups parsed from the same sequence.
