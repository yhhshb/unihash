#!/usr/bin/env bash
# Download the Milestone 6 benchmark datasets into data/.
# Run this on a machine with access to NCBI, UCSC and ENA.
#
#   tools/fetch_data.sh            # E. coli + chr21 (small, ~15 MB compressed)
#   tools/fetch_data.sh --reads    # also SRR5833294 and SRR5901135 (large)
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data
cd data

fetch() { # url out
    if [[ -s "$2" ]]; then echo "have $2"; return; fi
    echo "fetching $2"
    curl -fL --retry 3 -o "$2.part" "$1" && mv "$2.part" "$2"
}

# E. coli K-12 MG1655 (NC_000913.3, ~4.6 Mb)
fetch "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=fasta&retmode=text" \
      ecoli_NC_000913.3.fa

# Human chromosome 21, GRCh38 (soft-masked; N runs are split out by the reader)
fetch "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz" chr21.fa.gz

if [[ "${1:-}" == "--reads" ]]; then
    for acc in SRR5833294 SRR5901135; do
        urls=$(curl -fsS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${acc}&result=read_run&fields=fastq_ftp" \
               | tail -n +2 | cut -f2 | tr ';' '\n')
        [[ -n "$urls" ]] || { echo "no FASTQ URLs found for $acc" >&2; exit 1; }
        for u in $urls; do
            fetch "https://${u}" "$(basename "$u")"
        done
    done
fi
ls -lh
