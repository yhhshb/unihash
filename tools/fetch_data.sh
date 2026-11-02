#!/usr/bin/env bash
# Download the benchmark datasets into data/.
# Run this on a machine with access to NCBI and UCSC (and ENA for READS_SOURCE=ena).
#
#   bash tools/fetch_data.sh            # E. coli + chr21 (small, ~15 MB compressed)
#   bash tools/fetch_data.sh --reads    # also SRR5833294 and SRR5901135
#
# Read sets. The benchmarks insert the first 2 M records of SRR5833294 and
# query the next 2 M (held out), so only the first READS_RECORDS records
# (default 4000000; 0 = the whole run) are downloaded: the transfer stops
# there instead of fetching the whole run. The record count is written to
# data/SRR5833294.fastq.gz.records; a file with fewer records than requested
# is fetched again. SRR5901135 (paired, low-hit query sets) is fetched whole.
#
# READS_SOURCE selects where the reads come from:
#   sra (default)  NCBI SRA with fastq-dump from the SRA Toolkit, installed
#                  into external/sratoolkit (version SRA_TOOLKIT_VERSION,
#                  default 3.1.1; falls back to NCBI's "current" release).
#                  fastq-dump -X stops after the requested spots.
#   ena            the gzipped FASTQ files of ENA over HTTPS, decompressed
#                  on the fly and cut after the requested records (ENA's
#                  servers can be slow: tens of KB/s at times).
# Both give the same reads in the same order (the headers differ, which the
# benchmarks ignore).
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT=$(pwd)
READS_RECORDS=${READS_RECORDS:-4000000}
READS_SOURCE=${READS_SOURCE:-sra}
mkdir -p data
cd data

fetch() { # url out
    if [[ -s "$2" ]]; then echo "have $2"; return; fi
    echo "fetching $2"
    curl -fL --retry 3 -o "$2.part" "$1" && mv "$2.part" "$2"
}

# Records of a FASTQ file written by this script (none: unknown, assumed whole run).
records_of() { cat "$1.records" 2>/dev/null || echo whole; }

# Is $1 present with at least READS_RECORDS records (or the whole run)?
have_reads() {
    [[ -s "$1" ]] || return 1
    local have
    have=$(records_of "$1")
    [[ "$have" == "whole" ]] && return 0
    [[ "$READS_RECORDS" != "0" && "$have" -ge "$READS_RECORDS" ]]
}

# Record count of a gzipped 4-line FASTQ file.
count_records() { echo $(( $(gzip -dc "$1" | wc -l) / 4 )); }

# ------------------------------------------------------------ SRA Toolkit
SRA_DIR=$ROOT/external/sratoolkit
sra_setup() {
    if [[ -x "$SRA_DIR/bin/fastq-dump" ]]; then return; fi
    local v=${SRA_TOOLKIT_VERSION:-3.1.1}
    local os
    case "$(uname -s)-$(uname -m)" in
        Linux-x86_64)  os=ubuntu64 ;;
        Linux-aarch64) os=alma_linux64 ;;
        Darwin-*)      os=mac-$(uname -m) ;;
        *) echo "no SRA Toolkit build for $(uname -sm); use READS_SOURCE=ena" >&2; exit 1 ;;
    esac
    local tmp
    tmp=$(mktemp -d)
    echo "installing the SRA Toolkit $v into external/sratoolkit"
    if ! curl -fsSL --retry 3 -o "$tmp/sra.tar.gz" \
            "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/$v/sratoolkit.$v-$os.tar.gz"; then
        echo "version $v not found; using NCBI's current release"
        curl -fsSL --retry 3 -o "$tmp/sra.tar.gz" \
            "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-$os.tar.gz"
    fi
    tar -xzf "$tmp/sra.tar.gz" -C "$tmp"
    rm -rf "$SRA_DIR"
    mkdir -p "$ROOT/external"
    mv "$tmp"/sratoolkit.*/ "$SRA_DIR"
    rm -rf "$tmp"
    sra_env
    "$SRA_DIR/bin/fastq-dump" --version | grep -m1 . || true
}
# Private configuration: settings and cache stay inside external/ (never ~/.ncbi).
sra_env() {
    export NCBI_SETTINGS="$SRA_DIR/user-settings.mkfg"
    if [[ ! -s "$NCBI_SETTINGS" ]]; then
        mkdir -p "$SRA_DIR/cache"
        printf '/LIBS/GUID = "%s"\n/repository/user/main/public/root = "%s"\n/repository/user/cache-disabled = "true"\n' \
            "$(cat /proc/sys/kernel/random/uuid 2>/dev/null || echo 00000000-0000-0000-0000-000000000000)" \
            "$SRA_DIR/cache" > "$NCBI_SETTINGS"
    fi
}

# fastq-dump of an accession; first N spots (0 = all); --split-files for paired runs.
sra_dump() { # acc spots split(0/1)
    sra_setup
    sra_env
    local acc=$1 n=$2 split=$3
    local lim=()
    [[ "$n" != "0" ]] && lim=(-X "$n")
    if [[ "$split" == "1" ]]; then
        echo "fetching $acc from NCBI SRA (both mates)"
        rm -f "${acc}"_[12].fastq.gz
        "$SRA_DIR/bin/fastq-dump" --split-files --skip-technical --gzip "${lim[@]}" -O . "$acc"
    else
        echo "fetching ${n/#0/all} records of $acc from NCBI SRA"
        "$SRA_DIR/bin/fastq-dump" --stdout --skip-technical "${lim[@]}" "$acc" | gzip -c > "$acc.fastq.gz.part"
        mv "$acc.fastq.gz.part" "$acc.fastq.gz"
    fi
}

# ------------------------------------------------------------ ENA
ena_urls() { # acc → FASTQ URLs (one per line)
    curl -fsS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$1&result=read_run&fields=fastq_ftp" \
        | tail -n +2 | cut -f2 | tr ';' '\n' | grep .
}
# Stream an ENA gz file, keep the first N records (0 = whole file).
ena_head() { # url out records
    local url=https://$1 out=$2 n=$3
    if [[ "$n" == "0" ]]; then
        rm -f "$out"
        fetch "$url" "$out"
        return
    fi
    echo "fetching the first $n records of $out from ENA"
    set +o pipefail
    curl -fsSL --retry 3 "$url" | gzip -dc 2>/dev/null | head -n $((4 * n)) | gzip -c > "$out.part"
    local st=("${PIPESTATUS[@]}")
    set -o pipefail
    # curl exits with 23 (write error) or 141 (SIGPIPE) when head stops
    # reading: expected.
    if [[ "${st[0]}" != "0" && "${st[0]}" != "23" && "${st[0]}" != "141" ]]; then
        echo "download of $url failed (curl exit ${st[0]})" >&2
        exit 1
    fi
    mv "$out.part" "$out"
}

# ------------------------------------------------------------ data sets

# E. coli K-12 MG1655 (NC_000913.3, ~4.6 Mb)
fetch "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=fasta&retmode=text" \
      ecoli_NC_000913.3.fa

# Human chromosome 21, GRCh38 (soft-masked; N runs are split out by the reader)
fetch "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz" chr21.fa.gz

if [[ "${1:-}" == "--reads" ]]; then
    case "$READS_SOURCE" in sra|ena) ;; *) echo "READS_SOURCE must be sra or ena" >&2; exit 1 ;; esac

    # High-hit set SRR5833294 (single-end): first READS_RECORDS records.
    hi=SRR5833294.fastq.gz
    if have_reads "$hi"; then
        r=$(records_of "$hi")
        if [[ "$r" == "whole" ]]; then echo "have $hi (whole run)"; else echo "have $hi ($r records)"; fi
    else
        rm -f "$hi" "$hi.records" "$hi.part"
        if [[ "$READS_SOURCE" == "sra" ]]; then
            sra_dump SRR5833294 "$READS_RECORDS" 0
        else
            urls=$(ena_urls SRR5833294)
            [[ $(wc -l <<< "$urls") == 1 ]] || { echo "SRR5833294: expected one FASTQ file on ENA, got: $urls" >&2; exit 1; }
            ena_head "$urls" "$hi" "$READS_RECORDS"
        fi
        got=$(count_records "$hi")
        if [[ "$READS_RECORDS" == "0" || "$got" -lt "$READS_RECORDS" ]]; then
            # The whole run (shorter than requested): never fetched again.
            [[ "$READS_RECORDS" == "0" ]] || echo "note: SRR5833294 has only $got records (asked for $READS_RECORDS)"
            echo whole > "$hi.records"
        else
            echo "$got" > "$hi.records"
        fi
        echo "$hi: $got records"
    fi

    # Low-hit query sets SRR5901135 (paired): both mates, whole run.
    if compgen -G "SRR5901135_1.fastq.gz" >/dev/null && compgen -G "SRR5901135_2.fastq.gz" >/dev/null; then
        echo "have SRR5901135_{1,2}.fastq.gz"
    elif [[ "$READS_SOURCE" == "sra" ]]; then
        sra_dump SRR5901135 0 1
    else
        for u in $(ena_urls SRR5901135); do
            fetch "https://${u}" "$(basename "$u")"
        done
    fi
fi
ls -lh
