#!/usr/bin/env bash
# Deferred vs eager run indexing on the real data sets (see
# tools/fetch_data.sh), REPS runs each, alternating the two modes.
# Full outputs: bench/results/eager/<data>_<mode>_<rep>.txt; summary on stdout.
#
# Environment: REPS (default 3), MAX_READS (records of the read set,
# default 2000000), RESULTS_DIR (default bench/results; outputs go to its
# eager/ subfolder).
set -euo pipefail
cd "$(dirname "$0")/.."
cargo build --release --quiet
B=./target/release/unihash-bench
REPS=${REPS:-3}
MAX_READS=${MAX_READS:-2000000}
out=${RESULTS_DIR:-bench/results}/eager
mkdir -p "$out"

sets=()
[[ -s data/ecoli_NC_000913.3.fa ]] && sets+=("ecoli:data/ecoli_NC_000913.3.fa")
[[ -s data/chr21.fa.gz ]] && sets+=("chr21:data/chr21.fa.gz")
hi=$(ls data/SRR5833294*.fastq.gz 2>/dev/null | head -1 || true)
[[ -n "$hi" ]] && sets+=("SRR5833294:$hi")

for s in "${sets[@]}"; do
    name=${s%%:*}
    file=${s#*:}
    for rep in $(seq 1 "$REPS"); do
        for mode in deferred eager; do
            f="$out/${name}_${mode}_${rep}.txt"
            "$B" "$file" --max-records "$MAX_READS" --no-oracle --query-kmers 5000000 \
                --run-indexing "$mode" >"$f"
            printf '%-11s %-8s rep %s  insert %s  %s  pieces %s  stored %s  redundant %s\n' \
                "$name" "$mode" "$rep" \
                "$(grep -m1 -E '^table ' "$f" | awk '{print $2" s", $4" M/s"}')" \
                "$(grep -m1 'table total, accounted' "$f" | awk '{print $6" B/k-mer"}')" \
                "$(grep -m1 -E '^  pieces ' "$f" | awk '{print $2}')" \
                "$(grep -m1 'stored k-mers' "$f" | awk '{print $3}')" \
                "$(grep -m1 'redundant stored k-mers' "$f" | awk '{print $4}')"
        done
    done
done
