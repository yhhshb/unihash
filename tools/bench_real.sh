#!/usr/bin/env bash
# Benchmark matrix on real data (see tools/fetch_data.sh).
# Results are written to bench/results/<name>.txt.
#
# Environment overrides: K (31), L (17), MAX_READS (records per read set,
# default 2000000 — remove the cap with MAX_READS=0, which also disables the
# held-out query set).
set -euo pipefail
cd "$(dirname "$0")/.."
cargo build --release --quiet
B=./target/release/unihash-bench
K=${K:-31}
L=${L:-17}
MAX_READS=${MAX_READS:-2000000}
mkdir -p bench/results

run() { # name args...
    local name=$1; shift
    echo "== $name: $*"
    "$B" --k "$K" --l "$L" "$@" | tee "bench/results/$name.txt" \
        | grep -E "^-- |^(table|oracle)  |^(kmers|seq)-|^workload|answered by extension|table total|table heap|heap peak|^RSS|oracle HashSet|MISMATCH|FALSE|peak RSS"
}

if [[ -s data/ecoli_NC_000913.3.fa ]]; then
    run ecoli data/ecoli_NC_000913.3.fa
fi

if [[ -s data/chr21.fa.gz ]]; then
    run chr21 data/chr21.fa.gz
fi

cap=()
[[ "$MAX_READS" != "0" ]] && cap=(--max-records "$MAX_READS")
# High-hit read set: the first MAX_READS records are inserted and the next
# MAX_READS are queried as a held-out set (--holdout); every file of the
# low-hit set (both mates if paired) is queried too. Each query source is
# sampled uniformly (--query-kmers k-mers per source).
hi=$(ls data/SRR5833294*.fastq.gz 2>/dev/null | head -1 || true)
q=()
for f in data/SRR5901135*.fastq.gz; do
    [[ -s "$f" ]] && q+=(--query "$f")
done
if [[ -n "$hi" ]]; then
    hold=()
    [[ "$MAX_READS" != "0" ]] && hold=(--holdout)
    run SRR5833294 "$hi" "${q[@]}" "${cap[@]}" "${hold[@]}"
fi
