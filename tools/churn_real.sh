#!/usr/bin/env bash
# Deletion / re-insertion churn with several defragmentation policies on the
# real data sets (see tools/fetch_data.sh). Results: bench/results/churn/.
#
# Each run builds the table, defragments it once, then performs OPS
# operations "delete one unit, re-insert the unit deleted LAG operations
# earlier" under each policy:
#   never
#   overhead:X        full defrag() when fragmentation().overhead(k) > X
#   inc:X+full:Y      defrag_incremental() when incremental_overhead() > X,
#                     full defrag() when piece_overhead(k) > Y
#
# Environment: MAX_READS (records of the read set, default 2000000),
# CHECK=1 to verify every policy against the HashSet oracle (slower).
set -euo pipefail
cd "$(dirname "$0")/.."
cargo build --release --quiet
C=./target/release/unihash-churn
MAX_READS=${MAX_READS:-2000000}
mkdir -p bench/results/churn
P=(--policy never --policy overhead:0.1 --policy inc:0.02+full:0.05 --policy inc:0.02+full:0.1 --policy inc:0.05+full:0.2)
chk=()
[[ "${CHECK:-0}" == "1" ]] && chk=(--check)

run() { # name args...
    local name=$1; shift
    echo "== $name: $*"
    "$C" "$@" "${P[@]}" "${chk[@]}" | tee "bench/results/churn/$name.txt"
}

if [[ -s data/ecoli_NC_000913.3.fa ]]; then
    run ecoli data/ecoli_NC_000913.3.fa --unit 300 --ops 40000 --lag 1000
fi
if [[ -s data/chr21.fa.gz ]]; then
    run chr21 data/chr21.fa.gz --unit 1000 --ops 30000 --lag 2000
fi
hi=$(ls data/SRR5833294*.fastq.gz 2>/dev/null | head -1 || true)
if [[ -n "$hi" ]]; then
    run SRR5833294 "$hi" --max-records "$MAX_READS" --ops 500000 --lag 50000
fi
