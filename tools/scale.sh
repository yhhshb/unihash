#!/usr/bin/env bash
# Thread scaling and input-size scaling on the real data sets (see
# tools/fetch_data.sh). Results: bench/results/scale/<data>.txt.
#
# For each data set, unihash-scale builds the single-threaded KmerTable (the
# reference), then a ConcurrentKmerTable with 1, 2, 4, … threads (each
# checked against the reference: same k-mers, none stored twice, exact
# index), measures queries with as many threads, and builds growing
# prefixes of the input.
#
# Environment: THREADS (default: powers of two up to the core count, and
# the core count), SIZES (default 0.125,0.25,0.5,1), REPS (default 3),
# MAX_READS (records of the read set, default 2000000), CAPACITY=1 to size
# the concurrent tables with with_capacity (results in
# bench/results/scale/<data>_capacity.txt), DELETE (fraction of the work
# items deleted concurrently after each build, default 0.25; 0 = off), EXTRA
# (more unihash-scale options, e.g. "--shards 8192").
set -euo pipefail
cd "$(dirname "$0")/.."
cargo build --release --quiet
S=./target/release/unihash-scale
REPS=${REPS:-3}
SIZES=${SIZES:-0.125,0.25,0.5,1}
MAX_READS=${MAX_READS:-2000000}
out=bench/results/scale
mkdir -p "$out"
thr=()
[[ -n "${THREADS:-}" ]] && thr=(--threads "$THREADS")
suffix=""
if [[ "${CAPACITY:-0}" == "1" ]]; then
    thr+=(--capacity)
    suffix="_capacity"
fi
thr+=(--delete-frac "${DELETE:-0.25}")
read -r -a extra <<< "${EXTRA:-}"
thr+=("${extra[@]}")

run() { # name args...
    local name=$1; shift
    echo "== $name: $*"
    "$S" "$@" "${thr[@]}" --reps "$REPS" --sizes "$SIZES" --check | tee "$out/$name$suffix.txt"
}

[[ -s data/ecoli_NC_000913.3.fa ]] && run ecoli data/ecoli_NC_000913.3.fa
[[ -s data/chr21.fa.gz ]] && run chr21 data/chr21.fa.gz
hi=$(ls data/SRR5833294*.fastq.gz 2>/dev/null | head -1 || true)
[[ -n "$hi" ]] && run SRR5833294 "$hi" --max-records "$MAX_READS"
exit 0
