#!/usr/bin/env bash
# Head-to-head comparison of unihash with kache-hash
# (Khan, Patro & Pandey, "kache-hash: A dynamic, concurrent, and cache-efficient
# hash table for streaming k-mer operations", https://github.com/jamshed/kache-hash).
#
#   bash tools/kache_compare.sh                 # E. coli, chr21, SRR5833294 (as bench_real.sh)
#   bash tools/kache_compare.sh NAME ARGS...    # one single-threaded run: NAME + unihash-bench arguments
#   bash tools/kache_compare.sh --setup         # only fetch and build kache-hash and its drivers
#
# kache-hash is cloned into external/kache-hash (pinned commit) and the drivers
# are built into external/bin/; external/ is not part of the repository.
# tools/kache/kache-hash-lock.patch is applied first: it fixes a race in the
# bucket lock of kache-hash's concurrent table that can store a k-mer twice
# (found by the size check below; KACHE_PATCH=0 builds the upstream code).
#
# Two comparisons per data set.
#
# Single-threaded (all query patterns, held-out reads and extra query sets):
#   1. unihash-bench inserts the input, measures every query pattern, and writes
#      the workloads (window sequences, in-order and shuffled k-mers, negatives)
#      plus expected hit counts to a dump directory (--dump).
#   2. tools/kache/kache_bench builds kache-hash (single-threaded hash set, same
#      k = 31, ℓ = 17) from the same input files, runs the same workloads, and
#      checks its answers against the expected hit counts.
#   3. tools/kache/summarize.py prints a side-by-side table.
# Outputs: $RESULTS_DIR/kache/<NAME>_{unihash,kache,summary}.{txt,md}
#
# Thread scaling (single-threaded and concurrent versions of each structure):
#   1. unihash-scale --hashset builds KmerTable and std HashSet<u64> on one
#      thread, then ConcurrentKmerTable and a sharded Mutex<HashSet<u64>> with
#      each thread count (pre-sized for the distinct k-mers), queries the
#      single-threaded tables and the concurrent unihash table with as many
#      threads, and writes the work items and query workloads (--dump).
#   2. tools/kache/kache_scale does the same with kache-hash: the
#      single-threaded table (growing), queried from each thread count, and
#      the concurrent table (Streaming_Kmer_Hash_Table<k, true>, one token per
#      thread, created for the distinct k-mers) built with each thread count.
#      Every size and answer is checked against unihash's.
#   3. tools/kache/summarize_scale.py prints the tables.
# Outputs: $RESULTS_DIR/kache/<NAME>_scale_{unihash,kache,summary}.{txt,md}
#
# Environment: CXX (g++), KACHE_INIT_CAP (1048576), MAX_READS (2000000),
# DUMP_ROOT (target/kache-dump), KEEP_DUMP=1 to keep the dumps; SINGLE=0 or
# SCALE=0 to skip a comparison; THREADS (default: powers of two up to the
# core count, and the core count) and REPS (median of N builds, default 3)
# for the scaling comparison; DATA (comma-separated subset of
# ecoli,chr21,SRR5833294; default all available); RESULTS_DIR (default
# bench/results).
#
# Build note: kache-hash miscompiles with GCC -O3 unless strict aliasing is
# disabled (its answers become wrong: streaming queries find ~1 % of present
# k-mers). We therefore compile it with -fno-strict-aliasing; every run checks
# kache-hash's answers against unihash's oracle and the script fails on any
# mismatch.
#
# kache-hash's RW_Lock (used by its concurrent table) calls
# std::atomic_flag::test with memory_order_acq_rel, which is not a valid order
# for a load; GCC warns (-Winvalid-memory-model) and uses seq_cst, a stronger
# order, so behaviour is unaffected. The warning is silenced for the build.
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT=$(pwd)

KACHE_COMMIT=f7cfa5be4bcd28ffa061588c7c80445cbd5690d0
KDIR=external/kache-hash
BIN=external/bin
CXX=${CXX:-g++}

setup_kache() {
    if [[ ! -d "$KDIR/.git" ]]; then
        mkdir -p external
        echo "cloning kache-hash into $KDIR"
        git clone --quiet https://github.com/jamshed/kache-hash "$KDIR"
    fi
    if [[ "$(git -C "$KDIR" rev-parse HEAD)" != "$KACHE_COMMIT" ]]; then
        git -C "$KDIR" fetch --quiet origin || true
    fi
    # Pinned commit, then the lock fix (see tools/kache/kache-hash-lock.patch).
    git -C "$KDIR" checkout --quiet --force "$KACHE_COMMIT"
    local patch=none
    if [[ "${KACHE_PATCH:-1}" == "1" ]]; then
        git -C "$KDIR" apply "$ROOT/tools/kache/kache-hash-lock.patch"
        patch=$(cksum < tools/kache/kache-hash-lock.patch | cut -d' ' -f1)
        echo "kache-hash $KACHE_COMMIT with tools/kache/kache-hash-lock.patch"
    else
        echo "kache-hash $KACHE_COMMIT unpatched (KACHE_PATCH=0)"
    fi
    mkdir -p "$BIN"
    local L1
    L1=$(getconf LEVEL1_DCACHE_LINESIZE 2>/dev/null || echo 64)
    [[ "$L1" =~ ^[0-9]+$ && "$L1" -gt 0 ]] || L1=64
    # Rebuild when the driver, the pinned commit or the patch changed
    # (-march=native: binaries copied from another machine are rebuilt too).
    local stamp="$BIN/.kache-build"
    local want host
    host=$(uname -n)
    want="$KACHE_COMMIT patch=$patch $CXX L1=$L1 $host"
    for tool in kache_bench kache_scale kache_lock_stress; do
        if [[ ! "$BIN/$tool" -nt tools/kache/$tool.cpp || "$(cat "$stamp" 2>/dev/null)" != "$want" ]]; then
            echo "building $BIN/$tool with $CXX"
            "$CXX" -std=c++20 -O3 -funroll-loops -march=native -fno-strict-aliasing -Wno-invalid-memory-model \
                -DNDEBUG -DQUIET \
                -DXXH_INLINE_ALL -DXXH_STATIC_LINKING_ONLY -DL1_CACHE_LINE_SIZE="$L1" \
                -I "$KDIR/include" tools/kache/$tool.cpp -o "$BIN/$tool" -lz -lpthread
        fi
    done
    echo "$want" > "$stamp"
}

setup_kache
if [[ "${1:-}" == "--setup" ]]; then
    exit 0
fi

cargo build --release --quiet
RES=${RESULTS_DIR:-bench/results}
mkdir -p target "$RES/kache"

DUMP_ROOT=${DUMP_ROOT:-target/kache-dump}
KACHE_INIT_CAP=${KACHE_INIT_CAP:-1048576}

run() { # name unihash-bench-args...
    local name=$1; shift
    local dump="$DUMP_ROOT/$name"
    local out=$RES/kache/$name
    rm -rf "$dump"
    echo "== $name: unihash-bench $*"
    ./target/release/unihash-bench "$@" --dump "$dump" > "${out}_unihash.txt"
    echo "== $name: kache-hash"
    "$BIN/kache_bench" "$dump" "$KACHE_INIT_CAP" > "${out}_kache.txt" || {
        cat "${out}_kache.txt"; echo "kache-hash run failed or disagreed with the oracle" >&2; exit 1; }
    python3 tools/kache/summarize.py "${out}_unihash.txt" "${out}_kache.txt" | tee "${out}_summary.md"
    [[ "${KEEP_DUMP:-0}" == "1" ]] || rm -rf "$dump"
}

REPS=${REPS:-3}
scale() { # name unihash-scale-input-args...
    local name=$1; shift
    local dump="$DUMP_ROOT/${name}_scale"
    local out=$RES/kache/${name}_scale
    local thr=()
    [[ -n "${THREADS:-}" ]] && thr=(--threads "$THREADS")
    rm -rf "$dump"
    echo "== $name, thread scaling: unihash-scale $*"
    ./target/release/unihash-scale "$@" "${thr[@]}" --capacity --check --hashset --reps "$REPS" \
        --dump "$dump" > "${out}_unihash.txt" || {
        cat "${out}_unihash.txt"; echo "unihash-scale failed" >&2; exit 1; }
    echo "== $name, thread scaling: kache-hash"
    "$BIN/kache_scale" "$dump" "$KACHE_INIT_CAP" > "${out}_kache.txt" || {
        cat "${out}_kache.txt"; echo "kache-hash run failed or disagreed with unihash" >&2; exit 1; }
    python3 tools/kache/summarize_scale.py "${out}_unihash.txt" "${out}_kache.txt" | tee "${out}_summary.md"
    [[ "${KEEP_DUMP:-0}" == "1" ]] || rm -rf "$dump"
}

if [[ $# -ge 1 ]]; then
    run "$@"
    exit 0
fi

MAX_READS=${MAX_READS:-2000000}
hi=$(ls data/SRR5833294*.fastq.gz 2>/dev/null | head -1 || true)
want() { [[ -z "${DATA:-}" || ",$DATA," == *",$1,"* ]]; }
if [[ "${SINGLE:-1}" == "1" ]]; then
    if want ecoli && [[ -s data/ecoli_NC_000913.3.fa ]]; then
        run ecoli data/ecoli_NC_000913.3.fa
    fi
    if want chr21 && [[ -s data/chr21.fa.gz ]]; then
        run chr21 data/chr21.fa.gz
    fi
    if want SRR5833294 && [[ -n "$hi" ]]; then
        q=()
        for f in data/SRR5901135*.fastq.gz; do
            [[ -s "$f" ]] && q+=(--query "$f")
        done
        run SRR5833294 "$hi" "${q[@]}" --max-records "$MAX_READS" --holdout
    fi
fi
if [[ "${SCALE:-1}" == "1" ]]; then
    if want ecoli && [[ -s data/ecoli_NC_000913.3.fa ]]; then
        scale ecoli data/ecoli_NC_000913.3.fa
    fi
    if want chr21 && [[ -s data/chr21.fa.gz ]]; then
        scale chr21 data/chr21.fa.gz
    fi
    if want SRR5833294 && [[ -n "$hi" ]]; then
        scale SRR5833294 "$hi" --max-records "$MAX_READS"
    fi
fi
exit 0
