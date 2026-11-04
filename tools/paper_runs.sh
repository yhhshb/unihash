#!/usr/bin/env bash
# Every measurement reported in the unihash paper, on the current machine.
#
#   bash tools/paper_runs.sh               # everything (several hours with the read sets)
#   STEPS=kache,scale bash tools/paper_runs.sh
#
# The script
#   1. checks the build tools (cargo, a C++20 compiler, zlib, git, curl,
#      python3); if cargo is missing, it installs a private Rust toolchain
#      into external/rust (nothing outside the repository is touched);
#   2. downloads the data sets into data/ with tools/fetch_data.sh --reads
#      when any of them is missing (only the 2 x MAX_READS records of
#      SRR5833294 that the runs use);
#   3. clones kache-hash (pinned commit) into external/kache-hash and builds
#      its drivers into external/bin (tools/kache_compare.sh --setup);
#   4. records the machine (CPU, memory, OS, compilers, governor, commit);
#   5. runs the steps below, writing everything to one results folder.
# data/ and external/ are ignored by git: only the scripts that rebuild them
# are part of the repository.
#
# Steps (STEPS, comma-separated, default all, in this order):
#   tests     cargo test --release (the build is correct on this machine)
#   kache     tools/kache_compare.sh: single-threaded comparison with
#             kache-hash and HashSet (memory, throughput, per-query
#             counters) and thread scaling of the concurrent structures
#             (insertion and queries)
#   scale     tools/scale.sh CAPACITY=1 DELETE=0.25 with --defrag-compare:
#             ConcurrentKmerTable scaling for insertion, queries and
#             deletion, with checks; the same table queried with shard read
#             locks and through a read snapshot (no lock per query: the cost
#             of the read locks); defrag / defrag_incremental after deletion
#             in place vs a round trip through KmerTable
#   eager     tools/eager_compare.sh: deferred vs eager run indexing
#   ablation  hot-minimizer threshold sweep (chr21, E. coli), candidate scan
#             ties vs range (synthetic reads), max_piece_len 4096 vs 1 Mb
#             (chr21)
#   defrag    unihash-bench --defrag after bulk loading (SRR5833294, chr21)
#   churn     tools/churn_real.sh CHECK=1: deletion / re-insertion churn
#             under each defragmentation policy, checked against the oracle
#
# Environment:
#   OUT        results folder (default bench/results/<host>-<YYYY-MM-DD>)
#   STEPS      subset of the steps above
#   FORCE=1    re-run steps already completed in OUT (default: resume,
#              skipping steps marked done)
#   THREADS    thread counts for the scaling runs (default: every count
#              from 1 to the hardware threads, e.g. 1,2,...,12 on 6 cores
#              with SMT: physical cores alone, then hyperthreads)
#   REPS       repetitions for medians (default 3)
#   MAX_READS  records of SRR5833294 inserted (default 2000000, as in the paper)
#   CXX        C++ compiler (default g++, needs C++20 and <atomic> wait/notify)
#   NO_FETCH=1 fail instead of downloading missing data
#   READS_SOURCE  sra (default, NCBI SRA via the SRA Toolkit installed into
#              external/sratoolkit) or ena; see tools/fetch_data.sh
#
# For stable numbers: plug in, close other programs, and set the
# performance governor (sudo cpupower frequency-set -g performance). The
# script only reports the governor and load; it never changes system settings.
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT=$(pwd)

ALL_STEPS=tests,kache,scale,eager,ablation,defrag,churn
STEPS=${STEPS:-$ALL_STEPS}
REPS=${REPS:-3}
MAX_READS=${MAX_READS:-2000000}
export REPS MAX_READS
# Every thread count from 1 to the hardware threads, unless given.
THREADS=${THREADS:-$(seq -s, 1 "$(getconf _NPROCESSORS_ONLN)")}
export THREADS
[[ -n "${CXX:-}" ]] && export CXX

host=$(uname -n | cut -d. -f1)
OUT=${OUT:-bench/results/${host}-$(date +%F)}
mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)
LOG="$OUT/paper_runs.log"
export RESULTS_DIR="$OUT"

say() { printf '\n\033[1m== %s\033[0m\n' "$*" | tee -a "$LOG"; }
note() { echo "$*" | tee -a "$LOG"; }
die() { echo "error: $*" | tee -a "$LOG" >&2; exit 1; }
want() { [[ ",$STEPS," == *",$1,"* ]]; }

for s in ${STEPS//,/ }; do
    [[ ",$ALL_STEPS," == *",$s,"* ]] || die "unknown step '$s' (steps: $ALL_STEPS)"
done

note "paper_runs.sh started $(date -Is) on $host; results in $OUT; steps $STEPS"

# ---------------------------------------------------------------- 1. tools
say "checking build tools"
missing=()
for c in git curl python3 getconf; do
    command -v "$c" >/dev/null || missing+=("$c")
done
CXX_BIN=${CXX:-g++}
command -v "$CXX_BIN" >/dev/null || missing+=("$CXX_BIN")
if [[ ${#missing[@]} -gt 0 ]]; then
    die "missing: ${missing[*]} (Debian/Ubuntu: sudo apt install git curl python3 g++ zlib1g-dev; Fedora: sudo dnf install git curl python3 gcc-c++ zlib-devel)"
fi

# The C++ drivers need C++20 atomics and zlib.
probe=$(mktemp -d)
cat > "$probe/p.cpp" <<'CPP'
#include <atomic>
#include <zlib.h>
int main() { std::atomic_flag f; f.test(); return zlibVersion()[0] == 0; }
CPP
if ! "$CXX_BIN" -std=c++20 "$probe/p.cpp" -o "$probe/p" -lz 2>"$probe/err"; then
    cat "$probe/err" >&2
    rm -rf "$probe"
    die "$CXX_BIN cannot build a C++20 program against zlib (need GCC >= 11 and zlib headers: zlib1g-dev / zlib-devel)"
fi
rm -rf "$probe"

if ! command -v cargo >/dev/null; then
    export RUSTUP_HOME="$ROOT/external/rust/rustup" CARGO_HOME="$ROOT/external/rust/cargo"
    export PATH="$CARGO_HOME/bin:$PATH"
    if ! command -v cargo >/dev/null; then
        note "cargo not found: installing a private Rust toolchain into external/rust"
        mkdir -p external/rust
        curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
            | sh -s -- -y --no-modify-path --profile minimal --default-toolchain stable
    fi
fi
note "cargo: $(command -v cargo) ($(cargo -V))"

# ---------------------------------------------------------------- 2. data
say "checking data sets"
# Only the records of SRR5833294 that the runs use are downloaded: the first
# MAX_READS are inserted and the next MAX_READS queried as held-out reads.
if [[ "$MAX_READS" == "0" ]]; then READS_RECORDS=0; else READS_RECORDS=$((2 * MAX_READS)); fi
export READS_RECORDS
have_data() {
    [[ -s data/ecoli_NC_000913.3.fa && -s data/chr21.fa.gz ]] || return 1
    compgen -G "data/SRR5833294*.fastq.gz" >/dev/null || return 1
    compgen -G "data/SRR5901135*.fastq.gz" >/dev/null || return 1
    # A partial download made by fetch_data.sh must cover READS_RECORDS
    # (no .records file: a whole run fetched otherwise).
    local have
    have=$(cat data/SRR5833294.fastq.gz.records 2>/dev/null || echo whole)
    [[ "$have" == "whole" ]] || [[ "$READS_RECORDS" != "0" && "$have" -ge "$READS_RECORDS" ]]
}
if ! have_data; then
    [[ "${NO_FETCH:-0}" == "1" ]] && die "data sets missing in data/ and NO_FETCH=1"
    note "data sets missing: running tools/fetch_data.sh --reads (READS_RECORDS=$READS_RECORDS, READS_SOURCE=${READS_SOURCE:-sra})"
    df -h . | tee -a "$LOG"
    bash tools/fetch_data.sh --reads 2>&1 | tee -a "$LOG"
    have_data || die "tools/fetch_data.sh finished but data sets are still missing"
fi
ls -l data/ | tee -a "$LOG"

# ---------------------------------------------------------------- 3. external tools
say "building unihash and kache-hash (external/)"
cargo build --release --quiet
bash tools/kache_compare.sh --setup 2>&1 | tee -a "$LOG"

# ---------------------------------------------------------------- 4. machine record
say "recording the machine"
{
    echo "# Machine record, $(date -Is)"
    echo "host: $(uname -n)"
    echo "kernel: $(uname -srmo)"
    [[ -r /etc/os-release ]] && grep -E '^PRETTY_NAME=' /etc/os-release | sed 's/PRETTY_NAME=/os: /'
    echo "hardware threads: $(getconf _NPROCESSORS_ONLN)"
    echo
    echo "## lscpu"
    lscpu 2>/dev/null || cat /proc/cpuinfo
    echo
    echo "## memory"
    free -h 2>/dev/null || grep -E 'Mem(Total|Available)' /proc/meminfo
    echo
    echo "## frequency governor and power"
    for g in /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor; do
        [[ -r "$g" ]] && cat "$g"
    done | sort | uniq -c || true
    [[ -r /sys/devices/system/cpu/intel_pstate/no_turbo ]] && echo "intel no_turbo: $(cat /sys/devices/system/cpu/intel_pstate/no_turbo)"
    [[ -r /sys/devices/system/cpu/cpufreq/boost ]] && echo "cpufreq boost: $(cat /sys/devices/system/cpu/cpufreq/boost)"
    for p in /sys/class/power_supply/*/online; do
        [[ -r "$p" ]] && echo "$(basename "$(dirname "$p")") online: $(cat "$p")"
    done
    echo
    echo "## toolchain"
    rustc -vV
    cargo -V
    "$CXX_BIN" --version | head -1
    python3 --version
    echo "kache-hash: $(git -C external/kache-hash rev-parse HEAD)$( [[ "${KACHE_PATCH:-1}" == "1" ]] && echo " + tools/kache/kache-hash-lock.patch")"
    [[ -x external/sratoolkit/bin/fastq-dump ]] && echo "SRA Toolkit: $(external/sratoolkit/bin/fastq-dump --version | grep -m1 .)"
    echo "SRR5833294 records in data/: $(cat data/SRR5833294.fastq.gz.records 2>/dev/null || echo 'whole run')"
    echo
    echo "## unihash source"
    echo "commit: $(git rev-parse HEAD 2>/dev/null || echo 'not a git checkout')"
    echo "uncommitted changes:"
    git status --porcelain 2>/dev/null | grep -v -E "^\?\? ${OUT#"$ROOT"/}" || echo "  none"
    echo
    echo "## settings"
    echo "STEPS=$STEPS REPS=$REPS MAX_READS=$MAX_READS THREADS=$THREADS"
} > "$OUT/machine.txt" 2>&1
{ sed -n '2,5p' "$OUT/machine.txt"; grep -E '^(Model name|Mem:|commit: )' "$OUT/machine.txt"; } | tee -a "$LOG" || true
gov=$(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo unknown)
[[ "$gov" == "performance" ]] || note "note: CPU governor is '$gov' (performance is steadier: sudo cpupower frequency-set -g performance)"
if [[ -n "$(git status --porcelain -- src tools Cargo.toml Cargo.lock 2>/dev/null)" ]]; then
    note "note: src/, tools/ or Cargo files have uncommitted changes; commit first so the results map to a commit"
fi

# ---------------------------------------------------------------- 5. steps
B=./target/release/unihash-bench
SRR=$(ls data/SRR5833294*.fastq.gz | head -1)
ECOLI=data/ecoli_NC_000913.3.fa
CHR21=data/chr21.fa.gz

step() { # name: runs step_<name> unless done
    local name=$1
    want "$name" || return 0
    if [[ -e "$OUT/.done-$name" && "${FORCE:-0}" != "1" ]]; then
        note "skipping $name (done $(cat "$OUT/.done-$name"); FORCE=1 to re-run)"
        return 0
    fi
    say "step $name  (load average: $(cut -d' ' -f1-3 /proc/loadavg))"
    local t0=$SECONDS
    "step_$name" 2>&1 | tee -a "$LOG"
    local dt=$((SECONDS - t0))
    date -Is > "$OUT/.done-$name"
    printf '%-9s %6d s  finished %s\n' "$name" "$dt" "$(date -Is)" >> "$OUT/timings.txt"
    note "step $name finished in ${dt} s"
}

step_tests() {
    cargo test --release --quiet
}

step_kache() {
    # Tables: memory, single-thread throughput, per-query counters,
    # concurrent insertion, concurrent queries.
    SINGLE=1 SCALE=1 bash tools/kache_compare.sh
}

step_scale() {
    # Table: ConcurrentKmerTable scaling (insertion, deletion, speed-ups);
    # KmerTable deletion rate; concurrent table memory and into_table;
    # read locks vs read snapshot; in-place defrag vs round trip.
    CAPACITY=1 DELETE=0.25 EXTRA="--defrag-compare ${EXTRA:-}" bash tools/scale.sh
}

step_eager() {
    mkdir -p "$OUT/eager"
    bash tools/eager_compare.sh | tee "$OUT/eager/summary.txt"
}

summ() { # file: one summary line from a unihash-bench output
    # Rates in M k-mers/s on the inserted input (first query section) and
    # on random negatives; candidates = candidate k-mers per in-order query.
    local f=$1
    rate() { awk -v p="$1" '$1 == p { print $4; exit }' "$f"; }
    printf 'insert %s  %s  pieces %s  entries %s  hot %s  in-order %s  random %s  streaming %s  negatives %s  candidates %s\n' \
        "$(grep -m1 -E '^table ' "$f" | awk '{print $2" s", $4" M/s"}')" \
        "$(grep -m1 'table total, accounted' "$f" | awk '{print $6" B/k-mer"}')" \
        "$(grep -m1 -E '^  pieces ' "$f" | awk '{print $2}')" \
        "$(grep -m1 -E '^  index entries ' "$f" | awk '{print $3}')" \
        "$(grep -m1 'hot minimizers' "$f" | awk '{print $6"/"$8}')" \
        "$(rate kmers-in-order)" "$(rate kmers-random-order)" "$(rate seq-streaming)" \
        "$(sed -n '/^-- negative queries/,/^--/p' "$f" | awk '$1 == "table" { print $4; exit }')" \
        "$(sed -n '/^k-mers in order/,/^[a-z]/p' "$f" | awk '/candidate k-mers/ { print $5; exit }')"
}

step_ablation() {
    local d="$OUT/ablation"
    mkdir -p "$d"
    : > "$d/summary.txt"
    # Hot-minimizer routing: T = 0 (off), 16, 64 (default), 256.
    for t in 0 16 64 256; do
        echo "-- chr21 --hot-threshold $t"
        "$B" "$CHR21" --no-oracle --hot-threshold "$t" > "$d/chr21_hot$t.txt"
        echo "chr21 T=$t  $(summ "$d/chr21_hot$t.txt")" | tee -a "$d/summary.txt"
    done
    for t in 0 64; do
        echo "-- E. coli --hot-threshold $t"
        "$B" "$ECOLI" --no-oracle --hot-threshold "$t" > "$d/ecoli_hot$t.txt"
        echo "ecoli T=$t  $(summ "$d/ecoli_hot$t.txt")" | tee -a "$d/summary.txt"
    done
    # Candidate scan: tie offsets (default) vs the full range, synthetic reads,
    # checked against the oracle.
    for s in ties range; do
        echo "-- synthetic reads --scan $s"
        "$B" --synthetic 4000000 --reads 400000 --scan "$s" > "$d/reads_scan_$s.txt"
        echo "reads scan=$s  $(summ "$d/reads_scan_$s.txt")" | tee -a "$d/summary.txt"
    done
    # Memory cost of bounding pieces at 4096 bases (default) vs 1 Mb.
    for m in 4096 1048576; do
        echo "-- chr21 --max-piece-len $m"
        "$B" "$CHR21" --no-oracle --max-piece-len "$m" > "$d/chr21_maxpiece$m.txt"
        echo "chr21 max_piece_len=$m  $(summ "$d/chr21_maxpiece$m.txt")" | tee -a "$d/summary.txt"
    done
    grep -h -E "MISMATCHES|FALSE POSITIVES" "$d"/reads_scan_*.txt || true
}

step_defrag() {
    # defrag() after bulk loading (paper: deletion and defragmentation).
    local d="$OUT/defrag"
    mkdir -p "$d"
    q=()
    for f in data/SRR5901135*.fastq.gz; do q+=(--query "$f"); done
    echo "-- SRR5833294 --defrag"
    "$B" "$SRR" "${q[@]}" --max-records "$MAX_READS" --holdout --defrag > "$d/SRR5833294_defrag.txt"
    echo "-- chr21 --defrag"
    "$B" "$CHR21" --defrag > "$d/chr21_defrag.txt"
    grep -h -E "^defrag|MISMATCHES|FALSE POSITIVES" "$d"/*.txt
}

step_churn() {
    CHECK=1 bash tools/churn_real.sh
}

for s in tests kache scale eager ablation defrag churn; do
    step "$s"
done

# ---------------------------------------------------------------- summary
say "summary"
{
    echo "# Paper measurements on $host ($(date -Is))"
    echo
    echo "Machine: see machine.txt. Commit: $(git rev-parse --short HEAD 2>/dev/null || echo '?')."
    echo
    echo "| paper table / text | files |"
    echo "|---|---|"
    echo "| memory, single-thread throughput, per-query counters | kache/<data>_{summary.md,unihash.txt,kache.txt} |"
    echo "| concurrent insertion, concurrent queries | kache/<data>_scale_{summary.md,unihash.txt,kache.txt} |"
    echo "| ConcurrentKmerTable scaling and deletion, KmerTable deletion | scale/<data>_capacity.txt |"
    echo "| cost of the read locks (locked queries vs read snapshot, same table) | scale/<data>_capacity.txt, section \"read locks\" |"
    echo "| concurrent defrag in place vs round trip through KmerTable | scale/<data>_capacity.txt, section \"defragmentation after deletion\" |"
    echo "| eager vs deferred run indexing | eager/summary.txt |"
    echo "| hot threshold, candidate scan, max_piece_len | ablation/summary.txt |"
    echo "| defrag after bulk loading | defrag/*.txt |"
    echo "| churn | churn/<data>.txt |"
    echo
    echo "Not re-run by this script: the lexicographic minimizer order (no longer an"
    echo "option of the code), the 64-byte block size sweep and the A/B of the"
    echo "insertion shortcuts (comparisons against earlier commits), and the"
    echo "lock-free read path that was tried and dropped (not in the code; the"
    echo "read-lock section gives the most any lock-free read path could gain)."
    echo
    echo "## Step timings"
    echo
    echo '```'
    cat "$OUT/timings.txt" 2>/dev/null || true
    echo '```'
    for f in "$OUT"/kache/*summary.md; do
        [[ -e "$f" ]] || continue
        echo
        echo "## $(basename "$f")"
        echo
        cat "$f"
    done
    for f in "$OUT"/eager/summary.txt "$OUT"/ablation/summary.txt; do
        [[ -e "$f" ]] || continue
        echo
        echo "## ${f#"$OUT"/}"
        echo
        echo '```'
        cat "$f"
        echo '```'
    done
} > "$OUT/SUMMARY.md"

# unihash-bench reports wrong answers without failing; the other tools exit
# with an error on any failed check, which stops this script.
bad='MISMATCHES vs .*[[:space:]][1-9][0-9]* \(|MISMATCHES vs unihash oracle [1-9]|FALSE POSITIVES +[1-9]|MISMATCH:'
if grep -r -l -E "$bad" "$OUT" --include='*.txt' 2>/dev/null; then
    note "WARNING: the files above report mismatches or failed checks"
fi
note "done $(date -Is): $OUT/SUMMARY.md"
