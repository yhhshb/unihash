#!/usr/bin/env bash
# Local CI for unihash.
#
# Memory and undefined-behaviour safety is checked as follows (in place of
# ASan/UBSan):
#   * the crate contains no `unsafe` code (enforced by #![forbid(unsafe_code)]),
#     so out-of-bounds accesses are caught by bounds checks (panics);
#   * tests run with overflow-checks and debug-assertions ON, so any shift by
#     ≥ 64 or arithmetic overflow panics instead of being undefined;
#   * optionally, with a nightly toolchain, the suite also runs under Miri
#     (UB interpreter) and AddressSanitizer — set MIRI=1 / ASAN=1.
#
# Usage: bash tools/ci.sh            # fmt, clippy, tests (checked + release), smoke runs
#        MIRI=1 ASAN=1 bash tools/ci.sh
set -euo pipefail
cd "$(dirname "$0")/.."

step() { printf '\n\033[1m== %s\033[0m\n' "$*"; }

step "rustfmt"
cargo fmt --all -- --check

step "clippy (deny warnings)"
cargo clippy --all-targets -- -D warnings

step "tests: dev profile (overflow checks + debug assertions, unoptimised)"
cargo test --quiet --lib
cargo test --quiet --doc
cargo test --quiet --test m1_codec --test m2_minimizer --test m3_store_index --test m7_hardening

step "tests: test-opt profile (overflow checks + debug assertions, -O2) — full suite"
cargo test --quiet --profile test-opt

step "tests: release profile"
cargo test --quiet --release

step "benchmark smoke run (must report zero mismatches)"
cargo run --quiet --release --bin unihash-bench -- \
    --synthetic-repeats 300000 --reads 20000 --bucket-bits 14 --neg 100000 \
    | grep -E "MISMATCHES|FALSE POSITIVES"

step "churn smoke run (deletions, re-insertions, defrag; checked against the oracle)"
cargo run --quiet --release --bin unihash-churn -- \
    --synthetic-repeats 200000 --unit 300 --ops 3000 --lag 200 --check \
    --policy never --policy overhead:0.05 --policy inc:0.02+full:0.1 --policy every:500 \
    | grep -E "check vs oracle|estimate"

if [[ "${MIRI:-0}" == "1" ]]; then
    step "Miri (nightly) — unit + small tests"
    cargo +nightly miri test --lib --test m1_codec --test m3_store_index --test m7_hardening
fi

if [[ "${ASAN:-0}" == "1" ]]; then
    step "AddressSanitizer (nightly)"
    host=$(rustc -vV | sed -n 's/host: //p')
    RUSTFLAGS="-Zsanitizer=address" RUSTDOCFLAGS="-Zsanitizer=address" \
        cargo +nightly test --quiet --profile test-opt -Zbuild-std --target "$host" --tests
fi

step "all green"
