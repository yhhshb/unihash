# unihash

An exact, single-threaded set of canonical k-mers. The k-mers are stored as sequence
**pieces**, not as separate integers. Membership queries go through a **minimizer
index**. It supports streaming insertion, exact deletion and user-controlled
defragmentation.

```rust
use unihash::{Config, KmerTable};

let mut t = KmerTable::new(Config { k: 31, l: 17, ..Config::default() })?;
t.insert_read(b"ACGT...");            // FASTA/FASTQ sequence bytes; N etc. split the read

// Long records (chromosomes) can be streamed in chunks cut anywhere; the
// result is identical to insert_read on the concatenation.
t.insert_chunk(b"ACGTAC");
t.insert_chunk(b"GTTNNACG...");
t.end_read();
let hit: bool = t.contains_seq(b"ACGTACGTACGTACGTACGTACGTACGTACG")?; // either strand
let hit: bool = t.contains(packed_kmer);                              // 2-bit packed, either strand

// In-order queries of the successive k-mers of a sequence (chunks may be cut
// anywhere); much faster than independent lookups, same answers.
let mut q = t.query_stream();
q.push(b"ACGTTGCA...", |kmer, present| { /* canonical k-mer, answer */ });
q.end();                                                              // between sequences

// Deletion: every k-mer of a sequence (absent ones are ignored), or one k-mer.
let removed: u64 = t.delete_read(b"ACGTTGCA...");
let was_present: bool = t.delete(packed_kmer);

// Defragmentation is never automatic. Cheap and frequent: slide the base
// array over the bases freed by deletions (no index work). Rare: merge pieces
// that overlap by k-1 and rebuild the index.
let f = t.fragmentation();
if f.piece_overhead(t.k()) > 0.10 {
    println!("{}", t.defrag());
} else if f.incremental_overhead() > 0.02 {
    println!("{}", t.defrag_incremental());
}
println!("{}", t.stats());
```

## Build and test

```bash
cargo build --release
cargo test --profile test-opt     # full suite, overflow checks on, -O2 (about 5 s)
bash tools/ci.sh                  # fmt + clippy -D warnings + tests (dev, test-opt, release) + smoke runs
MIRI=1 ASAN=1 bash tools/ci.sh    # also Miri and AddressSanitizer (needs a nightly toolchain)
```

Replay a randomized test with one seed: `UNIHASH_SEED=1234 cargo test --profile test-opt <name>`.
Every assertion message includes its seed.

**Sanitizer equivalent.** The crate has `#![forbid(unsafe_code)]`, so out-of-bounds
accesses panic on bounds checks. All test profiles keep `overflow-checks` and
`debug-assertions` on, so a shift by 64 or an integer overflow panics instead of
being UB. `k = 32` and `ℓ = 32` are tested explicitly.

## Benchmarks

```bash
cargo run --release --bin unihash-bench -- --help
cargo run --release --bin unihash-bench -- --synthetic-repeats 4600000
bash tools/fetch_data.sh [--reads]  # E. coli, chr21 [+ SRR5833294, SRR5901135] into data/
bash tools/bench_real.sh            # benchmark matrix on real data → bench/results/
bash tools/kache_compare.sh         # same workloads on kache-hash → bench/results/kache/
bash tools/churn_real.sh            # deletion/re-insertion churn × defrag policies → bench/results/churn/
cargo run --release --bin unihash-churn -- --synthetic-repeats 1000000 --check
```

**Churn benchmark.** `unihash-churn` builds a table, then repeats "delete one
unit (a `--unit`-base window, or a whole read), re-insert the unit deleted
`--lag` operations earlier" under several defragmentation policies (`never`,
`every:N` operations, `overhead:X` = `defrag` when `fragmentation().overhead(k) > X`,
`inc:X+full:Y` = `defrag_incremental` when `incremental_overhead() > X` and `defrag`
when `piece_overhead(k) > Y`). For
each policy it reports edit time, defragmentation time, mean / max / final table
size and query speed in the final state; `--check` compares every policy with
the oracle — exact set equality in both directions (every stored k-mer is in the
HashSet, every HashSet k-mer is found) — and with `verify_index`. `unihash-bench --defrag` defragments right after
insertion.

**kache-hash comparison.** `tools/kache_compare.sh` clones
[kache-hash](https://github.com/jamshed/kache-hash) (pinned commit) into
`third_party/`, builds `tools/kache/kache_bench.cpp` against it (g++ ≥ 11, zlib,
AVX2), and for each data set runs `unihash-bench --dump DIR`, then `kache_bench DIR`
on the same inputs and the same query k-mers, then prints a side-by-side table
(`tools/kache/summarize.py`). kache-hash is built single-threaded with the same
k = 31 and ℓ = 17, and its answers are checked against the hit counts from the
unihash oracle. kache-hash must be compiled with `-fno-strict-aliasing`: at `-O3`
without that flag, GCC miscompiles it and it finds about 1 % of the present k-mers.
`bash tools/kache_compare.sh NAME ARGS…` runs a single comparison.

The bench tool reads FASTA/FASTQ, optionally gzipped, and reports:

* **Insertion** throughput, and the **table's memory on its own**, measured before
  the `HashSet<u64>` oracle is built in a second pass. Heap bytes are accounted
  from the structures' capacities and cross-checked against a counting global
  allocator (they agree to within a few hundred bytes). The output also gives the
  peak heap during insertion and RSS before and after.
* **Queries on uniform workloads.** Each query source (the inserted input, the
  held-out records with `--holdout`, and every `--query` file) is sampled
  uniformly: windows of `--window` consecutive bases spread at a fixed stride over
  the whole input, or every s-th record for reads, up to `--query-kmers` k-mers.
  All patterns use the same k-mers and are timed the same way for the table and the
  oracle:
  * pre-parsed k-mers in sequence order, and shuffled;
  * the sequence itself, as parsing only, independent lookups, and the streaming
    query.
* **Random negatives**, per-query counters (entries scanned, tag matches, candidate
  k-mers), a histogram of entries per minimizer, the hottest minimizers, and bucket
  occupancy.

It exits with status 1 if any answer differs from the oracle. FASTA records are
streamed in chunks, so a chromosome is never held in memory whole. Results are in
[`bench/RESULTS.md`](bench/RESULTS.md).

## Parameters

| Field | Default | Range | Meaning |
|---|---|---|---|
| `k` | 31 | 2 ≤ k ≤ 32 | k-mer length; a canonical k-mer fits in a `u64` |
| `l` | 17 | 1 ≤ ℓ ≤ k | minimizer length |
| `bucket_bits` | 10 | 0 ≤ b ≤ 30 | initial bucket count 2^b. The index grows by 1.5× whenever its bucket heads would be more than 85 % full, so it needs no sizing |
| `max_piece_len` | 4096 | k ≤ L ≤ 2^24−1 | longer runs are split into pieces that overlap by k−1 bases, and `defrag` merges up to this length. Bounds the cost of a deletion (the right part of a cut piece is re-filed) and the insertion buffers (about 14 B × L). A split costs ≈ 15 B, < 0.1 % of the table at 4096 |
| `scan` | `TieOffsets` | `TieOffsets` / `Range` | candidate positions checked per index entry |
| `hot_threshold` | 64 | 0 = off | entries per ℓ-minimizer (and per ℓ2-minimizer in layer 2) above which the minimizer is routed to the next layer |
| `l2` | 0 = min(k, ℓ+8) | ℓ < ℓ2 ≤ k | minimizer length of the layer-2 index |
| `spill_dir` | `None` | existing directory | where index growth writes the old entries (8 B each) while the new block array is built; `None` keeps them in a compact in-memory buffer |

Invalid parameters are rejected by `KmerTable::new` with a descriptive `ConfigError`.

## Correctness contract

After any sequence of `insert_read` calls (or `insert_chunk` … `end_read` sequences),
`delete_read` / `delete` calls and `defrag` calls, `contains(x)` is `true` exactly when the
canonical form of `x` occurs as a k-mer of some inserted read and has not been
deleted since its last insertion, where the k-mer
consists only of `ACGTacgt`. There are no false positives and no false negatives.
Answers do not depend on the order of insertion.

* Reads are split at every non-`ACGTacgt` character. Lowercase (soft-masked) bases
  are accepted and stored as uppercase. Segments shorter than `k` add nothing.
* While a read is open (`insert_chunk` called, `end_read` not yet), its last
  pending run of absent k-mers is not visible to queries.
* `contains` accepts either orientation. A packed value with bits set above `2k`
  returns `false`. `contains_seq` returns an error for a wrong length or a
  non-ACGT character.
* `delete_read`, `delete` and `defrag` panic while a read is open (call
  `end_read` first).

### How it works

* **Insertion**. Insertion is streaming: the rolling k-mer, the
  rolling ℓ-mer and a monotone deque of at most `w = k − ℓ + 1` ℓ-mer values carry
  over between chunks, and the only buffer that grows is the pending run (at most
  `max_piece_len` bases). The ℓ-mer values of a piece are computed when it is
  emitted. So memory does not depend on the length of a read or chromosome.
  Each segment is scanned with rolling k-mers. Present
  k-mers are skipped. Every maximal run of absent k-mers becomes a piece and is
  indexed before the scan continues, so later k-mers in the same read can see it.
  As a result each distinct k-mer lives in exactly one piece. The exception is a
  k-mer repeated *inside one pending run*, which is stored twice in the same piece.
  This redundancy is accepted and reported as
  `stored_kmers − distinct`. When a run reaches `max_piece_len` it is closed
  *before* the next lookup. This keeps a repeat from being duplicated across the
  split.
* **Minimizers**. `m(x)` is the minimum order value over the canonical
  ℓ-mers of `x`. The order value is a fixed bijective 64-bit mix of the canonical
  ℓ-mer (SplitMix64 finalizer), so ties are exactly those of the plain value and
  the hot minimizers can still be decoded back to their ℓ-mer. **Every** tied position is indexed. The order value depends only
  on the canonical ℓ-mer, so `m(x) = m(revcomp x)` and the tie sets mirror each
  other.
* **Pieces** are packed at 2 bits per base, back to back in one `u64` array. A
  piece is a view `(start, len)` into that array, packed into one `u64` (40-bit
  start, 24-bit length), so the cost is still 8 B per piece. Any k-mer (k ≤ 32)
  is read with two word loads and one shift, wherever it falls relative to word
  boundaries.
* **Index**. The index is one flat array of fixed-capacity 80 B blocks. Each block
  holds 8 entries as separate arrays of 8 tags, 8 piece ids and 8 offsets (9 B per
  entry), plus a fill count and an overflow link. There are no per-bucket
  allocations or headers.
  * **Placement.** Each minimizer has two candidate buckets, and a new entry goes
    to the less-full head ("two choices"). When both heads are full, the entry
    goes to an overflow chain.
  * **Lookup.** A lookup scans both chains, comparing 8-bit tags. A tag
    collision costs one exact comparison and never changes an answer.
  * **Growth.** Entries do not store `m`. When the heads pass 85 % load, the
    index grows by 1.5×, re-placing each entry with `m` recomputed from its piece.
    The old and new block arrays never coexist: the old entries are written as
    8 B `(pid, q)` records to a spill (a file in `spill_dir`, or a compact
    in-memory buffer), the old blocks are freed, the new ones allocated, and the
    records streamed back in. The spill file is deleted after each rebuild.
  * **Compaction.** After growth and in `shrink_to_fit`, overflow chains are laid
    out contiguously. Heads stay in place; only the overflow area is copied.
* **Hot minimizers** (three layers). Repeats (Alu, homopolymers, satellites) give
  some minimizers thousands of entries, and every query on such a k-mer would scan
  all of them.
  * **Layer 1** is the ℓ index above.
  * When an insertion overflows both heads of `m`, the entries of `m` are counted.
    At `hot_threshold` they are drained and `m` is marked hot (a flag bit on its
    first head, plus a small set). The k-mers those entries covered are re-indexed
    in **layer 2**: a second index of the same kind with longer minimizers
    (ℓ2 = ℓ+8 by default) over the **same pieces**. No sequence is copied, and
    contiguous runs stay pieces.
  * Hot ℓ2-minimizers in layer 2 go the same way into **layer 3**, a
    `HashSet<u64>` of full canonical k-mers.
  * A query computes `m`, then checks one bit. Cold minimizers take the usual
    path. Hot ones look up the k-mer's ℓ2-minimizer in layer 2, and, if that one
    is hot too, layer 3. The streaming query does the same; it caches candidates
    only for cold minimizers.
* **Query**. For each tag-matching entry `(pid, q)`:
  * `Range` checks every start in `[q−(k−ℓ), q]` that fits in the piece.
  * `TieOffsets` checks only `s = q − p` for `p ∈ Q(x) ∪ {(k−ℓ) − p'}`, where
    `p' ∈ Q(x)`. If `x` (or its reverse complement) occurs at `s`, then `q − s`
    is one of those offsets, so this mode is exact too. It usually needs 1–2
    comparisons instead of `k−ℓ+1`.
* **Streaming (in-order) queries** (`KmerTable::query_stream`). Successive k-mers
  share `k − 1` bases, so three things are reused from one k-mer to the next:
  1. The minimizer and its tie mask come from a rolling ℓ-mer and a monotone
     deque (O(1) amortised per base), not from the k ℓ-mers of each k-mer.
  2. If the previous k-mer was found at `(piece, s)`, the next one is at `s + 1`
     (same orientation) or `s − 1` (reverse complement) exactly when one stored
     base matches the incoming base or its complement. This single comparison is
     exact and needs no index access.
  3. While the minimizer stays the same, the tag-matching candidates of its
     bucket are cached, so the bucket is scanned once per minimizer run. Only the
     exact candidate comparisons remain.

  Answers are identical to `contains` for every k-mer (tested against the oracle,
  with random chunking, cap splits, collisions and small k).

* **Deletion** (`delete_read`, `delete`). The k-mers of the query are looked up
  in order; successive k-mers found at successive positions of one piece (either
  strand, checked with the streaming query's one-base extension) form a run, and
  each run is cut out of its piece in one operation:
  * the piece `P[0, L)` minus the k-mers at starts `lo..=hi` becomes a **left
    part** `P[0, lo + k − 1)`, which keeps `P`'s id (only its length shrinks, so
    its index entries stay valid), and a **right part** `P[hi + 1, L)`, a new
    piece *viewing the same bases*. Each part keeps the k − 1 bases that flank
    the run, so every surrounding k-mer survives. Nothing is copied; the bases
    of the run that neither part covers become dead space;
  * index maintenance is exact. The index items of a window (layer-1 entries,
    layer-2 entries, or its layer-3 location) are a function of its bases and
    of which minimizers are hot, so they are recomputed for the windows around
    the cut: items only deleted windows needed are removed; the right part's
    items are re-filed under its id in place (the bucket depends only on the
    minimizer), or copied when the left part needs them too. Layer 3 now maps
    each k-mer to its piece location, so a layer-3 k-mer can be cut as well;
  * a k-mer repeated inside one inserted run is stored twice, so `delete_read`
    re-checks its input with a streaming query and deletes leftover copies.

  The cost of a cut is proportional to the right part (≤ `max_piece_len`), plus
  the lookups of the query's k-mers.
* **Defragmentation** (`defrag`, never automatic). Pieces whose ends overlap by
  exactly k − 1 bases, in either orientation, are paired (canonical (k−1)-mer of
  each end, sorted; union-find keeps chains acyclic and within
  `max_piece_len`), and each chain is written once into a new base array
  without dead space. Merging two such pieces adds no k-mer: the junction
  windows are the pieces' own end windows. The indexes are freed first and
  rebuilt from the new pieces with the insertion code path; hot minimizers stay
  hot. `fragmentation()` reports dead and shared bases, extra pieces since the
  last defrag, and `overhead(k)`: the estimated fraction of the table `defrag`
  would reclaim, including index slack (block bytes beyond the bytes per entry
  of the last index rebuild). The indexes are freed before the piece ends are
  paired, so the defrag peak stays at about the table's size.
* **Incremental defragmentation** (`defrag_incremental`, never automatic).
  Most of what deletions free is dead bases, and those do not need the index:
  pieces keep their ids and lengths, so sliding the base array over the dead
  bases (32-base word moves, in place) only changes piece starts. It also drops
  unlinked overflow blocks and re-sizes an index layer whose slack reached 1/16
  of it (O(entries) rebuild, no windows recomputed). ≈ 12 ms on chr21 against
  2 s for `defrag`. Merging is left to `defrag`:
  `defrag_incremental_with(max_mean_len)` can also join the pieces changed
  since the last defragmentation with their k − 1-overlapping neighbours
  (found with 8 index lookups per piece), but a join costs ≈ 30 µs for ≈ 28 B,
  5× more per reclaimed byte than a full `defrag`.
* **When to call which** (`bench/RESULTS.md`). `fragmentation()` splits the
  estimate: `incremental_overhead()` (dead bases, re-sizable index slack) and
  `piece_overhead(k)` (extra pieces, shared bases, emptied records). Measured
  guidance: `defrag` once after bulk-loading reads (−18 % on SRR5833294); then
  `defrag_incremental` when `incremental_overhead() > 0.02` and `defrag` when
  `piece_overhead(k) > 0.10`. On chr21 that kept the same mean memory as
  calling `defrag` at 10 % overhead, for 0.7 % of the edit time instead of
  10.5 %.

## Layout

```
src/codec.rs          2-bit codec, rolling fwd/rc k-mers, reverse complement
src/piece_store.rs    packed pieces (2-bit bases + 8 B views), id reuse, base compaction
src/minimizer.rs      canonical ℓ-mer values, monotone-deque minima, tied positions
src/index.rs          flat fixed-capacity two-choice index with 8-bit tags, spilled
                      growth, hot-minimizer drain
src/table.rs          KmerTable: insert_read / insert_chunk / end_read / contains / stats,
                      layer-2 / layer-3 routing
src/edit.rs           delete_read / delete (piece cuts), defrag, defrag_incremental,
                      fragmentation, verify_index
src/stream_query.rs   StreamingQuery: in-order queries with hit extension and candidate cache
src/oracle.rs         HashSet<u64> oracle
src/fastx.rs          FASTA/FASTQ(.gz) reader;  src/rng.rs  SplitMix64 + synthetic genomes
src/bin/unihash_bench.rs  benchmark driver
src/bin/unihash_churn.rs  deletion / re-insertion churn × defragmentation policies
tests/m{1,...,7}_*.rs     codec, minimizers, store and index, table against the oracle
                          (invariance, small-k fuzz, chunked-insertion equivalence,
                          buffer bound), streaming queries, hot routing and disk
                          spill, hardening
tests/m8_delete.rs        deletion, alternating insert/delete, defrag, small-k exhaustive
tools/              ci.sh, fetch_data.sh, bench_real.sh, kache_compare.sh, churn_real.sh
tools/kache/        kache_bench.cpp (kache-hash driver), summarize.py
```

## Design notes

* **Testing.** The built-in test harness checks every operation against a
  `HashSet<u64>` oracle; `unihash-bench` measures end-to-end throughput and
  memory rather than micro-benchmarks. Overflow checks and the absence of
  `unsafe` stand in for ASan and UBSan; Miri and ASan can still be run on
  nightly.
* **Piece storage.** Pieces are views into one packed 2-bit array rather than a
  string each; packing was adopted once the random-order measurements justified
  it.
* **Index layout.** Flat fixed-capacity two-choice buckets with automatic
  growth, rather than a `Vec` per bucket and a fixed bucket count. Growth
  recomputes minimizers from pieces, so entries stay 9 B.
* **Minimizer order.** ℓ-mers are ranked by a bijective hash of their canonical
  value, not by the value itself. Lexicographic order was benchmarked and
  dropped: on chr21 it made insertion 4.8× and positive queries 2.5× slower (see
  `bench/RESULTS.md`).
* **Streaming insertion.** `insert_chunk` / `end_read` never allocate arrays the
  size of a segment. Peak RSS for chr21 insertion (without the oracle) fell from
  493 MiB to 162 MiB, for a 112 MiB table, when this was introduced.
* **Candidate scan.** `TieOffsets` is the default; the exhaustive `Range` scan
  is kept for comparison.
* **Deletion.** Pieces are views into the base array so that a cut copies
  nothing, and `max_piece_len` defaults to 4096 to bound the cost of a cut.
* **Piece minimizer positions.** Each piece's `(m, q)` set is computed with one
  deque pass and a mark array, not a sort followed by dedup.

## Known limits

* Single-threaded.
* Deletion leaves dead bases and extra pieces until `defrag` is called; `defrag`
  is O(table) and briefly holds the old and new base arrays (the indexes are
  freed before rebuilding).
* Index growth rebuilds the whole index. It no longer holds two block arrays,
  but it still costs 8 B per entry of spill (in memory, unless `spill_dir` is
  set) and a copy of the overflow area during compaction. Peak heap on chr21 is
  1.2× the final table with a disk spill and 1.4× with the in-memory spill.
* At most 2^32−1 pieces, each at most `max_piece_len` ≤ 2^24−1 bases, and at most
  2^40 bases in the base array (dead bases included until `defrag`). Exceeding
  any of these panics with a message.
* Hot routing is decided during insertion and never undone (deletion does not
  demote a minimizer). A minimizer that is hot early stays in layer 2 even if it
  would fit layer 1 after growth. Layer 3 stores whole k-mers with their location
  (about 25 B each); it is small on real data (23 k k-mers on chr21).
* Memory. Pieces cost 2 bits per base plus 8 B per piece, and index entries cost
  about 10–14 B each at 75–85 % head load (overflow blocks included). When pieces are short (reads), the
  index and the per-piece records dominate; `defrag` after loading reads merges
  overlapping pieces (−20 % on the synthetic read set).

## License

MIT — see [`LICENSE`](LICENSE).
