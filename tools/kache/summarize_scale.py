#!/usr/bin/env python3
"""Thread-scaling summary of unihash, kache-hash and std HashSet<u64>.

usage: summarize_scale.py UNIHASH_SCALE_OUTPUT KACHE_SCALE_OUTPUT

UNIHASH_SCALE_OUTPUT: `unihash-scale --hashset --dump DIR ...`
KACHE_SCALE_OUTPUT:   `kache_scale DIR` on the same dump
Prints Markdown tables: insertion and memory, then one table per query
pattern (positive, random, streaming) by thread count.
"""
import re
import sys

ROW = re.compile(r"^\s*(\d+) \|(.*)$")


CPU_NOTE = re.compile(r"CPU check ([0-9.]+)×")
# Rows measured while the one-thread CPU speed was this far below the run's
# median are marked: the machine was throttled or busy.
SLOW = 0.85


def rows_after(lines, header_start, cpu):
    """Rows `threads | ...` of the table whose header line starts with
    `header_start` (first such table). Each row's CPU check (last column when
    the header has one, else the diagnostics line under the row) goes to
    `cpu[(header_start, threads)]`."""
    out = {}
    for i, line in enumerate(lines):
        if line.startswith(header_start):
            last_col = line.rstrip().endswith("CPU check")
            t = None
            for row in lines[i + 1:]:
                m = ROW.match(row)
                if not m:
                    if row.startswith("  "):  # diagnostics under a row
                        c = CPU_NOTE.search(row)
                        if c and t is not None:
                            cpu[(header_start, t)] = float(c.group(1))
                        continue
                    break
                cells = [c.strip() for c in m.group(2).split("|")]
                if any("skipped" in c for c in cells):
                    t = None
                    continue
                t = int(m.group(1))
                out[t] = [float(c.rstrip(" %")) for c in cells]
                if last_col:
                    cpu[(header_start, t)] = out[t][-1]
            return out
    return out


def single(text, name):
    m = re.search(
        name + r": (\d+) k-mers in ([0-9.]+) s, ([0-9.]+) M k-mers/s, (\d+) distinct k-mers, ([0-9.]+) B/k-mer",
        text)
    if not m:
        sys.exit(f"summarize_scale: no '{name}' line")
    return {"time": float(m.group(2)), "rate": float(m.group(3)), "distinct": int(m.group(4)),
            "bytes": float(m.group(5))}


def fmt(x, d=2):
    return "—" if x is None else f"{x:.{d}f}"


def main():
    u_text = open(sys.argv[1]).read()
    k_text = open(sys.argv[2]).read()
    u, k = u_text.splitlines(), k_text.splitlines()
    if "check: all kache-hash sizes and answers equal unihash's" not in k_text:
        sys.exit("summarize_scale: kache-hash answers were not all checked equal")

    kt = single(u_text, "single-threaded KmerTable")
    hs = single(u_text, "single-threaded HashSet<u64>")
    ks = single(k_text, "single-threaded kache-hash")
    if not kt["distinct"] == hs["distinct"] == ks["distinct"]:
        sys.exit("summarize_scale: the three structures hold different numbers of k-mers")

    # threads -> columns
    cpu = {}
    H_CONC = "threads | insert s | M k-mers/s | speed-up vs 1 thread | vs KmerTable"
    H_SHARD = "threads | insert s | M k-mers/s | speed-up vs 1 thread | B/k-mer (heap) | peak"
    H_KCONC = "threads | insert s | M k-mers/s | speed-up vs 1 thread | B/k-mer (heap)"
    H_QU, H_QH, H_QK = "threads | KmerTable positive", "threads | HashSet positive", "threads | kache positive"
    conc = rows_after(u, H_CONC, cpu)
    shard = rows_after(u, H_SHARD, cpu)
    kconc = rows_after(k, H_KCONC, cpu)
    q_u = rows_after(u, H_QU, cpu)
    q_h = rows_after(u, H_QH, cpu)
    q_k = rows_after(k, H_QK, cpu)
    threads = sorted(set(conc) | set(kconc) | set(shard))
    checks = sorted(cpu.values())
    ref = checks[len(checks) // 2] if checks else None
    marked = []

    def mark(header, t):
        c = cpu.get((header, t))
        if ref and c is not None and c < SLOW * ref:
            marked.append(c / ref)
            return " †"
        return ""

    def ins(row, header, t):
        return None if row is None else f"{row[1]:.2f} ({row[0]:.2f} s){mark(header, t)}"

    print("Insertion, M k-mers/s (seconds); every structure holds the same "
          f"{kt['distinct']} distinct k-mers.\n")
    print("| threads | unihash | kache-hash | HashSet<u64> |")
    print("|---|---|---|---|")
    print(f"| 1, single-threaded structure | {kt['rate']:.2f} ({kt['time']:.2f} s) | "
          f"{ks['rate']:.2f} ({ks['time']:.2f} s) | {hs['rate']:.2f} ({hs['time']:.2f} s) |")
    for t in threads:
        print(f"| {t}, concurrent structure | {ins(conc.get(t), H_CONC, t) or '—'} | "
              f"{ins(kconc.get(t), H_KCONC, t) or '—'} | {ins(shard.get(t), H_SHARD, t) or '—'} |")
    print()
    print("Memory, B per distinct k-mer (heap after the build):\n")
    print("| | unihash | kache-hash | HashSet<u64> |")
    print("|---|---|---|---|")
    print(f"| single-threaded structure | {kt['bytes']:.3f} | {ks['bytes']:.3f} | {hs['bytes']:.3f} |")
    t_max = max(threads)
    c, kc, sc = conc.get(t_max), kconc.get(t_max), shard.get(t_max)
    print(f"| concurrent structure ({t_max} threads) | {fmt(c and c[5], 3)} | {fmt(kc and kc[3], 3)} | "
          f"{fmt(sc and sc[3], 3)} |")
    if c:
        print(f"| unihash concurrent table after `into_table` | {c[8]:.3f} | | |")
    print()

    for i, name in enumerate(["positive k-mers (shuffled)", "uniform random k-mers", "streaming (sequences)"]):
        print(f"Queries, {name}, M queries/s:\n")
        print("| threads | KmerTable | ConcurrentKmerTable | kache-hash | HashSet<u64> |")
        print("|---|---|---|---|---|")
        for t in sorted(set(q_u) | set(q_k) | set(q_h)):
            a = q_u.get(t)
            b = q_k.get(t)
            h = q_h.get(t)
            mu, mk, mh = mark(H_QU, t), mark(H_QK, t), mark(H_QH, t)
            print(f"| {t} | {fmt(a and a[i])}{mu} | {fmt(a and a[3 + i])}{mu} | {fmt(b and b[i])}{mk} | "
                  f"{fmt(h and h[i])}{mh} |")
        print()
    if marked:
        print(f"† measured while the one-thread CPU speed was {100 * min(marked):.0f}–{100 * max(marked):.0f} % "
              "of the run's median (CPU check): the machine was throttled or busy; repeat before "
              "drawing conclusions from these cells.\n")
    presized = "with_capacity(" in u_text
    print("Single-threaded structures grow from empty (kache-hash from its initial capacity); "
          + ("the concurrent ones are created for the distinct k-mer count (unihash "
             "`with_capacity`, kache-hash's `max_sz`, `HashSet::with_capacity` per shard). "
             if presized else "the concurrent ones grow too. ")
          + "KmerTable, kache-hash and HashSet are queried after the single-threaded build, shared "
          "read-only (no locks); ConcurrentKmerTable takes its shard read locks. The concurrent "
          "HashSet is 1024 shards of Mutex<HashSet<u64>> filled through per-thread buffers. "
          "HashSet streaming: independent lookups of the k-mers parsed from the same sequences.")


if __name__ == "__main__":
    main()
