#!/usr/bin/env python3
"""Side-by-side summary of an unihash-bench run and a kache_bench run.

usage: summarize.py UNIHASH_OUTPUT KACHE_OUTPUT
Prints a Markdown table (throughput in M k-mers/s, memory in B/distinct k-mer).
"""
import re
import sys


def rate(tok_value, tok_unit):
    v = float(tok_value)
    return v if tok_unit.startswith("M") else v / 1000.0


RATE = r"([0-9.]+) ([Mk])/s"


def parse_unihash(text):
    r = {"queries": {}}
    section = None
    for line in text.splitlines():
        if line.startswith("-- "):
            section = line.strip("- ").strip()
            continue
        m = re.match(r"^(table|oracle)\s+[0-9.]+ s\s+" + RATE, line)
        if m and section == "insertion":
            r["insert_" + m.group(1)] = rate(m.group(2), m.group(3))
        if m and section and section.startswith("negative"):
            r["neg_" + m.group(1)] = rate(m.group(2), m.group(3))
        m = re.match(r"^table heap, allocator\s+\d+ B\s+([0-9.]+) B/k-mer", line)
        if m:
            r["mem_table"] = float(m.group(1))
        m = re.match(r"^oracle HashSet<u64>, allocator\s+\d+ B\s+([0-9.]+) B/k-mer", line)
        if m:
            r["mem_oracle"] = float(m.group(1))
        m = re.match(r"^distinct k-mers\s+(\d+)", line)
        if m:
            r["distinct"] = int(m.group(1))
        if section and section.startswith("queries: "):
            label = section[len("queries: "):]
            q = r["queries"].setdefault(label, {})
            m = re.match(r"^(kmers-in-order|kmers-random-order|seq-independent|seq-streaming)\s+[0-9.]+ s\s+" + RATE + r"(?:\s+[0-9.]+ s\s+" + RATE + ")?", line)
            if m:
                q[m.group(1)] = rate(m.group(2), m.group(3))
                if m.group(4):
                    q[m.group(1) + "/oracle"] = rate(m.group(4), m.group(5))
    return r


def parse_kache(text):
    r = {"queries": {}}
    section = None
    for line in text.splitlines():
        if line.startswith("-- "):
            section = line.strip("- ").strip()
            continue
        m = re.match(r"^kache\s+[0-9.]+ s\s+" + RATE, line)
        if m and section == "insertion":
            r["insert"] = rate(m.group(1), m.group(2))
        if m and section and section.startswith("negative"):
            r["neg"] = rate(m.group(1), m.group(2))
        m = re.match(r"^kache heap \(glibc\)\s+(\d+) B", line)
        if m:
            r["heap"] = int(m.group(1))
        m = re.match(r"^records \d+, k-mers \d+, distinct (\d+)", line)
        if m:
            r["distinct"] = int(m.group(1))
        if section and section.startswith("queries: "):
            label = section[len("queries: "):]
            q = r["queries"].setdefault(label, {})
            m = re.match(r"^(kmers-in-order|kmers-random-order|seq-streaming)\s+[0-9.]+ s\s+" + RATE, line)
            if m:
                q[m.group(1)] = rate(m.group(2), m.group(3))
            if "MISMATCHES" in line:
                q["mismatches"] = int(line.split()[-1])
    return r


def fmt(x, digits=2):
    return "—" if x is None else f"{x:.{digits}f}"


def main():
    u = parse_unihash(open(sys.argv[1]).read())
    k = parse_kache(open(sys.argv[2]).read())
    kmem = k["heap"] / k["distinct"] if k.get("heap") and k.get("distinct") else None
    rows = [
        ("Memory, B per distinct k-mer", u.get("mem_table"), kmem, u.get("mem_oracle"), 2),
        ("Insert, M k-mers/s", u.get("insert_table"), k.get("insert"), u.get("insert_oracle"), 2),
    ]
    for label, q in u["queries"].items():
        kq = k["queries"].get(label, {})
        rows += [
            (f"{label}: k-mers in order", q.get("kmers-in-order"), kq.get("kmers-in-order"), q.get("kmers-in-order/oracle"), 2),
            (f"{label}: k-mers random order", q.get("kmers-random-order"), kq.get("kmers-random-order"), q.get("kmers-random-order/oracle"), 2),
            (f"{label}: sequence, streaming", q.get("seq-streaming"), kq.get("seq-streaming"), q.get("seq-independent/oracle"), 2),
        ]
        if kq.get("mismatches"):
            print(f"WARNING: kache-hash answers differ from the oracle on '{label}'", file=sys.stderr)
    rows.append(("Random negatives, M k-mers/s", u.get("neg_table"), k.get("neg"), u.get("neg_oracle"), 2))
    print("| | unihash | kache-hash | HashSet<u64> |")
    print("|---|---|---|---|")
    for name, a, b, c, d in rows:
        print(f"| {name} | {fmt(a, d)} | {fmt(b, d)} | {fmt(c, d)} |")
    print()
    print("Streaming row, HashSet column: independent lookups parsed from the same sequence.")


if __name__ == "__main__":
    main()
