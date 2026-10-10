#!/usr/bin/env python3
"""Split the example reference (one record) into two records at a point no feature spans, to
test that gene/region placement is multi-record-aware: writes FASTA and GenBank files with the two parts in order (partA,
partB) and swapped (partB, partA), and split.json with the split point and the parts' lengths.

Usage: make_split_reference.py EXAMPLE_FASTA EXAMPLE_GBK OUT_DIR [TARGET_POSITION]"""
import json
import os
import re
import sys

NUM = re.compile(r"\d+")


def read_fasta(path):
    with open(path) as fh:
        lines = fh.read().split("\n")
    return "".join(l.strip() for l in lines[1:] if l and not l.startswith(">"))


def features(gb_lines):
    """(header lines before FEATURES, [(key, location, lines)])."""
    k = next(i for i, l in enumerate(gb_lines) if l.startswith("FEATURES"))
    end = next(i for i, l in enumerate(gb_lines) if l.startswith("ORIGIN") or l.startswith("CONTIG") or l == "//")
    feats, cur = [], None
    for line in gb_lines[k + 1:end]:
        if len(line) > 5 and line[5] != " ":
            cur = [line]
            feats.append(cur)
        elif cur is not None:
            cur.append(line)
    out = []
    for f in feats:
        loc = f[0][21:].strip()
        for l in f[1:]:
            t = l.strip()
            if t.startswith("/"):
                break
            loc += t
        out.append((f[0][5:21].strip(), loc, f))
    return gb_lines[:k + 1], out


def shift(lines, loc, delta):
    """The feature's lines with the numbers of its location moved by delta."""
    new = []
    in_loc = True
    for i, l in enumerate(lines):
        if i > 0 and l.strip().startswith("/"):
            in_loc = False
        if in_loc:
            head, body = (l[:21], l[21:])
            body = NUM.sub(lambda m: str(int(m.group()) + delta), body)
            l = head + body
        new.append(l)
    return new


def record(name, length, header, feats):
    lines = [f"LOCUS       {name}            {length} bp    DNA     circular   BCT 05-OCT-2026",
             f"DEFINITION  {name}, part of the example reference split for the multi-record test.",
             f"ACCESSION   {name}", f"VERSION     {name}.1", "FEATURES             Location/Qualifiers",
             f"     source          1..{length}"]
    for f in feats:
        lines += f
    lines += ["ORIGIN", "//"]
    return lines


def wrap(seq, width=70):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def main():
    fasta, gbk, out = sys.argv[1:4]
    target = int(sys.argv[4]) if len(sys.argv) > 4 else 2_000_000
    seq = read_fasta(fasta)
    with open(gbk) as fh:
        gb_lines = fh.read().split("\n")
    _, feats = features(gb_lines)
    spans = []
    for key, loc, lines in feats:
        if key == "source":
            continue
        nums = [int(n) for n in NUM.findall(loc)]
        if nums:
            spans.append((min(nums), max(nums), key, loc, lines))
    # The first gap no feature covers, at or after the target
    covered_to = 0
    split = None
    for lo, hi, *_ in sorted(spans):
        if lo > covered_to + 1 and covered_to + 1 >= target:
            split = covered_to + (lo - covered_to) // 2  # middle of the gap; partA ends here
            break
        covered_to = max(covered_to, hi)
    if split is None:
        sys.exit("no gap found")
    a_feats = [f[4] for f in spans if f[1] <= split]
    b_feats = [shift(f[4], f[3], -split) for f in spans if f[0] > split]
    la, lb = split, len(seq) - split
    os.makedirs(out, exist_ok=True)
    rec_a = record("partA", la, None, a_feats)
    rec_b = record("partB", lb, None, b_feats)
    for name, order in (("split", [("partA", seq[:split], rec_a), ("partB", seq[split:], rec_b)]),
                        ("swapped", [("partB", seq[split:], rec_b), ("partA", seq[:split], rec_a)])):
        with open(os.path.join(out, name + ".fasta"), "w") as fh:
            fh.write("".join(f">{n}.1 example split\n{wrap(s)}\n" for n, s, _ in order))
        with open(os.path.join(out, name + ".gb"), "w") as fh:
            fh.write("\n".join(l for _, _, r in order for l in r) + "\n")
    json.dump({"split": split, "length_A": la, "length_B": lb, "features_A": len(a_feats),
               "features_B": len(b_feats)}, open(os.path.join(out, "split.json"), "w"), indent=1)
    print(json.dumps({"split": split, "length_A": la, "length_B": lb}))


if __name__ == "__main__":
    main()
