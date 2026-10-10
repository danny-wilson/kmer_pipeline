#!/usr/bin/env python3
"""Split-reference test: compare the genes and regions each k-mer was
placed on with the unsplit example reference and with it split into two records (in order or
swapped). Regions are compared by their positions mapped back to the unsplit reference, not by
name. A k-mer may differ only where the split must change things: regions within WINDOW bases
of the split point or of the origin, and the wrap-round regions (whose extents change by
construction). Every other difference is a failure.

Run with the image's Python and PYTHONPATH set to a staging directory (for reference.py).
Usage: split_compare.py UNSPLIT_STAGE7 UNSPLIT_GB SPLIT_STAGE7 SPLIT_GB split.json split|swapped"""
import gzip
import glob
import json
import os
import sys

import reference

WINDOW = 5000


def placements(stage7):
    lookup_file = glob.glob(os.path.join(stage7, "nucleotidekmer31_kmergenealign", "*_gene_id_name_lookup.txt"))[0]
    names = {}
    for line in open(lookup_file):
        n, i = line.rstrip("\n").split("\t")
        names[i] = n
    merged = glob.glob(os.path.join(stage7, "*.kmeralignmerge.txt.gz"))[0]
    out = {}
    with gzip.open(merged, "rt") as fh:
        for line in fh:
            k, g = line.strip().split(",")
            out.setdefault(int(k), set()).add(names[g])
    return out


def region_ranges(gb):
    ref = reference.genes(gb)
    regs = reference.regions(ref, reference.records(gb))
    return {r.name: (r.kind, r.ranges) for r in regs}


def main():
    u_dir, u_gb, s_dir, s_gb, split_json, mode = sys.argv[1:7]
    info = json.load(open(split_json))
    la, lb, total = info["length_A"], info["length_B"], info["length_A"] + info["length_B"]
    split = info["split"]

    def to_unsplit(p):
        if mode == "split":
            return p
        return p + la if p <= lb else p - lb  # swapped: partB first

    u_regions, s_regions = region_ranges(u_gb), region_ranges(s_gb)

    def canon(name, regions, mapping):
        kind, ranges = regions[name]
        return (kind if kind == "wrap" else "region",
                tuple(sorted((min(mapping(a), mapping(b)), max(mapping(a), mapping(b))) for a, b in ranges)))

    def near(c):
        kind, ranges = c
        if kind == "wrap":
            return True
        return any(a - WINDOW <= split + 1 and b + WINDOW >= split or a <= WINDOW or b >= total - WINDOW
                   for a, b in ranges)

    u, s = placements(u_dir), placements(s_dir)
    stats = {"kmers": len(set(u) | set(s)), "same": 0, "differ_near_split_or_wrap": 0, "differ_elsewhere": 0}
    examples = []
    for k in sorted(set(u) | set(s)):
        cu = {canon(n, u_regions, lambda p: p) for n in u.get(k, ())}
        cs = {canon(n, s_regions, to_unsplit) for n in s.get(k, ())}
        if cu == cs:
            stats["same"] += 1
        elif all(near(c) for c in cu ^ cs):
            stats["differ_near_split_or_wrap"] += 1
        else:
            stats["differ_elsewhere"] += 1
            if len(examples) < 10:
                examples.append((k, sorted(u.get(k, ())), sorted(s.get(k, ()))))
    print(json.dumps(stats))
    for e in examples:
        print("DIFFERS:", e)
    sys.exit(0 if stats["differ_elsewhere"] == 0 else 1)


if __name__ == "__main__":
    main()
