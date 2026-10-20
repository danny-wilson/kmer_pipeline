#!/usr/bin/env python3
"""Compare two analysis directories according to the output manifest (manifest.tsv).

Usage: compare.py DIR_A DIR_B [--manifest manifest.tsv]

Each file's path relative to its directory is matched against the manifest
patterns in order; the first match decides its class:
  exact   identical content (after gzip decompression for .gz)
  path    identical after replacing a run's own base_dir with the literal
          string "@BASE_DIR@" on each side (--prefix-a/--prefix-b); a fixed,
          literal substring substitution, not a path-canonicalization pass.
          --extra-a/--extra-b (each "OLD=NEW") apply one further literal
          substring substitution per side, after the @BASE_DIR@ one -- for
          comparing two runs that also differ in a data-directory name below
          base_dir itself (e.g. the launcher's own clone layout vs the e2e
          harness's), not just in base_dir's own spelling (Part D, D3).
  strip   identical after deleting lines that match the pattern's regex
  numeric identical apart from numbers that differ by at most the relative
          tolerance in the regex column
  visual  must exist in both; content reviewed by eye
  ignore  not compared, and may exist on one side only (logs vary in number)
Unmatched files are an error, so new outputs can't slip past unclassified.
Exit status 0 only if everything matches.
"""
import argparse
import fnmatch
import gzip
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def read_manifest(path):
    rules = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            pattern, cls = fields[0], fields[1]
            if cls == "numeric":
                regex = float(fields[2])
            else:
                regex = re.compile(fields[2]) if len(fields) > 2 and fields[2] else None
            if cls not in ("exact", "path", "strip", "numeric", "visual", "ignore"):
                sys.exit(f"bad class {cls!r} in manifest line: {line}")
            rules.append((pattern, cls, regex))
    return rules


def list_files(root):
    out = set()
    for d, _, files in os.walk(root):
        for f in files:
            out.add(os.path.relpath(os.path.join(d, f), root))
    return out


def content(path):
    with open(path, "rb") as fh:
        data = fh.read()
    if path.endswith(".gz"):
        data = gzip.decompress(data)
    return data


def strip_lines(data, regex):
    return b"\n".join(l for l in data.split(b"\n") if not regex.search(l.decode("utf-8", "replace")))


NUM_RE = re.compile(rb"-?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")


def numeric_equal(a, b, rtol):
    """Equal apart from numbers within rtol (relative): the texts with numbers
    replaced must match, and each pair of numbers must agree."""
    if NUM_RE.sub(b"#", a) != NUM_RE.sub(b"#", b):
        return False
    for x, y in zip(NUM_RE.findall(a), NUM_RE.findall(b)):
        if x != y:
            fx, fy = float(x), float(y)
            if abs(fx - fy) > rtol * max(abs(fx), abs(fy)):
                return False
    return True


def first_difference(a, b):
    la, lb = a.split(b"\n"), b.split(b"\n")
    for i, (x, y) in enumerate(zip(la, lb)):
        if x != y:
            return f"line {i + 1}: {x[:120]!r} vs {y[:120]!r}"
    return f"line counts {len(la)} vs {len(lb)}"


def classify(rel, rules):
    return next((r for r in rules if fnmatch.fnmatch(rel, r[0])), None)


def normalize_path_class(data, prefix, extra=None):
    """Replace a literal base_dir prefix with "@BASE_DIR@", then optionally one
    further literal (old, new) substring substitution (extra). Both are fixed
    substring substitutions only -- no canonicalization, so an un-collapsed
    double slash in the golden (e.g. "kmergwas//") survives unchanged on both
    sides."""
    if prefix:
        data = data.replace(prefix.encode(), b"@BASE_DIR@")
    if extra:
        old, new = extra
        data = data.replace(old.encode(), new.encode())
    return data


def compare(dir_a, dir_b, manifest_path, prefix_a=None, prefix_b=None, extra_a=None, extra_b=None):
    rules = read_manifest(manifest_path)
    files_a, files_b = list_files(dir_a), list_files(dir_b)
    problems, counts = [], {}

    for rel in sorted(files_a | files_b):
        rule = classify(rel, rules)
        if rule is None:
            problems.append(f"UNCLASSIFIED {rel}")
            continue
        _, cls, regex = rule
        counts[cls] = counts.get(cls, 0) + 1
        if cls == "ignore":
            continue
        if rel not in files_a or rel not in files_b:
            problems.append(f"MISSING in {'A' if rel not in files_a else 'B'}: {rel}")
            continue
        if cls == "visual":
            continue
        a = content(os.path.join(dir_a, rel))
        b = content(os.path.join(dir_b, rel))
        if cls == "strip":
            a, b = strip_lines(a, regex), strip_lines(b, regex)
        if cls == "path":
            a = normalize_path_class(a, prefix_a, extra_a)
            b = normalize_path_class(b, prefix_b, extra_b)
        if cls == "numeric":
            if not numeric_equal(a, b, regex):
                problems.append(f"DIFFERS ({cls}) {rel}: {first_difference(a, b)}")
        elif a != b:
            problems.append(f"DIFFERS ({cls}) {rel}: {first_difference(a, b)}")

    return problems, counts, len(files_a | files_b)


def main():
    ap = argparse.ArgumentParser(allow_abbrev=False)
    ap.add_argument("dir_a")
    ap.add_argument("dir_b")
    ap.add_argument("--manifest", default=os.path.join(HERE, "manifest.tsv"))
    ap.add_argument("--prefix-a", default=None, help="dir_a's own base_dir, for the 'path' class")
    ap.add_argument("--prefix-b", default=None, help="dir_b's own base_dir, for the 'path' class")
    ap.add_argument("--extra-a", default=None, metavar="OLD=NEW",
                     help="one further literal substitution for dir_a's 'path' class files, after @BASE_DIR@")
    ap.add_argument("--extra-b", default=None, metavar="OLD=NEW",
                     help="one further literal substitution for dir_b's 'path' class files, after @BASE_DIR@")
    args = ap.parse_args()
    extra_a = tuple(args.extra_a.split("=", 1)) if args.extra_a else None
    extra_b = tuple(args.extra_b.split("=", 1)) if args.extra_b else None
    problems, counts, total = compare(args.dir_a, args.dir_b, args.manifest, args.prefix_a, args.prefix_b,
                                       extra_a, extra_b)

    for p in problems:
        print(p)
    summary = ", ".join(f"{k} {v}" for k, v in sorted(counts.items()))
    print(f"{total} files ({summary}); {len(problems)} problems")
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
