#!/usr/bin/env python3
"""Golden-output hashes: a committed TSV standing in for the full goldens
(which stay local only, for file-size reasons). Shares its manifest reader and
per-class normalization with compare.py, so a hash and a direct comparison can
never disagree about what counts as a difference.

Usage:
  golden_hashes.py generate --run NAME DIR [--run NAME DIR ...] [--manifest FILE] > golden_hashes.tsv
  golden_hashes.py verify TSV --run NAME DIR [--run NAME DIR ...] [--manifest FILE] [--visual-advisory]

generate hashes every classified (non-ignored) file under each named run
directory: `exact`/`strip`/`numeric` classes hash the same normalized content
compare.py would compare (gzip-decompressed, strip-lines-removed); `visual`
(PNG) hashes the raw bytes.

verify re-hashes the given run directories and compares against the TSV. It
fails if: a TSV row's path is missing from the run; a classified, non-ignored
file in the run is missing from the TSV; or a hash differs -- for `visual`
rows only, --visual-advisory downgrades a hash mismatch to a warning (a path
listed but entirely missing from the run still fails even then).
"""
import argparse
import datetime
import hashlib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from compare import classify, content, list_files, read_manifest, strip_lines  # noqa: E402

HEADER = "run\tpath\tclass\tsha256"


def hash_one(run_dir, rel, cls, regex):
    if cls == "visual":
        with open(os.path.join(run_dir, rel), "rb") as fh:
            data = fh.read()
    else:
        data = content(os.path.join(run_dir, rel))
        if cls == "strip":
            data = strip_lines(data, regex)
        elif cls == "numeric":
            return "-"
    return hashlib.sha256(data).hexdigest()


def hash_run(run_dir, rules):
    """{rel_path: (class, sha256_or_'-')} for every classified, non-ignored file."""
    out = {}
    for rel in sorted(list_files(run_dir)):
        rule = classify(rel, rules)
        if rule is None:
            sys.exit(f"UNCLASSIFIED {rel} (in {run_dir})")
        _, cls, regex = rule
        if cls == "ignore":
            continue
        out[rel] = (cls, hash_one(run_dir, rel, cls, regex))
    return out


def generate(runs, manifest_path, commit=None, nextflow_version=None, image_digest=None):
    rules = read_manifest(manifest_path)
    rows = []
    for name, run_dir in runs:
        for rel, (cls, digest) in hash_run(run_dir, rules).items():
            rows.append((name, rel, cls, digest))
    rows.sort()
    print("# golden_hashes.tsv -- sha256, class exact/strip/numeric/visual")
    print(f"# generated {datetime.datetime.now(datetime.timezone.utc).isoformat()}")
    print(f"# commit: {commit or '(not given)'}")
    print(f"# nextflow version: {nextflow_version or '(not given)'}")
    print(f"# image digest: {image_digest or '(not given)'}")
    print("# regenerate with: golden_hashes.py generate "
          + " ".join(f"--run {n} $KMER_E2E_ROOT/goldens/<label>-<sha7>/maxp2/{n}/stage7" for n, _ in runs)
          + " --commit <sha> --nextflow-version <version> --image-digest <digest> > tests/e2e/golden_hashes.tsv")
    print(HEADER)
    for row in rows:
        print("\t".join(row))


def read_tsv(path):
    rows = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line == HEADER:
                continue
            run, rel, cls, digest = line.split("\t")
            rows[(run, rel)] = (cls, digest)
    return rows


def verify(tsv_path, runs, manifest_path, visual_advisory):
    rules = read_manifest(manifest_path)
    want = read_tsv(tsv_path)
    problems, warnings = [], []
    seen = set()
    for name, run_dir in runs:
        got = hash_run(run_dir, rules)
        for rel, (cls, digest) in got.items():
            seen.add((name, rel))
            if (name, rel) not in want:
                problems.append(f"NOT IN TSV: {name} {rel} ({cls})")
                continue
            want_cls, want_digest = want[(name, rel)]
            if want_cls != cls:
                problems.append(f"CLASS CHANGED: {name} {rel}: tsv={want_cls} now={cls}")
            elif want_digest != digest:
                msg = f"HASH MISMATCH: {name} {rel} ({cls})"
                if cls == "visual" and visual_advisory:
                    warnings.append(msg)
                else:
                    problems.append(msg)
    for (name, rel) in want:
        if name in {n for n, _ in runs} and (name, rel) not in seen:
            problems.append(f"MISSING: {name} {rel} (listed in TSV, absent from the run)")
    return problems, warnings


def parse_runs(pairs):
    if len(pairs) % 2 != 0:
        sys.exit("--run takes two values: NAME DIR")
    return [(pairs[i], pairs[i + 1]) for i in range(0, len(pairs), 2)]


def main():
    ap = argparse.ArgumentParser(allow_abbrev=False)
    ap.add_argument("mode", choices=["generate", "verify"])
    ap.add_argument("tsv", nargs="?", help="required for verify")
    ap.add_argument("--run", action="append", nargs=2, metavar=("NAME", "DIR"), default=[])
    ap.add_argument("--manifest", default=os.path.join(HERE, "manifest.tsv"))
    ap.add_argument("--visual-advisory", action="store_true")
    ap.add_argument("--commit", help="generate: the commit the goldens were made from, for the TSV header")
    ap.add_argument("--nextflow-version", help="generate: for the TSV header")
    ap.add_argument("--image-digest", help="generate: the toolchain image's public registry digest, for the TSV header")
    args = ap.parse_args()
    runs = [(n, d) for n, d in args.run]
    if not runs:
        sys.exit("at least one --run NAME DIR is required")

    if args.mode == "generate":
        generate(runs, args.manifest, args.commit, args.nextflow_version, args.image_digest)
        return

    if not args.tsv:
        sys.exit("verify needs the TSV path")
    problems, warnings = verify(args.tsv, runs, args.manifest, args.visual_advisory)
    for w in warnings:
        print("WARNING:", w)
    for p in problems:
        print(p)
    print(f"{len(problems)} problems, {len(warnings)} warnings")
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
