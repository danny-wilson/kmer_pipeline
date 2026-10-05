#!/usr/bin/env python3
"""pattern2presencecount.py: count the genomes each k-mer pattern is present in
(optionally only genomes with a non-NA phenotype). Port of pattern2presencecount.Rscript."""
import argparse
import gzip
import math
import os

import numpy as np

import rcompat
from rcompat import r_cat, r_stop

_LOGICAL = {"FALSE": False, "F": False, "false": False, "False": False,
            "TRUE": True, "T": True, "true": True, "True": True}

CHUNK = 200000  # patterns read at a time (R reads them in up to 20 batches)


def r_seq_length_out(from_, to, length_out):
    """seq(from, to, length.out =) for doubles, as R's seq.default."""
    length_out = math.ceil(length_out)
    if length_out == 0:
        return []
    if length_out > 2:
        if from_ == to:
            return [from_] * length_out
        n1 = length_out - 1
        d = (to - from_) / n1
        return [from_] + [from_ + k * d for k in range(1, length_out - 1)] + [to]
    return [from_, to][:length_out]


def format_count(v):
    """A presence count: a plain integer (D1a; R's cat() wrote 100000 as "1e+05")."""
    if v == v and v == int(v):
        return str(int(v))
    return rcompat.r_str(v, 7)


def presence_counts(lines, notNA):
    """get_pc for each pattern: the sum over the selected genomes of the digits of
    the pattern string; NA if any selected character is not a digit."""
    out = np.empty(len(lines), dtype=float)
    widths = {len(l) for l in lines}
    if len(widths) == 1 and lines:
        w = widths.pop()
        a = np.frombuffer("".join(lines).encode("latin-1"), dtype=np.uint8).reshape(len(lines), w)
        idx = [i for i in notNA if i < w]
        sel = a[:, idx].astype(np.int16) - 48
        bad = ((sel < 0) | (sel > 9)).any(axis=1) | (len(idx) < len(notNA))
        out[:] = sel.sum(axis=1)
        out[bad] = np.nan
        return out
    for k, l in enumerate(lines):  # patterns of unequal length (not written by the pipeline)
        vals = [float(c) if c.isdigit() else math.nan for c in l]
        out[k] = sum(vals[i] if i < len(vals) else math.nan for i in notNA)
    return out


def write_presence_counts(kmerfilePrefix, out_prefix, selected):
    """Write <out_prefix>.patternmerge.presenceCount.txt.gz: for each pattern of <kmerfilePrefix>,
    the number of the selected genomes (0-based columns) it is present in. Returns its path."""
    kmerKeySizeFile = kmerfilePrefix + ".patternmerge.patternKeySize.txt"
    kmerKey = kmerfilePrefix + ".patternmerge.patternKey.txt.gz"
    pheno_notNA = selected
    # Read in total number of kmer patterns
    sizes = [float(v) for v in open(kmerKeySizeFile).read().split()]
    nPatterns = sizes[0] if len(sizes) == 1 else sizes
    r_cat("Number of patterns:", nPatterns, "\n")
    if len(sizes) != 1:
        r_stop("Error: expected one number in ", kmerKeySizeFile)

    s = [float(round(v)) for v in r_seq_length_out(1.0, nPatterns, min(nPatterns, 21))]
    if not s:
        r_stop("Error: no patterns in ", kmerKeySizeFile)

    # The patterns are read in the batches R uses (lines s[1] .. end); every
    # pattern is counted once, so they are read here in chunks from line s[1]
    presencecount = []
    with rcompat.r_open(kmerKey) as f:
        for _ in range(int(s[0]) - 1):
            next(f)
        chunk = []
        for line in f:
            tokens = line.split()
            chunk.extend(tokens)  # scan(what = character()) reads whitespace-separated items
            if len(chunk) >= CHUNK:
                presencecount.append(presence_counts(chunk, pheno_notNA))
                chunk = []
        if chunk:
            presencecount.append(presence_counts(chunk, pheno_notNA))
    presencecount = np.concatenate(presencecount) if presencecount else np.array([])
    if len(presencecount) != nPatterns:
        r_stop("Error: length of presence count vector not equal to total number of patterns", "\n")

    # Written to a temporary name, then renamed (gzip would refuse to replace an existing file)
    presencecountfile = out_prefix + ".patternmerge.presenceCount.txt.gz"
    with gzip.open(presencecountfile + ".tmp", "wt") as f:
        f.write("".join(("NA" if v != v else format_count(v)) + "\n" for v in presencecount)
                if len(presencecount) else "\n")
    os.replace(presencecountfile + ".tmp", presencecountfile)
    r_cat("Written presence counts to file:", presencecountfile, "\n")
    return presencecountfile


def main():
    rcompat.script_setup(__file__)
    parser = argparse.ArgumentParser(description="pattern2presencecount.py count the genomes each pattern is present in",
                                     allow_abbrev=False)
    parser.add_argument("--kmerfile-prefix", required=True, help="prefix of the patternmerge files")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--id-file", required=True)
    parser.add_argument("--include-na", default="TRUE", help="count genomes with an NA phenotype (TRUE/FALSE)")
    args = parser.parse_args()

    # Initialize variables
    kmerfilePrefix = args.kmerfile_prefix
    kmerKeySizeFile = kmerfilePrefix + ".patternmerge.patternKeySize.txt"
    kmerKey = kmerfilePrefix + ".patternmerge.patternKey.txt.gz"
    output_dir = args.output_dir
    id_file = args.id_file
    includeNA = _LOGICAL.get(args.include_na)

    if not os.path.exists(output_dir):
        r_stop("Error: output directory doesn't exist", "\n")
    if not output_dir.endswith("/"):
        output_dir = output_dir + "/"
    if not os.path.exists(kmerKeySizeFile):
        r_stop("Error: kmer pattern key size file doesn't exist", "\n")
    if not os.path.exists(kmerKey):
        r_stop("Error: kmer pattern file doesn't exist", "\n")
    if not os.path.exists(id_file):
        r_stop("Error: sample id file does not exist", "\n")
    if includeNA is None:
        r_stop("Error: include NAs must be a logical", "\n")

    # Remove any directories in kmerfilePrefix
    parts = kmerfilePrefix.split("/")
    while parts and parts[-1] == "":  # strsplit drops a trailing empty piece
        parts.pop()
    kmerfilePrefix_noDir = parts[-1]

    # Read in ID file and phenotype
    id_table = rcompat.r_read_table(id_file, header=True, sep="\t")
    pheno = [rcompat.r_as_numeric_value(v) for v in id_table["pheno"]]
    # Which genomes are not NA for this pheno
    # If includeNA is TRUE, set this to include all samples
    if not includeNA and any(v is None for v in pheno):
        pheno_notNA = [i for i, v in enumerate(pheno) if v is not None]
        r_cat("Calculating presence counts for samples with non NA phenotypes:",
              str(len(pheno_notNA)) + "/" + str(len(pheno)), "samples", "\n")
    else:
        pheno_notNA = list(range(len(pheno)))
        r_cat("Calculating presence counts for all samples", "\n")

    write_presence_counts(kmerfilePrefix, output_dir + kmerfilePrefix_noDir, pheno_notNA)


if __name__ == "__main__":
    main()
