#!/usr/bin/env python3
"""prepare_gemma.py: step 4's preparation, run once before the GEMMA tasks (and, with --cleanup,
once after them).

- The analysed set: the genomes GEMMA analyses, those with a finite phenotype and a complete
  covariate row. Written to <gemma_dir>/<prefix>.analysed_phenotypes.txt (id, phenotype,
  analysed); every later use of the phenotype (presence counts, MAF/MAC denominator, figures,
  reports) reads it, so they all describe the same genomes as GEMMA.
- GEMMA's phenotype file, written once: analysed phenotypes at full precision, NA otherwise.
- The presence count of each pattern among the analysed genomes
  (<prefix>.patternmerge.presenceCount.txt.gz, made by step 3 before).
- The kinship matrix, decompressed once for all GEMMA tasks (removed by --cleanup).

Stops with a clear message if the analysed genomes cannot be fitted: too few genomes for the
covariates, fewer than two distinct phenotype values, or covariates that are linearly
dependent among the analysed genomes (GEMMA would abort or return NaN for every pattern)."""
import argparse
import gzip
import math
import os
import shutil
import sys

import numpy as np

import rcompat
from rcompat import r_cat, r_stop

ANALYSED_SUFFIX = ".analysed_phenotypes.txt"


def gemma_dir(analysis_dir, kmer_type, kmer_length):
    return os.path.join(analysis_dir, f"{kmer_type}kmer{kmer_length}_gemma")


def file_prefix(output_prefix, kmer_type, kmer_length):
    return f"{output_prefix}_{kmer_type}{kmer_length}"


def analysed_file(analysis_dir, output_prefix, kmer_type, kmer_length):
    return os.path.join(gemma_dir(analysis_dir, kmer_type, kmer_length),
                        file_prefix(output_prefix, kmer_type, kmer_length) + ANALYSED_SUFFIX)


def read_phenotypes(id_file):
    """(ids as the scripts read them, phenotypes as numbers or None)."""
    table = rcompat.r_read_table(id_file, header=True, sep="\t")
    ids = [rcompat.r_as_character(v) for v in table["id"]]
    pheno = [rcompat.r_as_numeric_value(v) for v in table["pheno"]]
    return ids, pheno


def read_covariates(covariate_file):
    """Rows of numbers (None for missing), as rungemma reads the file; [] without a file."""
    if not covariate_file:
        return []
    table = rcompat.r_read_table(covariate_file, header=False, sep="\t")
    return [[rcompat.r_as_numeric_value(v) for v in row] for row in table.itertuples(index=False, name=None)]


def finite(v):
    return v is not None and isinstance(v, (int, float)) and math.isfinite(v)


def analysed_set(pheno, covariates):
    """For each genome, whether GEMMA analyses it: a finite phenotype and, if there are
    covariates, a complete row of finite covariates."""
    out = []
    for k, v in enumerate(pheno):
        ok = finite(v)
        if covariates:
            ok = ok and k < len(covariates) and all(finite(c) for c in covariates[k])
        out.append(ok)
    return out


def check_analysed(pheno, analysed, covariates):
    """Errors that make the model impossible to fit on the analysed genomes."""
    errors = []
    values = [v for v, a in zip(pheno, analysed) if a]
    ncov = len(covariates[0]) if covariates else 1  # GEMMA's intercept
    if len(values) <= ncov + 1:
        errors.append(f"only {len(values)} genomes have a phenotype (and complete covariates): GEMMA needs more "
                      f"than {ncov + 1} with {ncov} covariate column(s)")
    elif len(set(values)) < 2:
        errors.append(f"all {len(values)} analysed genomes have the same phenotype ({values[0]:g})")
    if covariates and len(values) > ncov:
        x = np.array([row for row, a in zip(covariates, analysed) if a], dtype=float)
        rank = np.linalg.matrix_rank(x)
        if rank < x.shape[1]:
            errors.append(f"the covariates are linearly dependent among the {len(values)} analysed genomes (rank "
                          f"{rank} of {x.shape[1]} columns, including the intercept): for example a covariate that "
                          "takes one value for all of them, or a set of indicator columns plus the intercept")
    return errors


def write_analysed(path, ids, pheno, analysed):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        fh.write("id\tpheno\tanalysed\n")
        for i, v, a in zip(ids, pheno, analysed):
            fh.write(f"{i}\t{'NA' if v is None else repr(float(v))}\t{'TRUE' if a else 'FALSE'}\n")
    os.replace(tmp, path)


def read_analysed(path):
    """(ids, phenotypes of the analysed genomes with None for the others)."""
    ids, pheno = [], []
    with open(path) as fh:
        next(fh)
        for line in fh:
            i, v, a = line.rstrip("\n").split("\t")
            ids.append(i)
            pheno.append(float(v) if a == "TRUE" else None)
    return ids, pheno


def gemma_phenotype_text(v):
    """A phenotype as written to GEMMA's phenotype file (D3): finite values at full precision,
    anything else (missing, NaN, infinite) as NA."""
    if not finite(v):
        return "NA"
    return "%.17g" % v


def analysed_phenotypes(analysis_dir, output_prefix, kmer_type, kmer_length, id_file, covariate_file, gemma_logs):
    """The analysed phenotypes for steps 6-7: from the analysed-phenotype file, or, for an
    analysis made by an earlier release (no such file), rebuilt from id_file and the covariates
    and checked against the number of analysed genomes in GEMMA's logs."""
    path = analysed_file(analysis_dir, output_prefix, kmer_type, kmer_length)
    if os.path.exists(path):
        return read_analysed(path)[1]
    ids, pheno = read_phenotypes(id_file)
    analysed = analysed_set(pheno, read_covariates(covariate_file))
    n = sum(analysed)
    for log in gemma_logs:
        logged = gemma_log_individuals(log)
        if logged is not None and logged != n:
            r_stop("Error: GEMMA analysed ", logged, " genomes (", log, ") but ", n, " have a phenotype and complete "
                   "covariates in id_file and the covariate file: these are not the inputs step 4 was run with; "
                   "rerun step 4 (overwrite = true)", "\n")
    r_cat("Warning: no ", os.path.basename(path), " (analysis made by an earlier release): analysed genomes taken "
          "from id_file and the covariate file, consistent with the GEMMA logs", "\n")
    return [v if a else None for v, a in zip(pheno, analysed)]


def gemma_log_individuals(log):
    """GEMMA's "number of analyzed individuals" in a (gzipped) log, or None."""
    opener = gzip.open if log.endswith(".gz") else open
    with opener(log, "rt") as fh:
        for line in fh:
            if "number of analyzed individuals" in line:
                return int(line.split("=")[-1])
    return None


def main():
    rcompat.script_setup(__file__)
    parser = argparse.ArgumentParser(description="prepare_gemma.py prepare step 4 (GEMMA) once for all its tasks",
                                     allow_abbrev=False)
    parser.add_argument("--kmerfile-prefix", required=True, help="prefix of the pattern and kinship files")
    parser.add_argument("--id-file", required=True)
    parser.add_argument("--covariate-file", default=None, help="GEMMA covariate file (first column all 1s)")
    parser.add_argument("--analysis-dir", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--kmer-type", required=True)
    parser.add_argument("--kmer-length", required=True)
    parser.add_argument("--cleanup", action="store_true", help="remove the decompressed kinship matrix")
    args = parser.parse_args()

    kmer_type = args.kmer_type.lower()
    gdir = gemma_dir(args.analysis_dir, kmer_type, args.kmer_length)
    prefix = file_prefix(args.output_prefix, kmer_type, args.kmer_length)
    kinship_txt = os.path.join(gdir, prefix + ".kinship.txt")
    if args.cleanup:
        if os.path.lexists(kinship_txt):
            os.remove(kinship_txt)
        r_cat("Removed", kinship_txt, "\n")
        return

    os.makedirs(gdir, exist_ok=True)
    ids, pheno = read_phenotypes(args.id_file)
    covariates = read_covariates(args.covariate_file)
    if covariates and len(covariates) != len(pheno):
        r_stop("Error: covariate file has ", len(covariates), " rows but id_file has ", len(pheno), " genomes", "\n")
    analysed = analysed_set(pheno, covariates)
    r_cat("Genomes analysed:", sum(analysed), "of", len(pheno), "\n")
    for k, (v, a) in enumerate(zip(pheno, analysed)):
        if finite(v) and not a:
            r_cat("Warning: genome", ids[k], "has a missing covariate, so it is not analysed", "\n")
    errors = check_analysed(pheno, analysed, covariates)
    if errors:
        r_stop("Error: ", "; ".join(errors), "\n")

    write_analysed(analysed_file(args.analysis_dir, args.output_prefix, kmer_type, args.kmer_length),
                   ids, pheno, analysed)
    rcompat.r_cat_lines([gemma_phenotype_text(v) if a else "NA" for v, a in zip(pheno, analysed)],
                        os.path.join(gdir, prefix + "_gemma_formatted_phenotype.txt"))

    import pattern2presencecount
    pattern2presencecount.write_presence_counts(args.kmerfile_prefix,
                                                os.path.join(args.analysis_dir, prefix),
                                                [k for k, a in enumerate(analysed) if a])

    kinship_gz = args.kmerfile_prefix + ".kinshipmerge.kinship.txt.gz"
    with gzip.open(kinship_gz, "rb") as src, open(kinship_txt + ".tmp", "wb") as dst:
        shutil.copyfileobj(src, dst, 1 << 22)
    os.replace(kinship_txt + ".tmp", kinship_txt)
    r_cat("Decompressed kinship matrix:", kinship_txt, "\n")


if __name__ == "__main__":
    sys.exit(main())
