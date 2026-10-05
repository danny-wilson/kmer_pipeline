#!/usr/bin/env python3
"""rungemma.py: run GEMMA on a subset of the unique patterns.
Port of rungemma.Rscript. Daniel Wilson (2022)."""
import argparse
import math
import os
import sys
import time

import rcompat
from rcompat import r_cat, r_paste0, r_stop


def gemma_phenotype_text(v):
    """A phenotype as written to GEMMA's phenotype file (D3): finite values at full precision,
    anything else (missing, NaN, infinite) as NA. R wrote 7 significant digits."""
    if v is None or not math.isfinite(v):
        return "NA"
    return "%.17g" % v


def main():
    rcompat.script_setup(__file__)
    start_time = time.monotonic()
    parser = argparse.ArgumentParser(description="rungemma.py run gemma on a subset of the unique patterns. "
                                                 "Daniel Wilson (2022)", allow_abbrev=False)
    parser.add_argument("--task-id", required=True)
    parser.add_argument("--p", required=True, help="number of parallel tasks")
    parser.add_argument("--kmerfile-prefix", required=True, help="prefix of the pattern and kinship files")
    parser.add_argument("--id-file", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--analysis-dir", required=True)
    parser.add_argument("--kmertype", required=True, help="protein or nucleotide")
    parser.add_argument("--kmer-length", required=True)
    parser.add_argument("--software-file", required=True)
    parser.add_argument("--covariate-file", default=None, help="GEMMA covariate file (first column all 1s)")
    parser.add_argument("--prepared", action="store_true",
                        help="use the phenotype file and decompressed kinship matrix prepare_gemma.py wrote once "
                             "for all tasks (the workflow does), instead of making them in each task")
    args = parser.parse_args()

    # Initialize variables
    t = rcompat.r_as_integer(args.task_id)
    p = rcompat.r_as_integer(args.p)
    prefix = args.kmerfile_prefix
    id_file = args.id_file
    output_prefix = args.output_prefix
    output_dir = args.analysis_dir
    kmertype = args.kmertype.lower()
    kmerlen = rcompat.r_as_integer(args.kmer_length)
    software_file = args.software_file
    covariate_file = args.covariate_file

    if p is None:
        r_stop("Error: p must be an integer", "\n")
    if not os.path.exists(id_file):
        r_stop("Error: sample ID file doesn't exist", "\n")
    if not os.path.exists(output_dir):
        r_stop("Error: output directory doesn't exist", "\n")
    if not output_dir.endswith("/"):
        output_dir = output_dir + "/"
    if kmertype != "protein" and kmertype != "nucleotide":
        r_stop("Error: kmer type must be either 'protein' or 'nucleotide'", "\n")
    if kmerlen is None:
        r_stop("Error: kmer length must be an integer", "\n")
    if not os.path.exists(software_file):
        r_stop("Error: software file doesn't exist", "\n")
    if covariate_file is not None and not os.path.exists(covariate_file):
        r_stop("Error: covariate file doesn't exist", "\n")
    if args.prepared and covariate_file is not None:
        # prepare_gemma.py wrote the covariates in GEMMA's format and id_file's order (the file
        # given may have an id column instead, N5)
        covariate_file = r_paste0(output_dir, kmertype, "kmer", kmerlen, "_gemma/", output_prefix, "_", kmertype, kmerlen,
                                  "_gemma_covariates.txt")
        if not os.path.exists(covariate_file):
            r_stop("Error: covariate file from prepare_gemma.py doesn't exist: ", covariate_file, "\n")

    keyfile = prefix + ".patternmerge.patternKey.txt.gz"
    keySizefile = prefix + ".patternmerge.patternKeySize.txt"
    kinfile = prefix + ".kinshipmerge.kinship.txt.gz"

    if not os.path.exists(keyfile):
        r_stop("Error: pattern key file doesn't exist: ", keyfile, " \n")
    if not os.path.exists(keySizefile):
        r_stop("Error: pattern key size file doesn't exist: ", keySizefile, " \n")
    if not os.path.exists(kinfile):
        r_stop("Error: kinship file doesn't exist: ", kinfile, " \n")

    # Sanity check covariate file
    if covariate_file is not None:
        nsamples = [x for x in rcompat.r_system_intern("wc -l " + id_file)[0].split(" ") if x != ""]
        nsamples = float(nsamples[0]) - 1
        covariates = rcompat.r_read_table(covariate_file, header=False, sep="\t")
        if not all(rcompat.r_as_numeric_value(v) == 1 for v in covariates.iloc[:, 0]):
            r_stop("Error: first column of covariate file must be a column of 1s for the intercept", "\n")
        if nsamples != len(covariates):
            r_stop("Error: covariate file is not the same length as the number of samples in id_file", "\n")

    # Read in software file
    software_paths = rcompat.r_read_table(software_file, header=True, sep="\t", quote="")
    names = [rcompat.r_as_character(v) for v in software_paths["name"]]
    paths = [rcompat.r_as_character(v) for v in software_paths["path"]]
    # Required software and script paths
    required_software = ["gemma"]
    if any(r not in names for r in required_software):
        r_stop("Error: missing required software path in the software file - requires " + ", ".join(required_software), "\n")
    gemmapath = [pth for nm, pth in zip(names, paths) if nm.lower() == "gemma"][0]
    if not os.path.exists(gemmapath):
        r_stop("Error: GEMMA path specified in the software paths file doesn't exist", "\n")
    if any(nm.lower() == "gemma_libraries" for nm in names):
        gemma_libraries_path = [pth for nm, pth in zip(names, paths) if nm.lower() == "gemma_libraries"][0]
        if not os.path.isdir(gemma_libraries_path):
            r_stop("Error: directory for GEMMA libraries does not exist", "\n")
        # N3: keep the existing library path (R's original dropped the $)
        gemmapath = "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:" + gemma_libraries_path + "; " + gemmapath

    # Report variables
    r_cat("#############################################", "\n")
    r_cat("Running on host: ", rcompat.r_system_intern("hostname"), "\n")
    r_cat("Command line arguments", "\n")
    r_cat(sys.argv[1:], "\n\n")
    r_cat("Parameters:", "\n")
    r_cat("task_id:", t, "\n")
    r_cat("p:", p, "\n")
    r_cat("Kmer file prefix:", prefix, "\n")
    r_cat("ID file path:", id_file, "\n")
    r_cat("Output prefix:", output_prefix, "\n")
    r_cat("Analysis directory:", output_dir, "\n")
    r_cat("Kmer type:", kmertype, "\n")
    r_cat("Kmer length:", kmerlen, "\n")
    r_cat("Software file:", software_file, "\n")
    r_cat("GEMMA path:", gemmapath, "\n")
    r_cat("Covariate file:", [] if covariate_file is None else covariate_file, "\n")
    r_cat("Process:", t, "\n")
    r_cat("\n")
    r_cat("Kmer pattern key file:", keyfile, "\n")
    r_cat("Kmer pattern key size file:", keySizefile, "\n")
    r_cat("Kmer kinship file:", kinfile, "\n")
    r_cat("#############################################", "\n\n")

    # Create gemma directory if doesn't already exist
    gemma_dir = output_dir + "/" + r_paste0(kmertype, "kmer", kmerlen, "_gemma", "/")  # file.path
    if not os.path.isdir(gemma_dir):
        rcompat.r_dir_create(gemma_dir)

    # Number of unique patterns: scan(what = integer(0)), not quiet
    toks = open(keySizefile).read().split()
    print(f"Read {len(toks)} item{'' if len(toks) == 1 else 's'}", file=sys.stderr, flush=True)
    if len(toks) != 1:
        r_stop("Could not read ", keySizefile)
    try:
        n = int(toks[0])
    except ValueError:
        r_stop("scan() expected 'an integer', got '", toks[0], "'")
    if n < 1:
        r_stop("Read n<1 from ", keySizefile)
    # Read in pheno file
    id_table = rcompat.r_read_table(id_file, header=True, sep="\t")
    pheno = [rcompat.r_as_numeric_value(v) for v in id_table["pheno"]]
    if args.prepared:  # N6: written once by prepare_gemma.py
        phenofile = r_paste0(gemma_dir, output_prefix, "_", kmertype, kmerlen, "_gemma_formatted_phenotype.txt")
        if not os.path.exists(phenofile):
            r_stop("Error: phenotype file from prepare_gemma.py doesn't exist: ", phenofile, "\n")
    else:
        phenofile = r_paste0(gemma_dir, output_prefix, "_", kmertype, kmerlen, "_gemma_formatted_phenotype_process", t,
                             ".txt")
        r_cat("Writing phenotype to gemma formatted file:", phenofile, "\n")
        rcompat.r_cat_lines([gemma_phenotype_text(v) for v in pheno], phenofile)

    # Compute other variables
    b = math.ceil(n / p)
    if p != math.ceil(n / b):
        r_cat("Warning: adjusting number of processes to equal ceiling(n/b)\n")
        p = math.ceil(n / b)

    if t > p:
        r_cat("Task", t, "not required\n")
        return

    # Create a temporary bimbam file from the unique patterns
    beg = b * (t - 1) + 1  # integers, written as plain integers (D1a)
    end = min(b * t, n)
    if end < beg:
        r_stop("Problem with input arguments, please check")
    genofile_prefix = r_paste0(gemma_dir, output_prefix, "_", kmertype, kmerlen, ".", beg, "-", end, ".prefix.bimbam.txt")
    genofile = r_paste0(gemma_dir, output_prefix, "_", kmertype, kmerlen, ".", beg, "-", end, ".bimbam.txt")
    # No directory for output file as gemma will always put files in subdirectory 'output'
    outfile_prefix = r_paste0(output_prefix, "_", kmertype, kmerlen, ".", beg, "-", end, "")
    # Pattern index (GEMMA's rs), then the two allele columns; indices as plain integers (D1a)
    with open(genofile_prefix, "w") as f:
        f.write("".join("%d\t1\t0\n" % k for k in rcompat.r_colon(beg, end)))
    rcompat.r_system(r_paste0("zcat ", keyfile, " | head -n ", end, " | tail -n ", end - beg + 1,
                              " | sed 's/./&\t/g' | paste ", genofile_prefix, " /dev/stdin > ", genofile))

    # If it does not exist, create output directory with correct permissions
    if not os.path.isdir(gemma_dir + "/output"):
        rcompat.r_dir_create(gemma_dir + "/output")
    # Move to directory above output directory as gemma will always put file in subdirectory output
    os.chdir(gemma_dir)

    # Run gemma
    if args.prepared:  # N6: decompressed once by prepare_gemma.py, removed after the last task
        kinfile_txt = r_paste0(gemma_dir, output_prefix, "_", kmertype, kmerlen, ".kinship.txt")
        if not os.path.exists(kinfile_txt):
            r_stop("Error: decompressed kinship matrix from prepare_gemma.py doesn't exist: ", kinfile_txt, "\n")
    else:
        kinfile_txt = r_paste0(gemma_dir, output_prefix, "_", kmertype, kmerlen, ".kinship.", beg, "-", end, ".txt")
        if kinfile_txt != kinfile:
            rcompat.r_system("zcat " + kinfile + " > " + kinfile_txt)
    if covariate_file is None:
        rcompat.r_system(gemmapath + " -g " + genofile + " -p " + phenofile + " -k " + kinfile_txt
                         + " -lmm 4 -maf 0 -o " + outfile_prefix)
    else:
        rcompat.r_system(gemmapath + " -g " + genofile + " -p " + phenofile + " -k " + kinfile_txt + " -c "
                         + covariate_file + " -lmm 4 -maf 0 -o " + outfile_prefix)

    # Delete temporary files
    r_cat("Deleting intermediate files", "\n")
    temporary = [genofile_prefix, genofile] + ([] if args.prepared else [phenofile] + ([kinfile_txt] if kinfile_txt != kinfile else []))
    for f in temporary:
        if os.path.lexists(f):
            os.remove(f)

    # The LRT p-value (column 10) of each tested pattern, labelled by its index (rs, column 2;
    # N2: GEMMA leaves out patterns that don't vary among the analysed genomes)
    assoc_file = gemma_dir + "output/" + outfile_prefix + ".assoc.txt"
    log_file = gemma_dir + "output/" + outfile_prefix + ".log.txt"
    pval_file = gemma_dir + "output/" + outfile_prefix + ".pval.txt.gz"
    rcompat.r_system("cut -f2,10 " + assoc_file + " | gzip -c > " + pval_file)
    rcompat.r_system("gzip " + assoc_file)
    rcompat.r_system("gzip " + log_file)

    r_cat("Finished in", (time.monotonic() - start_time) / 60, "minutes\n")


if __name__ == "__main__":
    main()
