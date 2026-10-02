#!/usr/bin/env python3
"""nucleotidekmermerge.py: merge the nucleotide k-mer files of all samples into a
sorted list of unique k-mers, as a pyramid of merges shared between tasks.
Port of nucleotidekmermerge.Rscript. Tasks wait for the files of the previous
round by polling for *.completed.txt files."""
import argparse
import math
import os
import re
import time

import rcompat
from rcompat import r_cat, r_paste0, r_stop

###################################################################################################
## Functions and software paths
###################################################################################################


def create_final_file(outfile, output_dir, output_prefix, kmer_length):
    # Check outfile isn't empty
    outfile_size = rcompat.r_system_intern("ls -l " + outfile + " | cut -d ' ' -f5")
    if outfile_size == ["0"]:
        r_stop(outfile, " file is empty", "\n")
    cmd = "mv " + outfile + " " + r_paste0(output_dir, output_prefix, "_nucleotide", kmer_length, ".kmermerge.txt")
    rcompat.r_system2("/bin/bash", "-c '" + cmd + "'")
    cmd = "gzip " + r_paste0(output_dir, output_prefix, "_nucleotide", kmer_length, ".kmermerge.txt")
    rcompat.r_system2("/bin/bash", "-c '" + cmd + "'")

    # Remove all completed files
    completed_files = rcompat.r_dir(output_dir, glob=r_paste0(output_prefix, ".nucleotide", kmer_length, "*.completed.txt"),
                                    full_names=True)
    cmd = " ".join(["rm", " ".join(completed_files)])
    rcompat.r_system2("/bin/bash", "-c '" + cmd + "'")


###################################################################################################


def main():
    rcompat.script_setup(__file__)
    start_time = time.monotonic()
    parser = argparse.ArgumentParser(description="nucleotidekmermerge.py merge kmer files. Daniel Wilson (2018)",
                                     allow_abbrev=False)
    parser.add_argument("--n", required=True, help="number of samples")
    parser.add_argument("--p", required=True, help="number of parallel tasks")
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--input-dir", required=True, help="directory of the per-sample k-mer files")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--id-file", required=True)
    parser.add_argument("--kmer-length", required=True)
    parser.add_argument("--process", required=True, help="task number (from 1)")
    args = parser.parse_args()

    # Initialize variables
    n = rcompat.r_as_integer(args.n)
    p = rcompat.r_as_integer(args.p)
    output_prefix = args.output_prefix
    input_dir = args.input_dir
    output_dir = args.output_dir
    id_file = args.id_file
    kmer_length = rcompat.r_as_numeric(args.kmer_length)
    process = rcompat.r_as_integer(args.process)

    if n is None:
        r_stop("Error: n must be an integer", "\n")
    if p is None:
        r_stop("Error: p must be an integer", "\n")
    if not os.path.exists(input_dir):
        r_stop("Error: input directory doesn't exist", "\n")
    if not input_dir.endswith("/"):
        input_dir = input_dir + "/"
    if not os.path.exists(output_dir):
        r_stop("Error: output directory doesn't exist", "\n")
    if not output_dir.endswith("/"):
        output_dir = output_dir + "/"
    if not os.path.exists(id_file):
        r_stop("Error: sample ID file doesn't exist", "\n")
    if process is None:
        r_stop("Error: process must be an integer", "\n")

    # Read in sample IDs
    id_table = rcompat.r_read_table(id_file, header=True, sep="\t")
    sample_id = [rcompat.r_as_character(v) for v in id_table["id"]]
    kmerpaths = [r_paste0(input_dir, s, ".kmer", kmer_length, ".txt.gz") for s in sample_id]

    b = math.ceil(n / p)
    if b == 1:
        r_stop("Cannot have batchsize = 1. Try p < n/2")
    p = math.ceil(n / b)
    t = process
    imax = math.ceil(math.log(n) / math.log(b))

    if t > p:
        r_cat("Task", t, "not required\n")
        return

    # Merge
    i = 0.0  # a double in R, so b^(i-1) and the file-name pieces follow R's types
    outfile = None
    while True:
        i = i + 1
        if not ((t % b ** int(i - 1)) == 0 or (t == p and i <= imax)):
            break
        r_cat("t=", t, "i=", i, "\n")
        if i == 1:
            # First round: merge source files
            outfile = r_paste0(output_dir, output_prefix, ".nucleotide", kmer_length, ".j.", i, ".", t, ".txt")
            r_cat("Creating temp file:", outfile, "\n")
            outfile_completed = r_paste0(output_dir, output_prefix, ".nucleotide", kmer_length, ".j.", i, ".", t,
                                         ".completed.txt")
            beg = float(b * (t - 1) + 1)
            end = min(b * t, n)
            r_cat("Beg:", beg, "End:", end, "\n")
            if end < beg:
                r_stop("Problem with input arguments, please check")
            infiles = rcompat.r_index(kmerpaths, rcompat.r_colon(int(beg), end))
            infiles = ["NA" if f is None else f for f in infiles]
            infiles_size = [float(rcompat.r_system_intern("zcat " + x + " | wc -c")[0]) for x in infiles]
            if not all(os.path.exists(f) for f in infiles) or not all(s > 0 for s in infiles_size):
                r_stop("Could not find files or files empty", " ".join(infiles))

            # Unzip the input files and remove the counts
            tmpinfiles = [re.sub(".txt.gz", r_paste0(".sorted.j.", i, ".", t, ".txt"), f) for f in infiles]
            subcmds = ["zcat " + f + " | cut -d \" \" -f1 > " + tf for f, tf in zip(infiles, tmpinfiles)]
            for subcmd in subcmds:
                rcompat.r_system2("/bin/bash", "-c '" + subcmd + "'")

            if len(infiles) == 1:
                cmd = "cat " + tmpinfiles[0] + " > " + outfile
            elif len(infiles) == 2:
                cmd = "LC_ALL=C sort -um " + tmpinfiles[0] + " " + tmpinfiles[1] + " > " + outfile
            else:
                cmd = ("LC_ALL=C sort -um " + tmpinfiles[0] + " " + tmpinfiles[1]
                       + "".join(" | LC_ALL=C sort -um - " + f for f in tmpinfiles[2:]) + " > " + outfile)
            rcompat.r_system2("/bin/bash", "-c '" + cmd + "'")
            rcompat.r_system2("/bin/bash", "-c 'touch " + outfile_completed + "'")

            # Remove the temporary input files
            for f in tmpinfiles:
                if os.path.lexists(f):
                    os.remove(f)
        else:
            # Subsequent rounds: merge merged files
            te = math.ceil(t / (b ** (i - 1))) * int(b ** (i - 1))
            outfile = r_paste0(output_dir, output_prefix, ".nucleotide", kmer_length, ".j.", i, ".", te, ".txt")
            r_cat("Creating temp file:", outfile, "\n")
            outfile_completed = r_paste0(output_dir, output_prefix, ".nucleotide", kmer_length, ".j.", i, ".", te,
                                         ".completed.txt")
            beg = te - float(b) ** (i - 1) + float(b) ** (i - 2)
            end = min(te, math.ceil(t / b ** (i - 2)) * int(b ** (i - 2)))
            r_cat(r_paste0("Beg", i, ":"), beg, r_paste0("End", i, ":"), end, "\n")
            if end < beg:
                r_stop("Problem with input arguments, please check")
            inc = float(b) ** (i - 2)
            js = rcompat.r_seq(beg, end, inc)
            infiles = [r_paste0(output_dir, output_prefix, ".nucleotide", kmer_length, ".j.", i - 1, ".", j, ".txt")
                       for j in js]
            infiles_completed = [r_paste0(output_dir, output_prefix, ".nucleotide", kmer_length, ".j.", i - 1, ".", j,
                                          ".completed.txt") for j in js]
            nattempts = 0
            while not all(os.path.exists(f) for f in infiles_completed) or not all(os.path.exists(f) for f in infiles):
                nattempts = nattempts + 1
                if nattempts > 100:
                    r_stop("Could not find files", "".join(infiles_completed))
                time.sleep(60)

            infiles_size = [float(rcompat.r_system_intern("ls -l " + x + " | cut -d ' ' -f5")[0]) for x in infiles]
            if any(s == 0 for s in infiles_size):
                r_stop("One or more file size is zero ", "".join(infiles), "\n")

            if len(infiles) == 1:
                cmd = "mv " + infiles[0] + " " + outfile
            elif len(infiles) == 2:
                cmd = "LC_ALL=C sort -um " + infiles[0] + " " + infiles[1] + " > " + outfile
            else:
                cmd = ("LC_ALL=C sort -um " + infiles[0] + " " + infiles[1]
                       + "".join(" | LC_ALL=C sort -um - " + f for f in infiles[2:]) + " > " + outfile)
            rcompat.r_system2("/bin/bash", "-c '" + cmd + "'")
            # If length of infiles is one, it has been moved and doesn't exist. If more than one, delete the temp files.
            if len(infiles) > 1:
                cmd = " ".join(["rm", " ".join(infiles)])
                rcompat.r_system2("/bin/bash", "-c '" + cmd + "'")
            rcompat.r_system2("/bin/bash", "-c 'touch " + outfile_completed + "'")
    if t == p:
        ## Copy final file to end location
        create_final_file(outfile=outfile, output_dir=output_dir, output_prefix=output_prefix, kmer_length=kmer_length)

    r_cat("Finished in", (time.monotonic() - start_time) / 60, "minutes\n")


if __name__ == "__main__":
    main()
