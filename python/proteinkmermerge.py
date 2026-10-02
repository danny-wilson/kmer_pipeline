#!/usr/bin/env python3
"""proteinkmermerge.py: merge the protein k-mer files of all samples into a sorted
list of unique k-mers, as a pyramid of merges shared between tasks.
Port of proteinkmermerge.Rscript. Tasks wait for the files of the previous round
by polling for *.completed.txt files."""
import argparse
import math
import os
import time

import rcompat
from rcompat import r_cat, r_paste0, r_stop

###################################################################################################
## Functions and software paths
###################################################################################################


def sort_final_file(output_dir, output_prefix, sort_strings, kmer_length):
    # Read in kmers
    kmerFile = r_paste0(output_dir, output_prefix, "_protein", kmer_length, ".kmermerge.unsorted.txt.gz")
    kmers = rcompat.r_system_intern("zcat " + kmerFile)
    kmer_output_file = r_paste0(output_dir, output_prefix, "_protein", kmer_length, ".kmermerge.wdummycount.txt")
    with open(kmer_output_file, "w") as f:  # write.table of cbind(kmers, 1)
        f.write("".join(k + "\t1\n" for k in kmers))
    # Gzip file
    rcompat.r_system("gzip " + kmer_output_file)
    # Run sort strings
    kmer_output_file_gz = kmer_output_file + ".gz"
    final_kmer_txt_gz = r_paste0(output_dir, output_prefix, "_protein", kmer_length, ".kmermerge.sorted.wdummycount.txt.gz")
    sortCommand = " ".join([sort_strings, kmer_output_file_gz, "| gzip -c >", final_kmer_txt_gz])
    rcompat.r_system(sortCommand)
    # Remove column of dummy counts
    final_kmer_txt_gz_sorted = r_paste0(output_dir, output_prefix, "_protein", kmer_length, ".kmermerge.txt.gz")
    rcompat.r_system("zcat " + final_kmer_txt_gz + " | cut -f1 | gzip -c > " + final_kmer_txt_gz_sorted)
    # Tests
    kmers_new = rcompat.r_system_intern("zcat " + final_kmer_txt_gz_sorted)
    if len(kmers) != len(kmers_new):
        r_stop("Error: issue when sorting the final kmer file, number of kmers does not match", "\n")
    if not set(kmers) <= set(kmers_new):
        r_stop("Error: issue when sorting the final kmer file, the kmers do not match", "\n")
    if kmers != kmers_new:
        r_stop("Error: issue when sorting the final kmer file, the kmers are not in the same order", "\n")
    # Remove temp dummy count files
    rcompat.r_system("rm " + kmer_output_file_gz)
    rcompat.r_system("rm " + final_kmer_txt_gz)
    # Remove the unsorted file
    rcompat.r_system("rm " + r_paste0(output_dir, output_prefix, "_protein", kmer_length, ".kmermerge.unsorted.txt.gz"))


def create_final_file(outfile, output_dir, output_prefix, kmer_length):
    # Check outfile isn't empty
    outfile_size = rcompat.r_system_intern("ls -l " + outfile + " | cut -d ' ' -f5")
    if outfile_size == ["0"]:
        r_stop(outfile, " file is empty", "\n")
    cmd = "mv " + outfile + " " + r_paste0(output_dir, output_prefix, "_protein", kmer_length, ".kmermerge.unsorted.txt")
    rcompat.r_system2("/bin/bash", "-c '" + cmd + "'")
    cmd = "gzip " + r_paste0(output_dir, output_prefix, "_protein", kmer_length, ".kmermerge.unsorted.txt")
    rcompat.r_system2("/bin/bash", "-c '" + cmd + "'")

    # Remove all completed files
    completed_files = rcompat.r_system_intern("ls " + r_paste0(output_dir, output_prefix, ".protein", kmer_length,
                                                                "*.completed.txt"))
    cmd = " ".join(["rm", " ".join(completed_files)])
    rcompat.r_system2("/bin/bash", "-c '" + cmd + "'")


def dedupe(outfile):
    """bash sort ignores *s, so keeps duplicate k-mers: keep the first of each
    (scan, unique, cat as R)."""
    kmers_i = rcompat.r_unique(rcompat.r_scan_lines(outfile, quiet=True))
    rcompat.r_cat_lines(kmers_i, outfile)


###################################################################################################


def main():
    rcompat.script_setup(__file__)
    start_time = time.monotonic()
    parser = argparse.ArgumentParser(
        description="proteinkmermerge.py merge protein kmer files. Daniel Wilson (2018) kmermerge.Rscript "
                    "modified to proteinkmermerge.Rscript by Sarah Earle (2019)", allow_abbrev=False)
    parser.add_argument("--n", required=True, help="number of samples")
    parser.add_argument("--p", required=True, help="number of parallel tasks")
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--input-dir", required=True, help="directory of the per-sample k-mer files")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--id-file", required=True)
    parser.add_argument("--kmer-length", required=True)
    parser.add_argument("--software-file", required=True)
    parser.add_argument("--process", required=True, help="task number (from 1)")
    args = parser.parse_args()

    # Initialize variables
    n = rcompat.r_as_integer(args.n)
    p = rcompat.r_as_integer(args.p)
    output_prefix = args.output_prefix
    input_dir = args.input_dir
    output_dir = args.output_dir
    id_file = args.id_file
    kmer_length = rcompat.r_as_integer(args.kmer_length)
    software_file = args.software_file
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
    if not os.path.exists(software_file):
        r_stop("Error: software file doesn't exist", "\n")
    if process is None:
        r_stop("Error: process must be an integer", "\n")

    # Read in software file
    software_paths = rcompat.r_read_table(software_file, header=True, sep="\t", quote="")
    names = [rcompat.r_as_character(v) for v in software_paths["name"]]
    paths = [rcompat.r_as_character(v) for v in software_paths["path"]]
    # Required software and script paths
    required_software = ["scriptpath"]
    if any(r not in names for r in required_software):
        r_stop("Error: missing required software path in the software file - requires " + ", ".join(required_software), "\n")
    script_location = [pth for nm, pth in zip(names, paths) if nm.lower() == "scriptpath"][0]
    if not os.path.isdir(script_location):
        r_stop("Error: script location directory specified in the software paths file doesn't exist", "\n")
    sort_strings = script_location + "/sort_strings"
    if not os.path.exists(sort_strings):
        r_stop("Error: sort_strings path doesn't exist - check pipeline script location in the software file", "\n")

    # Read in sample IDs
    id_table = rcompat.r_read_table(id_file, header=True, sep="\t")
    sample_id = [rcompat.r_as_character(v) for v in id_table["id"]]

    if n != len(sample_id):
        r_stop("n does not equal the number of samples in the sample ID file", "\n")

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
    i = 0.0  # a double in R
    outfile = None
    while True:
        i = i + 1
        if not ((t % b ** int(i - 1)) == 0 or (t == p and i <= imax)):
            break
        r_cat("t =", t, "i =", i, "\n")
        if i == 1:
            # First round: merge source files
            outfile = r_paste0(output_dir, output_prefix, ".protein", kmer_length, ".j.", i, ".", t, ".txt")
            r_cat("Creating temp file:", outfile, "\n")
            outfile_completed = r_paste0(output_dir, output_prefix, ".protein", kmer_length, ".j.", i, ".", t,
                                         ".completed.txt")
            beg = float(b * (t - 1) + 1)
            end = min(b * t, n)
            r_cat("Beg:", beg, " End:", end, "\n")
            if end < beg:
                r_stop("Problem with input arguments, please check")
            ids = ["NA" if s is None else s for s in rcompat.r_index(sample_id, rcompat.r_colon(int(beg), end))]
            infiles = [r_paste0(input_dir, s, ".kmer", kmer_length, ".txt.gz") for s in ids]
            infiles_size = [float(rcompat.r_system_intern("zcat " + x + " | wc -c")[0]) for x in infiles]
            if not all(os.path.exists(f) for f in infiles) or not all(s > 0 for s in infiles_size):
                r_stop("Could not find files or files empty ", " ".join(infiles))

            if len(infiles) == 1:
                cmd = "zcat " + infiles[0] + " | cut -f1 > " + outfile
            elif len(infiles) == 2:
                cmd = "LC_ALL=C sort -um <(zcat " + infiles[0] + " | cut -f1) <(zcat " + infiles[1] + " | cut -f1) > " + outfile
            else:
                cmd = ("LC_ALL=C sort -um <(zcat " + infiles[0] + " | cut -f1) <(zcat " + infiles[1] + " | cut -f1)"
                       + "".join(" | LC_ALL=C sort -um - <(zcat " + f + " | cut -f1)" for f in infiles[2:])
                       + " > " + outfile)
            rcompat.r_system2("/bin/bash", "-c '" + cmd + "'")
            if len(infiles) > 1:
                dedupe(outfile)
            rcompat.r_system2("/bin/bash", "-c 'touch " + outfile_completed + "'")
        else:
            # Subsequent rounds: merge merged files
            te = math.ceil(t / (b ** (i - 1))) * int(b ** (i - 1))
            outfile = r_paste0(output_dir, output_prefix, ".protein", kmer_length, ".j.", i, ".", te, ".txt")
            r_cat("Creating temp file:", outfile, "\n")
            outfile_completed = r_paste0(output_dir, output_prefix, ".protein", kmer_length, ".j.", i, ".", te,
                                         ".completed.txt")
            beg = te - float(b) ** (i - 1) + float(b) ** (i - 2)
            end = min(te, math.ceil(t / b ** (i - 2)) * int(b ** (i - 2)))
            r_cat(r_paste0("Beg", i, ":"), beg, r_paste0("End", i, ":"), end, "\n")
            if end < beg:
                r_stop("Problem with input arguments, please check")
            inc = float(b) ** (i - 2)
            js = rcompat.r_seq(beg, end, inc)
            infiles = [r_paste0(output_dir, output_prefix, ".protein", kmer_length, ".j.", i - 1, ".", j, ".txt")
                       for j in js]
            infiles_completed = [r_paste0(output_dir, output_prefix, ".protein", kmer_length, ".j.", i - 1, ".", j,
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
            if len(infiles) > 1:
                dedupe(outfile)
            # If length of infiles is one, it has been moved and doesn't exist. If more than one, delete the temp files.
            if len(infiles) > 1:
                cmd = " ".join(["rm", " ".join(infiles)])
                rcompat.r_system2("/bin/bash", "-c '" + cmd + "'")
            rcompat.r_system2("/bin/bash", "-c 'touch " + outfile_completed + "'")
    if t == p:
        ## Copy final file to end location
        create_final_file(outfile=outfile, output_dir=output_dir, output_prefix=output_prefix, kmer_length=kmer_length)
        ## Sort final kmer file
        sort_final_file(output_dir=output_dir, output_prefix=output_prefix, sort_strings=sort_strings,
                        kmer_length=kmer_length)

    r_cat("Finished in", (time.monotonic() - start_time) / 60, "minutes\n")


if __name__ == "__main__":
    main()
