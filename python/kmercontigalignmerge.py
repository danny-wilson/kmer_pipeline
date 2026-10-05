#!/usr/bin/env python3
"""kmercontigalignmerge.py: merge the k-mer/gene combinations of all samples into a
sorted list, then count, in batches shared between tasks, the samples each
combination is present in. Port of kmercontigalignmerge.Rscript (adapted from
Daniel Wilson's kmermerge.Rscript, 2018)."""
import argparse
import math
import os
import sys
import time

import rcompat
from rcompat import r_cat, r_paste0, r_stop

###################################################################################################
## Functions and software paths
###################################################################################################


def get_count_batch_parameters(p, b, n):
    """Rows (beg, end) for each task, as integers, written as plain integers (D1a)."""
    return [((t - 1) * b + 1, min(t * b, n)) for t in range(1, int(p) + 1)]


def format_count(v):
    """A count: plain integer when it has an integer value (D1a; R's cat() wrote
    100000 as "1e+05"), otherwise as R's cat()."""
    if v == v and v == int(v):
        return str(int(v))
    return rcompat.r_str(v, 7)


###################################################################################################


def main():
    rcompat.script_setup(__file__)
    start_time = time.monotonic()
    parser = argparse.ArgumentParser(description="kmercontigalignmerge.py merge kmer/gene alignment combinations. "
                                                 "Adapted from Daniel Wilson (2018) kmermerge.Rscript",
                                     allow_abbrev=False)
    parser.add_argument("--task-id", required=True)
    parser.add_argument("--n", required=True, help="number of samples")
    parser.add_argument("--p", required=True, help="number of parallel tasks")
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--analysis-dir", required=True)
    parser.add_argument("--input-files", required=True, help="file listing the per-sample k-mer/gene files")
    parser.add_argument("--kmer-type", required=True, help="protein or nucleotide")
    parser.add_argument("--kmer-length", required=True)
    parser.add_argument("--ref-fa", required=True)
    parser.add_argument("--nucmerident", required=True)
    parser.add_argument("--software-file", required=True)
    parser.add_argument("--merge-wait-minutes", default="100",
                        help="minutes to wait for files written by other tasks before stopping (default 100)")
    args = parser.parse_args()
    rcompat.set_merge_wait(args.merge_wait_minutes)

    # Initialize variables
    process = rcompat.r_as_integer(args.task_id)
    n = rcompat.r_as_integer(args.n)
    p = rcompat.r_as_integer(args.p)
    out_prefix = args.output_prefix
    output_dir = args.analysis_dir
    input_files = args.input_files
    kmer_type = args.kmer_type.lower()
    kmer_length = rcompat.r_as_integer(args.kmer_length)
    ref_fa = args.ref_fa
    ident_threshold = rcompat.r_as_numeric(args.nucmerident)
    software_file = args.software_file

    if n is None:
        r_stop("Error: n must be an integer", "\n")
    if p is None:
        r_stop("Error: p must be an integer", "\n")
    if not os.path.exists(output_dir):
        r_stop("Error: output directory doesn't exist", "\n")
    if not output_dir.endswith("/"):
        output_dir = output_dir + "/"
    if not os.path.exists(input_files):
        r_stop("Error: input file paths file doesn't exist", "\n")
    if kmer_type != "protein" and kmer_type != "nucleotide":
        r_stop("Error: variant type must be either protein or nucleotide", "\n")
    if kmer_length is None:
        r_stop("Error: kmer length must be an integer", "\n")
    if not os.path.exists(ref_fa):
        r_stop("Error: reference fasta file doesn't exist", "\n")
    if ident_threshold is None or ident_threshold > 100 or ident_threshold < 0:
        r_stop("Error: nucmer identity threshold must be between 0-100", "\n")
    if not os.path.exists(software_file):
        r_stop("Error: software file doesn't exist", "\n")

    # Input files
    input_files_path = input_files
    # Contig align directory
    contigalign_dir = output_dir + "/" + r_paste0(kmer_type, "kmer", kmer_length, "_kmergenealign/")  # file.path
    if not os.path.isdir(contigalign_dir):
        r_stop("Error: contig align directory", contigalign_dir, "doesn't exist", "\n")

    # Read in software file
    software_paths = rcompat.r_read_table(software_file, header=True, sep="\t", quote="")
    names = [rcompat.r_as_character(v) for v in software_paths["name"]]
    paths = [rcompat.r_as_character(v) for v in software_paths["path"]]
    required_software = ["scriptpath"]
    if any(r not in names for r in required_software):
        r_stop("Error: missing required software path in the software file - requires " + ", ".join(required_software), "\n")
    script_location = [pth for nm, pth in zip(names, paths) if nm.lower() == "scriptpath"][0]
    if not os.path.isdir(script_location):
        r_stop("Error: script location directory specified in the software paths file doesn't exist", "\n")
    sortstringspath = script_location + "/sort_strings"
    if not os.path.exists(sortstringspath):
        r_stop("Error: sort_strings path doesn't exist - check pipeline script location in the software file", "\n")
    stringlist2countpath = script_location + "/stringlist2count"
    if not os.path.exists(stringlist2countpath):
        r_stop("Error: stringlist2count path doesn't exist - check pipeline script location in the software file", "\n")

    # Report variables
    r_cat("#############################################", "\n")
    r_cat("Running on host: ", rcompat.r_system_intern("hostname"), "\n")
    r_cat("Command line arguments", "\n")
    r_cat(sys.argv[1:], "\n\n")
    r_cat("Parameters:", "\n")
    r_cat("task_id:", process, "\n")
    r_cat("n:", n, "\n")
    r_cat("p:", p, "\n")
    r_cat("Output prefix:", out_prefix, "\n")
    r_cat("Analysis directory:", output_dir, "\n")
    r_cat("Input files path:", input_files, "\n")
    r_cat("Kmer type:", kmer_type, "\n")
    r_cat("Kmer length:", kmer_length, "\n")
    r_cat("Reference fasta file:", ref_fa, "\n")
    r_cat("Nucmer alignment minimum % identity:", ident_threshold, "\n")
    r_cat("Software file:", software_file, "\n")
    r_cat("Script location:", script_location, "\n")
    r_cat("#############################################", "\n\n")

    # Read in sample IDs
    input_files = rcompat.r_scan_lines(input_files, quiet=True)

    if n != len(input_files):
        r_stop("n does not equal the number of samples in the input path file", "\n")
    if not all(os.path.exists(f) for f in input_files):
        r_stop("Error: not all input files exist", "\n")

    # Get the reference name
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import kmercontigalignonly
    ref_name = kmercontigalignonly.read_reference_name(ref_fa)
    r_cat("Reference name:", ref_name, "\n")

    b = math.ceil(n / p)
    if b == 1:
        r_stop("Cannot have batchsize = 1. Try p < n/2")
    p = math.ceil(n / b)
    t = process
    imax = math.ceil(math.log(n) / math.log(b))

    r_cat("Parameters [n batchsize processes taskid imax]: ", n, b, p, t, imax, "\n")

    if t > p:
        r_cat("Task", t, "not required\n")
        return

    stem = r_paste0(contigalign_dir, out_prefix, "_", kmer_type, kmer_length)
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
            outfile = r_paste0(stem, ".j.", i, ".", t, ".txt")
            outfile_completed = r_paste0(stem, ".j.", i, ".", t, ".completed.txt")
            beg = b * (t - 1) + 1
            end = min(b * t, n)
            r_cat("Beg:", beg, " End:", end, "\n")
            if end < beg:
                r_stop("Problem with input arguments, please check")

            infiles = ["NA" if f is None else f for f in rcompat.r_index(input_files, rcompat.r_colon(int(beg), end))]

            infiles_size = [size_of(rcompat.r_system_intern("ls -l " + x + " | cut -d ' ' -f5")) for x in infiles]

            if not all(os.path.exists(f) for f in infiles) or not r_all_positive(infiles_size):
                r_stop("Could not find files or files empty ", " ".join(infiles))

            if len(infiles) == 1:
                cmd = "zcat " + infiles[0] + " | cut -f1 > " + outfile
            elif len(infiles) == 2:
                cmd = "sort -u <(zcat " + infiles[0] + " | cut -f1) <(zcat " + infiles[1] + " | cut -f1) > " + outfile
            else:
                cmd = ("sort -u <(zcat " + infiles[0] + " | cut -f1) <(zcat " + infiles[1] + " | cut -f1)"
                       + "".join(" | sort -u - <(zcat " + f + " | cut -f1)" for f in infiles[2:]) + " > " + outfile)
            rcompat.r_system2("/bin/bash", "-c '" + cmd + "'")

            rcompat.r_system2("/bin/bash", "-c 'touch " + outfile_completed + "'")
        else:
            # Subsequent rounds: merge merged files
            te = math.ceil(t / (b ** (i - 1))) * int(b ** (i - 1))
            outfile = r_paste0(stem, ".j.", i, ".", te, ".txt")
            outfile_completed = r_paste0(stem, ".j.", i, ".", te, ".completed.txt")
            beg = te - float(b) ** (i - 1) + float(b) ** (i - 2)
            end = min(te, math.ceil(t / b ** (i - 2)) * int(b ** (i - 2)))
            r_cat(r_paste0("Beg", i, ":"), beg, r_paste0("End", i, ":"), end, "\n")
            if end < beg:
                r_stop("Problem with input arguments, please check")
            inc = float(b) ** (i - 2)
            js = rcompat.r_seq(beg, end, inc)
            infiles = [r_paste0(stem, ".j.", i - 1, ".", x, ".txt") for x in js]
            infiles_completed = [r_paste0(stem, ".j.", i - 1, ".", x, ".completed.txt") for x in js]
            nattempts = 0
            while not all(os.path.exists(f) for f in infiles_completed) or not all(os.path.exists(f) for f in infiles):
                nattempts = nattempts + 1
                if rcompat.wait_exceeded(nattempts, 60):
                    r_stop("Could not find files ", " ".join(infiles_completed), rcompat.wait_message())
                time.sleep(60)

            infiles_size = [size_of(rcompat.r_system_intern("ls -l " + x + " | cut -d ' ' -f5")) for x in infiles]
            if r_any_zero(infiles_size):
                r_stop("One or more file size is zero ", "".join(infiles), "\n")

            if len(infiles) == 1:
                cmd = "mv " + infiles[0] + " " + outfile
            elif len(infiles) == 2:
                cmd = "sort -u " + infiles[0] + " " + infiles[1] + " > " + outfile
            else:
                cmd = ("sort -u " + infiles[0] + " " + infiles[1]
                       + "".join(" | sort -u - " + f for f in infiles[2:]) + " > " + outfile)
            rcompat.r_system2("/bin/bash", "-c '" + cmd + "'")
            # If length of infiles is one, it has been moved and doesn't exist. If more than one, delete the temp files.
            if len(infiles) > 1:
                cmd = " ".join(["rm", " ".join(infiles)])
                rcompat.r_system2("/bin/bash", "-c '" + cmd + "'")
            rcompat.r_system2("/bin/bash", "-c 'touch " + outfile_completed + "'")
    if t == p:
        # Check outfile isn't empty
        outfile_size = rcompat.r_system_intern("ls -l " + outfile + " | cut -d ' ' -f5")
        if outfile_size == ["0"]:
            r_stop(outfile, " file is empty", "\n")

        final_kmer_txt_gz = r_paste0(output_dir, out_prefix, "_", kmer_type, kmer_length, ".", ref_name, "_t",
                                     ident_threshold, ".kmeralignmerge.txt.gz")
        nKmersFinal = rcompat.r_system_intern("cat " + outfile + " | wc -l")
        dummycount = ["1"] * int(nKmersFinal[0])
        dummycountfile = stem + "_dummycount.txt"
        r_cat("Written dummy count", "\n")
        rcompat.r_cat_lines(dummycount, dummycountfile)
        outfile_count = stem + "_kmeralignmerge_dummycount.txt"
        r_cat("outfile:", outfile, "\n")
        r_cat("dummycountfile:", dummycountfile, "\n")
        r_cat("outfile_count:", outfile_count, "\n")
        rcompat.r_system("paste " + outfile + " " + dummycountfile + " > " + outfile_count)
        sortCommand = " ".join([sortstringspath, outfile_count, "| cut -f1 | gzip -c >", final_kmer_txt_gz])
        r_cat("Sort command:", "\n")
        rcompat.r_system(sortCommand)
        r_cat(sortCommand, "\n")

        # Remove all completed files
        completed_files = rcompat.r_system_intern("ls " + stem + "*.j*.completed.txt")
        completed_files = completed_files + rcompat.r_system_intern("ls " + stem + "*.kmercontigalign.completed.txt")
        cmd = " ".join(["rm", " ".join(completed_files)])
        rcompat.r_system2("/bin/bash", "-c '" + cmd + "'")
        rcompat.r_system("rm " + outfile)
        rcompat.r_system("rm " + stem + "_dummycount.txt")
        rcompat.r_system("rm " + stem + "_kmeralignmerge_dummycount.txt")
        # Write file for final file being completed
        outfile_completed = r_paste0(stem, ".", ref_name, "_t", ident_threshold, ".final.kmeralignmerge.completed.txt")
        rcompat.r_system2("/bin/bash", "-c 'touch " + outfile_completed + "'")

    r_cat("Process finished for kmeralignmerge", "\n")
    # Get counts for kmer/gene combinations
    r_cat("Getting counts for kmeralignmerge", "\n")

    infile_completed = r_paste0(stem, ".", ref_name, "_t", ident_threshold, ".final.kmeralignmerge.completed.txt")
    infile = r_paste0(output_dir, out_prefix, "_", kmer_type, kmer_length, ".", ref_name, "_t", ident_threshold,
                      ".kmeralignmerge.txt.gz")

    nattempts = 0
    while not os.path.exists(infile_completed) or not os.path.exists(infile):
        nattempts = nattempts + 1
        if rcompat.wait_exceeded(nattempts, 60):
            r_stop("Could not find file ", infile_completed, rcompat.wait_message())
        time.sleep(60)

    # Read kmers
    # Total number of kmers
    n = rcompat.r_as_integer(rcompat.r_pipe("zcat " + infile + " | wc -l").split()[0])
    if n < 1:
        r_stop("No kmer/gene combinations found in", infile)
    # Number of kmers (batch size) per process
    b = math.ceil(n / p)

    # Define kmers to process
    params = get_count_batch_parameters(p=p, n=n, b=b)
    beg = params[t - 1][0]
    end = params[t - 1][1]
    out_prefix_counts = r_paste0(stem, ".", beg, "-", end)

    kmersublistfile = r_paste0(out_prefix_counts, ".", t, ".temp_kmerlist.txt.gz")
    cmd = r_paste0("zcat ", infile, " | head -n ", end, " | tail -n ", end - beg + 1, " | gzip -c > ", kmersublistfile)
    rcompat.r_system(cmd)

    cmd = rcompat.r_paste(stringlist2countpath, kmersublistfile, input_files_path, out_prefix_counts, 0.0, 1.0)
    rcompat.r_system(cmd)

    cmd = rcompat.r_paste("rm", kmersublistfile)
    rcompat.r_system(cmd)

    # Check that files have been created and are not empty
    outfiles = [out_prefix_counts + ".count.txt.gz"]
    outfiles_size = [size_of(rcompat.r_system_intern("ls -l " + x + " | cut -d ' ' -f5")) for x in outfiles]
    if r_any_zero(outfiles_size):
        r_stop("One or more JOB_INDEX ", t, " ", out_prefix_counts, " kmer align count files are empty")
    # Write file for final file being completed
    outfile_completed = out_prefix_counts + ".count.completed.txt"
    rcompat.r_system2("/bin/bash", "-c 'touch " + outfile_completed + "'")

    # Merge files and remove intermediate files
    if t == p:
        count_batch_parameters = get_count_batch_parameters(p=p, n=n, b=b)
        alignCountFiles = [r_paste0(stem, ".", b0, "-", b1, ".count.txt.gz") for b0, b1 in count_batch_parameters]
        alignCountCompletedFiles = [r_paste0(stem, ".", b0, "-", b1, ".count.completed.txt")
                                    for b0, b1 in count_batch_parameters]
        nattempts = 0
        while (not all(os.path.exists(f) for f in alignCountCompletedFiles)
               or not all(os.path.exists(f) for f in alignCountFiles)):
            nattempts = nattempts + 1
            if rcompat.wait_exceeded(nattempts, 60):
                r_stop("Could not find file ", "".join(alignCountCompletedFiles), rcompat.wait_message())
            time.sleep(60)

        count = []
        for f in alignCountFiles:  # scan(gzfile(f), sep = "\n"): numbers, one per line
            for tok in rcompat.r_scan_lines(f, quiet=True):
                v = rcompat.r_as_numeric(tok)
                if v is None and tok.strip() != "NA":
                    r_stop("scan() expected 'a real', got '", tok, "'")
                count.append(math.nan if v is None else v)
        final_count_file = r_paste0(output_dir, out_prefix, "_", kmer_type, kmer_length, ".", ref_name, "_t",
                                    ident_threshold, ".kmeralignmerge.count.txt")
        with open(final_count_file, "w") as fh:
            fh.write("".join(("NA" if v != v else format_count(v)) + "\n" for v in count) if count else "\n")
        rcompat.r_system("gzip " + final_count_file)
        # Delete intermediate count files
        r_cat("Deleting intermediate files", "\n")
        cmd = " ".join(["rm", " ".join(alignCountFiles + alignCountCompletedFiles + [infile_completed])])
        rcompat.r_system2("/bin/bash", "-c '" + cmd + "'")

    r_cat("Finished in", (time.monotonic() - start_time) / 60, "minutes\n")


def size_of(out):
    """as.numeric() of the output of system("ls -l <path> | cut -d ' ' -f5", intern = T):
    the size, or None (R's numeric(0)) when ls printed nothing."""
    if not out:
        return None
    v = rcompat.r_as_numeric(out[0])
    return math.nan if v is None else v


def r_all_positive(sizes):
    """all(sizes > 0) on R's sapply result; numeric(0) elements make it NA."""
    if any(s is not None and s == s and not s > 0 for s in sizes):
        return False
    if any(s is None or s != s for s in sizes):
        raise rcompat.RError("missing value where TRUE/FALSE needed")
    return True


def r_any_zero(sizes):
    if any(s == 0 for s in sizes if s is not None):
        return True
    if any(s is None or s != s for s in sizes):
        raise rcompat.RError("missing value where TRUE/FALSE needed")
    return False


if __name__ == "__main__":
    main()
