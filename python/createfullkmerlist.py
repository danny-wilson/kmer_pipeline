#!/usr/bin/env python3
"""createfullkmerlist.py: merge k-mers into a list of unique k-mers present across
the dataset, by running proteinkmermerge.py or nucleotidekmermerge.py for this task.
Port of createfullkmerlist.Rscript."""
import argparse
import os
import sys
import time

import rcompat
from rcompat import r_cat, r_stop


def main():
    rcompat.script_setup(__file__)
    start_time = time.monotonic()  # noqa: F841 (R records it but never reports it)
    parser = argparse.ArgumentParser(
        description="createfullkmerlist.py merge kmers into a list of unique kmers present across the dataset",
        allow_abbrev=False)
    parser.add_argument("--task-id", required=True)
    parser.add_argument("--n", required=True, help="number of samples")
    parser.add_argument("--p", required=True, help="number of parallel tasks")
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--analysis-dir", required=True)
    parser.add_argument("--id-file", required=True)
    parser.add_argument("--kmer-type", required=True, help="protein or nucleotide")
    parser.add_argument("--kmer-length", required=True)
    parser.add_argument("--software-file", required=True)
    args = parser.parse_args()

    # Initialize variables
    process = rcompat.r_as_integer(args.task_id)
    n = rcompat.r_as_integer(args.n)
    p = rcompat.r_as_integer(args.p)
    output_prefix = args.output_prefix
    output_dir = args.analysis_dir
    id_file = args.id_file
    kmer_type = args.kmer_type.lower()
    kmer_length = rcompat.r_as_integer(args.kmer_length)
    software_file = args.software_file

    # Check input arguments
    if not os.path.exists(output_dir):
        r_stop("Error: output directory doesn't exist", "\n")
    if not output_dir.endswith("/"):
        output_dir = output_dir + "/"
    if not os.path.exists(id_file):
        r_stop("Error: sample ID file doesn't exist", "\n")
    if kmer_type != "protein" and kmer_type != "nucleotide":
        r_stop("Error: kmer type must be either 'protein' or 'nucleotide'", "\n")
    if kmer_length is None:
        r_stop("Error: kmer length must be an integer", "\n")
    # Get directory containing the kmers: file.path(output_dir, "<type>kmer<k>", "/"), keeping R's "//"
    input_dir = output_dir + "/" + rcompat.r_paste0(kmer_type, "kmer", kmer_length) + "/" + "/"
    if not os.path.exists(input_dir):
        r_stop("Error: input directory " + input_dir + " doesn't exist", "\n")

    # Read in software file
    software_paths = rcompat.r_read_table(software_file, header=True, sep="\t", quote="")
    names = [rcompat.r_as_character(v) for v in software_paths["name"]]
    paths = [rcompat.r_as_character(v) for v in software_paths["path"]]
    # The R original launched its children with the "R" entry's Rscript; the Python
    # children run with this interpreter, but software files keep the same entries.
    required_software = ["R", "scriptpath"]
    if any(r not in names for r in required_software):
        r_stop("Error: missing required software path in the software file - requires " + ", ".join(required_software), "\n")
    python_path = sys.executable
    script_location = [pth for nm, pth in zip(names, paths) if nm.lower() == "scriptpath"][0]
    if not os.path.isdir(script_location):
        r_stop("Error: script location directory specified in the software paths file doesn't exist", "\n")
    proteinkmermergescript = script_location + "/proteinkmermerge.py"
    if not os.path.exists(proteinkmermergescript):
        r_stop("Error: proteinkmermerge.py path doesn't exist - check pipeline script location in the software file", "\n")
    nucleotidekmermergescript = script_location + "/nucleotidekmermerge.py"
    if not os.path.exists(nucleotidekmermergescript):
        r_stop("Error: nucleotidekmermerge.py path doesn't exist - check pipeline script location in the software file", "\n")

    # Report variables
    r_cat("#############################################", "\n")
    r_cat("Running on host: ", rcompat.r_system_intern("hostname"), "\n")
    r_cat("Command line arguments", "\n")
    r_cat(sys.argv[1:], "\n\n")
    r_cat("Parameters:", "\n")
    r_cat("task_id:", process, "\n")
    r_cat("n:", n, "\n")
    r_cat("p:", p, "\n")
    r_cat("Output prefix:", output_prefix, "\n")
    r_cat("Analysis directory:", output_dir, "\n")
    r_cat("ID file path:", id_file, "\n")
    r_cat("Kmer type:", kmer_type, "\n")
    r_cat("Kmer length", kmer_length, "\n")
    r_cat("Software file:", software_file, "\n")
    r_cat("Script location:", script_location, "\n")
    r_cat("Python path:", python_path, "\n")
    r_cat("#############################################", "\n\n")

    def arg(v):
        return "NA" if v is None else rcompat.r_as_character(v)

    if kmer_type == "protein":
        rcompat.r_system(python_path + " " + proteinkmermergescript + " --n " + arg(n) + " --p " + arg(p)
                         + " --output-prefix " + output_prefix + " --input-dir " + input_dir + " --output-dir " + output_dir
                         + " --id-file " + id_file + " --kmer-length " + arg(kmer_length)
                         + " --software-file " + software_file + " --process " + arg(process))

    if kmer_type == "nucleotide":
        rcompat.r_system(python_path + " " + nucleotidekmermergescript + " --n " + arg(n) + " --p " + arg(p)
                         + " --output-prefix " + output_prefix + " --input-dir " + input_dir + " --output-dir " + output_dir
                         + " --id-file " + id_file + " --kmer-length " + arg(kmer_length)
                         + " --process " + arg(process))


if __name__ == "__main__":
    main()
