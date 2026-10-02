#!/usr/bin/env python3
"""runbowtie.py: map the k-mers to the reference genome with bowtie2 and keep the
mappings with MAPQ of at least samtools_filter. Port of runbowtie.Rscript. Not
called by kmer_pipeline.nf."""
import argparse
import os
import sys
import time

import rcompat
from rcompat import r_cat, r_stop


def main():
    rcompat.script_setup(__file__)
    start_time = time.monotonic()
    parser = argparse.ArgumentParser(description="runbowtie.py map kmers to the reference genome with bowtie2",
                                     allow_abbrev=False)
    for name in ("output-prefix", "analysis-dir", "kmerfile-prefix", "ref-fa", "kmer-type", "kmer-length",
                 "software-file"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--bowtie-parameters", default="--very-sensitive",
                        help="bowtie2 options, or a file containing them (default --very-sensitive); give options as --bowtie-parameters=VALUE")
    parser.add_argument("--samtools-filter", default="10", help="minimum MAPQ (default 10)")
    args = parser.parse_args()

    output_prefix = args.output_prefix
    output_dir = args.analysis_dir
    kmerfilePrefix = args.kmerfile_prefix
    ref_fa = args.ref_fa
    kmer_type = args.kmer_type
    kmer_length = rcompat.r_as_integer(args.kmer_length)
    software_file = args.software_file
    bowtie_parameters = args.bowtie_parameters
    samtools_filter = rcompat.r_as_integer(args.samtools_filter)

    if not os.path.exists(output_dir):
        r_stop("Error: output directory doesn't exist", "\n")
    if not output_dir.endswith("/"):
        output_dir = output_dir + "/"
    kmerSeqFile = kmerfilePrefix + ".kmermerge.txt.gz"
    if not os.path.exists(kmerSeqFile):
        r_stop("Error: kmer sequence file doesn't exist", "\n")
    if not os.path.exists(ref_fa):
        r_stop("Error: reference fasta file doesn't exist", "\n")
    if kmer_type != "protein" and kmer_type != "nucleotide":
        r_stop("Error: kmer type must be either 'protein' or 'nucleotide'", "\n")
    if kmer_length is None:
        r_stop("Error: kmer length must be an integer", "\n")
    if not os.path.exists(software_file):
        r_stop("Error: software file doesn't exist", "\n")
    if bowtie_parameters != "--very-sensitive":
        if not os.path.exists(bowtie_parameters):
            r_stop("Error: bowtie_parameters file doesn't exist", "\n")
        lines = rcompat.r_scan_lines(bowtie_parameters, quiet=False)
        if len(lines) != 1:
            r_stop("Error in system(...): 'command' must be a character string (", bowtie_parameters,
                   " has ", len(lines), " lines)")
        bowtie_parameters = lines[0]

    software_paths = rcompat.r_read_table(software_file, header=True, sep="\t", quote="")
    names = [rcompat.r_as_character(v) for v in software_paths["name"]]
    paths = [rcompat.r_as_character(v) for v in software_paths["path"]]
    required_software = ["bowtie2", "samtools"]
    if any(r not in names for r in required_software):
        r_stop("Error: missing required software path in the software file - requires " + ", ".join(required_software), "\n")

    def software(name):
        return [pth for nm, pth in zip(names, paths) if nm.lower() == name][0]
    bowtie2Dir = software("bowtie2")
    if not os.path.isdir(bowtie2Dir):
        r_stop("Error: bowtie2 installation directory specified in the software paths file doesn't exist", "\n")
    bowtie2path = bowtie2Dir + "/bowtie2"
    if not os.path.exists(bowtie2path):
        r_stop("Error: bowtie path", bowtie2path, " doesn't exist", "\n")
    bowtie2buildpath = bowtie2Dir + "/bowtie2-build"
    if not os.path.exists(bowtie2buildpath):
        r_stop("Error: bowtie2-build path", bowtie2buildpath, " doesn't exist", "\n")
    samtoolspath = software("samtools")
    if not os.path.exists(samtoolspath):
        r_stop("Error: samtoolspathpath specified in the software paths file doesn't exist", "\n")

    r_cat("#############################################", "\n")
    r_cat("Running on host: ", rcompat.r_system_intern("hostname"), "\n")
    r_cat("Command line arguments", "\n")
    r_cat(sys.argv[1:], "\n\n")
    r_cat("Parameters:", "\n")
    for label, v in (("Output prefix:", output_prefix), ("Analysis directory:", output_dir),
                     ("Kmer file prefix:", kmerfilePrefix), ("Reference fasta file:", ref_fa),
                     ("Kmer type:", kmer_type), ("Kmer length:", kmer_length), ("Software file:", software_file),
                     ("bowtie2 directory:", bowtie2Dir), ("Bowtie2 parameters:", bowtie_parameters),
                     ("Samtools path:", samtoolspath)):
        r_cat(label, v, "\n")
    r_cat("#############################################", "\n\n")

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import kmercontigalignonly
    ref_name = kmercontigalignonly.read_reference_name(ref_fa)

    bowtie_dir = output_dir + "/" + kmer_type + "kmer" + rcompat.r_as_character(kmer_length) + "_bowtie2mapping"
    if not os.path.isdir(bowtie_dir):
        rcompat.r_dir_create(bowtie_dir)

    os.chdir(bowtie_dir)
    rcompat.r_system(bowtie2buildpath + " -f " + ref_fa + " " + ref_name + "_bowtie_ref")

    r_cat("\n")
    r_cat("Running bowtie with the parameters:", bowtie_parameters, "\n")
    bowtie_outfile = output_prefix + "_" + kmer_type + rcompat.r_as_character(kmer_length) + "_map_to_" + ref_name
    rcompat.r_system(bowtie2path + " " + bowtie_parameters + " -r -x " + ref_name + "_bowtie_ref -U " + kmerSeqFile
                     + " -S " + bowtie_outfile)

    samtools_outputfile = rcompat.r_paste0(output_dir, output_prefix, "_", kmer_type, kmer_length, ".", ref_name, ".SAMq",
                                           samtools_filter, ".bowtie2map.txt")
    rcompat.r_system(samtoolspath + " view -q " + rcompat.r_as_character(samtools_filter) + " -S " + bowtie_outfile
                     + " > " + samtools_outputfile)
    rcompat.r_system("gzip " + samtools_outputfile)
    rcompat.r_system("gzip " + bowtie_outfile)

    r_cat("Completed in", (time.monotonic() - start_time) / 60, "minutes\n")


if __name__ == "__main__":
    main()
