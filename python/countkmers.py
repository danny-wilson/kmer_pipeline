#!/usr/bin/env python3
"""countkmers.py: count nucleotide or protein k-mers for one sample.
Port of countkmers.Rscript."""
import argparse
import collections
import os
import sys
import time

import rcompat
from rcompat import r_cat, r_stop

###################################################################################################
## Functions and software paths
###################################################################################################


def proteinKmerCount(protein, kmerLen):
    """All k-mers of a protein, or None if it is shorter than kmerLen."""
    if len(protein) >= kmerLen:
        return [protein[s:s + kmerLen] for s in range(len(protein) - (kmerLen - 1))]
    return None


def count_protein_kmers(fastaFile, writeToFile=False, kmerLen=31, kmerDir=None, id=None, oneLetterCodes=None):
    # Read in contigs
    proteins = rcompat.r_scan_lines(fastaFile, quiet=True)
    # Only keep the rows with the translated contigs in, not the fasta header rows
    proteins = [x for x in proteins if ">" not in x]
    allowed = set(oneLetterCodes.values())
    if not set("".join(proteins)) <= allowed:
        r_stop("Error: unexpected characters in proteins", "\n")
    # Count kmers for all contigs
    k = int(kmerLen)
    kmers = collections.Counter()
    for x in proteins:
        km = proteinKmerCount(x, k)
        if km is not None:
            kmers.update(km)
    # Write kmers to file
    if writeToFile:
        if not os.path.exists(kmerDir):
            r_stop("Error: kmer file directory doesn't exist", "\n")
        final_kmer_file = rcompat.r_paste0(kmerDir, id, ".kmer", kmerLen, ".unsorted.txt")
        # R writes table() order; the file is sorted by sort_strings next, so the order here does not matter
        with open(final_kmer_file, "w") as f:
            f.write("".join(f"{km}\t{n}\n" for km, n in kmers.items()))
        rcompat.r_system("gzip " + final_kmer_file)
        return final_kmer_file + ".gz"
    # Or return counted kmer sequences
    else:
        return kmers


def create_kmercount_dir(dir, kmertype, kmerlength):
    kmer_dir = dir + "/" + rcompat.r_paste0(kmertype, "kmer", kmerlength, "/")  # file.path
    if not os.path.isdir(kmer_dir):
        os.mkdir(kmer_dir)
    r_cat(rcompat.r_paste0("Counting ", kmertype, " kmers of length ", kmerlength), "\n")
    return kmer_dir


def nucleotide_kmer_counting(output_dir, kmerlength, dsk_path, dsk2ascii_path, contig_path, sample_id):
    # If it doesn't exist, create a directory to contain kmers for this kmer length
    kmer_dir = create_kmercount_dir(dir=output_dir, kmertype="nucleotide", kmerlength=kmerlength)
    os.chdir(kmer_dir)
    k = rcompat.r_as_character(kmerlength)
    # Run DSK
    dsktmpdir = sample_id + "-countkmers-tmpdir"
    os.mkdir(dsktmpdir)
    rcompat.r_system(dsk_path + " -file " + contig_path + " -kmer-size " + k + " -max-disk 0 -abundance-min 1 -out "
                     + sample_id + " -out-tmp " + dsktmpdir)
    import shutil
    shutil.rmtree(dsktmpdir, ignore_errors=True)
    # Convert DSK format into text file
    rcompat.r_system(dsk2ascii_path + " -file " + sample_id + ".h5 -out " + sample_id + ".kmer" + k + "_unsorted.txt")
    # Sort the kmers
    rcompat.r_system("sort -k 1 " + sample_id + ".kmer" + k + "_unsorted.txt > " + sample_id + ".kmer" + k + ".txt")
    # Write the number of kmers to file
    nKmers = float(rcompat.r_system_intern("wc -l " + sample_id + ".kmer" + k + ".txt")[0].split(" ")[0])
    rcompat.r_cat_lines([rcompat.r_paste_collapse(["Total", rcompat.r_as_character(nKmers)], "\t")],
                        sample_id + ".kmer" + k + ".total.txt")
    # Gzip the sorted kmer file
    rcompat.r_system("gzip " + sample_id + ".kmer" + k + ".txt")
    # Remove intermediate files
    rcompat.r_system("rm " + sample_id + ".h5")
    rcompat.r_system("rm " + sample_id + ".kmer" + k + "_unsorted.txt")
    kmerfile = kmer_dir + sample_id + ".kmer" + k
    r_cat("Counted nucleotide kmers length " + k + " for sample ID " + sample_id + ". Output files: "
          + kmerfile + ".txt.gz " + kmerfile + ".total.txt", "\n")
    return kmer_dir


def protein_kmer_counting(output_dir, kmerlength, translated_contigs_path, sample_id, oneLetterCodes, sort_strings):
    # If it doesn't exist, create a directory to contain kmers for this kmer length
    kmer_dir = create_kmercount_dir(dir=output_dir, kmertype="protein", kmerlength=kmerlength)
    k = rcompat.r_as_character(kmerlength)
    # Count protein kmers - unsorted
    kmerfile_unsorted = count_protein_kmers(fastaFile=translated_contigs_path, writeToFile=True, kmerLen=kmerlength,
                                            kmerDir=kmer_dir, id=sample_id, oneLetterCodes=oneLetterCodes)
    # Sort the kmer file
    sorted_kmerfile = kmer_dir + sample_id + ".kmer" + k + ".txt.gz"
    sortCommand = " ".join([sort_strings, kmerfile_unsorted, "| gzip -c >", sorted_kmerfile])
    rcompat.r_system(sortCommand)
    rcompat.r_system("rm " + kmerfile_unsorted)
    # Write the number of kmers to file
    nKmers = float(rcompat.r_system_intern("zcat " + sorted_kmerfile + " | wc -l")[0].split(" ")[0])
    rcompat.r_cat_lines([rcompat.r_paste_collapse(["Total", rcompat.r_as_character(nKmers)], "\t")],
                        kmer_dir + sample_id + ".kmer" + k + ".total.txt")
    r_cat("Counted protein kmers length " + k + " for sample ID " + sample_id + ". Output file: " + sorted_kmerfile, "\n")
    return kmer_dir


def write_kmer_filepaths_to_file(process, id_file, kmertype, kmerlength, output_dir, output_prefix, kmer_dir):
    # For the last file count, create file containing paths to all output files
    # Won't check if they are all created - flag warning
    if process == len(id_file):
        k = rcompat.r_as_character(kmerlength)
        all_outfiles = [kmer_dir + rcompat.r_as_character(i) + ".kmer" + k + ".txt.gz" for i in id_file["id"]]
        all_outfiles_path = output_dir + "/" + output_prefix + "_" + kmertype + k + "_kmers_filepaths.txt"
        r_cat("Writing file paths to all kmer counts for " + kmertype + " kmer length " + k,
              "to file (Warning: have not checked that all kmer counting is completed and all files exist):",
              all_outfiles_path, "\n")
        rcompat.r_cat_lines(all_outfiles, all_outfiles_path)


###################################################################################################


def main():
    rcompat.script_setup(__file__)
    parser = argparse.ArgumentParser(description="countkmers.py count nucleotide or protein kmers",
                                     allow_abbrev=False)
    parser.add_argument("--task-id", required=True, type=int, help="row of the ID file to process (from 1)")
    parser.add_argument("--id-file", required=True, help="sample ID file (columns id, paths, pheno)")
    parser.add_argument("--analysis-dir", required=True, help="analysis directory")
    parser.add_argument("--output-prefix", required=True, help="output prefix")
    parser.add_argument("--software-file", required=True, help="software paths file")
    parser.add_argument("--analyses-list", default=None,
                        help="file of k-mer types and lengths to count (default: nucleotide k-mers of length 31)")
    args = parser.parse_args()

    start_time = time.monotonic()

    process = args.task_id
    id_file = args.id_file
    output_dir = args.analysis_dir
    output_prefix = args.output_prefix
    software_file = args.software_file
    analyses_list = args.analyses_list

    # Check input arguments
    if not os.path.exists(output_dir):
        r_stop("Error: output directory doesn't exist", "\n")
    if output_dir.endswith("/"):
        output_dir = output_dir[:-1]
    if not os.path.exists(id_file):
        r_stop("Error: sample ID file doesn't exist", "\n")
    if not os.path.exists(software_file):
        r_stop("Error: software file doesn't exist", "\n")
    if analyses_list is not None:
        if not os.path.exists(analyses_list):
            r_stop("Error: analyses list file doesn't exist", "\n")
    else:
        r_cat("No kmer type input so counting nucleotide kmers of length 31bp", "\n")

    # Read in software file
    software_paths = rcompat.r_read_table(software_file, header=True, sep="\t", quote="")
    names = [rcompat.r_as_character(n) for n in software_paths["name"]]
    paths = [rcompat.r_as_character(p) for p in software_paths["path"]]

    def software(name):  # as.character(software_paths$path)[which(tolower(name) == name)]
        hits = [p for n, p in zip(names, paths) if n.lower() == name]
        if len(hits) != 1:
            r_stop("Error: expected one software path for ", name, " in the software file")
        return hits[0]

    # Required software and script paths
    # Begin with software
    required_software = ["scriptpath", "dsk", "dsk2ascii"]
    if any(r not in names for r in required_software):
        r_stop("Error: missing required software path in the software file - requires " + ", ".join(required_software), "\n")
    dsk_path = software("dsk")
    if not os.path.exists(dsk_path):
        r_stop("Error: dsk path doesn't exist", "\n")
    dsk2ascii_path = software("dsk2ascii")
    if not os.path.exists(dsk2ascii_path):
        r_stop("Error: dsk2ascii path doesn't exist", "\n")
    # Get the script location and then create missing paths
    script_location = software("scriptpath")
    if not os.path.isdir(script_location):
        r_stop("Error: script location directory specified in the software paths file doesn't exist", "\n")
    sort_strings = script_location + "/sort_strings"
    if not os.path.exists(sort_strings):
        r_stop("Error: sort_strings path doesn't exist - check pipeline script location in the software file", "\n")
    sequence_functions_file = script_location + "/sequence_functions.py"
    if not os.path.exists(sequence_functions_file):
        r_stop("Error: sequence_functions.py path doesn't exist - check pipeline script location in the software file", "\n")
    sys.path.insert(0, script_location)
    import sequence_functions

    # Report variables
    r_cat("#############################################", "\n")
    r_cat("Running on host: ", rcompat.r_system_intern("hostname"), "\n")
    r_cat("Command line arguments", "\n")
    r_cat(sys.argv[1:], "\n\n")
    r_cat("Parameters:", "\n")
    r_cat("task_id:", process, "\n")
    r_cat("ID file path:", id_file, "\n")
    r_cat("Analysis directory:", output_dir, "\n")
    r_cat("Output prefix:", output_prefix, "\n")
    r_cat("Software file:", software_file, "\n")
    r_cat("Script location:", script_location, "\n")
    r_cat("dsk path:", dsk_path, "\n")
    r_cat("dsk2ascii path:", dsk2ascii_path, "\n")
    r_cat("Analyses list file path:", [] if analyses_list is None else analyses_list, "\n")
    r_cat("#############################################", "\n\n")

    # Read in ID file
    id_file = rcompat.r_read_table(id_file, header=True, sep="\t")
    # Just keep the sample ID of the current process
    if not 1 <= process <= len(id_file):
        r_stop("Error: task_id ", process, " is not a row of the ID file")
    sample_id = rcompat.r_as_character(id_file["id"].iloc[process - 1])
    contig_path = rcompat.r_as_character(id_file["paths"].iloc[process - 1])
    # Check that the contig file exists
    if not os.path.exists(contig_path):
        r_stop("Error: contig path", contig_path, "doesn't exist", "\n")

    if not os.path.isdir(output_dir + "/translated_contigs"):
        os.mkdir(output_dir + "/translated_contigs")

    # Read in analyses to run
    if analyses_list is not None:
        kt = rcompat.r_read_table(analyses_list, header=True, sep="\t")
        types = [rcompat.r_as_character(t).lower() for t in kt.iloc[:, 0]]
        if any(t not in ("protein", "nucleotide") for t in types):
            r_stop("Error: kmer type must be either 'nucleotide' or 'protein'", "\n")
        lengths = []
        for v in kt.iloc[:, 1]:
            try:  # as.integer(): NA for anything that isn't a number
                lengths.append(float(rcompat.r_as_character(v)))
            except ValueError:
                r_stop("Error: kmer length must be an integer", "\n")
        kmertype = list(zip(types, lengths))
    else:
        kmertype = [("nucleotide", 31.0)]

    if any(t == "protein" for t, _ in kmertype):
        ## For each assembly, translate all contigs into 6 reading frames
        r_cat("Translating contigs for ID", sample_id, "\n")
        # Translate contigs
        translated_contigs_path = sequence_functions.translate_6_frames(
            contig_path=contig_path, id=sample_id, outDir=output_dir + "/translated_contigs/",
            oneLetterCodes=sequence_functions.oneLetterCodes, revcompl=sequence_functions.revcompl)
        r_cat("Translated contigs for ID " + sample_id + ". Output file: " + translated_contigs_path, "\n")

        # Protein kmer lengths to count
        for kl in [l for t, l in kmertype if t == "protein"]:
            kmer_dir = protein_kmer_counting(output_dir=output_dir, kmerlength=kl,
                                             translated_contigs_path=translated_contigs_path, sample_id=sample_id,
                                             oneLetterCodes=sequence_functions.oneLetterCodes,
                                             sort_strings=sort_strings)
            write_kmer_filepaths_to_file(process=process, id_file=id_file, kmertype="protein", kmerlength=kl,
                                         output_dir=output_dir, output_prefix=output_prefix, kmer_dir=kmer_dir)
        r_cat("\n")
        r_cat("#############################################", "\n\n")

    if any(t == "nucleotide" for t, _ in kmertype):
        for kl in [l for t, l in kmertype if t == "nucleotide"]:
            kmer_dir = nucleotide_kmer_counting(output_dir=output_dir, kmerlength=kl, dsk_path=dsk_path,
                                                dsk2ascii_path=dsk2ascii_path, contig_path=contig_path,
                                                sample_id=sample_id)
            write_kmer_filepaths_to_file(process=process, id_file=id_file, kmertype="nucleotide", kmerlength=kl,
                                         output_dir=output_dir, output_prefix=output_prefix, kmer_dir=kmer_dir)
        r_cat("#############################################", "\n\n")

    r_cat("Completed in", (time.monotonic() - start_time) / 60, "minutes\n")


if __name__ == "__main__":
    main()
