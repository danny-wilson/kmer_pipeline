#!/usr/bin/env python3
"""kmercontigalign.py: align contigs to the reference genome using nucmer, assign
k-mers to genes or intergenic regions, then (for the first floor(n/5) tasks) run
kmercontigalignmerge.py. Port of kmercontigalign.Rscript, which is
kmercontigalignonly.Rscript followed by the merge; the shared body is in
kmercontigalignonly.py. Not called by kmer_pipeline.nf."""
import math
import os
import sys
import time

import rcompat
from rcompat import r_cat, r_paste0, r_stop

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmercontigalignonly  # noqa: E402


def merge(state):
    process, n = state["process"], state["n"]
    contigalign_dir, prefix = state["contigalign_dir"], state["prefix"]
    kmer_type, kmer_length = state["kmer_type"], state["kmer_length"]
    ref_name, ident_threshold, ids = state["ref_name"], state["ident_threshold"], state["ids"]

    # For a subset of the processes, use to merge the batches created by the other processes
    # Which processes to keep (floor() gives an R double)
    p = float(math.floor(n / 5))
    p = 1.0 if p == 0 else p

    if process > p:
        r_cat("Task", process, "not required for kmer contig align merging", "\n")
        r_cat("Finished in", (time.monotonic() - state["start_time"]) / 3600, "hours\n")
        return
    r_cat("Running kmer contig align merge", "\n\n")

    input_files_filepath = r_paste0(contigalign_dir, prefix, "_", kmer_type, kmer_length, "_", ref_name,
                                    "_kmergenecombination_filepaths.txt")
    final_file_prefix = r_paste0(contigalign_dir, prefix, "_", kmer_type, kmer_length, "_", ref_name, "_t",
                                 ident_threshold, "_")
    input_files = [final_file_prefix + i + "_nucmeralign_kmer_list_gene_IDs.txt.gz" for i in ids]
    input_files_completed = [r_paste0(contigalign_dir, prefix, "_", kmer_type, kmer_length, "_", i,
                                      ".kmercontigalign.completed.txt") for i in ids]

    nattempts = 0
    while not all(os.path.exists(f) for f in input_files_completed) or not all(os.path.exists(f) for f in input_files):
        nattempts = nattempts + 1
        if rcompat.wait_exceeded(nattempts, 60):
            r_stop("Could not find files", "".join(input_files_completed), rcompat.wait_message())
        time.sleep(60)

    rcompat.r_system(state["python_path"] + " " + state["kmercontigalignmergepath"] + " --task-id " + str(process)
                     + " --n " + str(n) + " --p " + rcompat.r_as_character(p) + " --output-prefix " + prefix
                     + " --analysis-dir " + state["output_dir"] + " --input-files " + input_files_filepath
                     + " --kmer-type " + kmer_type + " --kmer-length " + rcompat.r_as_character(kmer_length)
                     + " --ref-fa " + state["ref_fa"] + " --nucmerident " + rcompat.r_as_character(ident_threshold)
                     + " --software-file " + state["software_file"])
    r_cat("Finished in", (time.monotonic() - state["start_time"]) / 60, "minutes\n")


def main():
    kmercontigalignonly.run(__file__, "kmercontigalign.py align contigs to the reference genome using nucmer and "
                                      "assign kmers to genes or intergenic regions", merge_hook=merge)


if __name__ == "__main__":
    main()
