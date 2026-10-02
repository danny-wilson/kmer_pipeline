#!/usr/bin/env python3
"""plotManhattanbowtie.py: QQ and Manhattan plots, top genes and close-up
alignment figures, with k-mer positions from the bowtie2 mapping (runbowtie.py)
instead of the contig alignment. Port of plotManhattanbowtie.Rscript; the body
shared with plotManhattan is plotManhattan.run(). Not called by kmer_pipeline.nf."""
import argparse
import os
import sys

import rcompat
from rcompat import r_cat, r_stop

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plotManhattan  # noqa: E402


def read_bowtie_pos(mappingFile, kmerIndex, ref_length):
    """Positions of the mapped k-mers (SAM column 4) and their indices (column 1,
    0-based), followed by the unmapped k-mers placed beyond the genome."""
    bowtie_pos = [float(x) for x in rcompat.r_system_intern("zcat " + mappingFile + " | cut -f4")]
    bowtie_index = [float(x) + 1 for x in rcompat.r_system_intern("zcat " + mappingFile + " | cut -f1")]
    present = set(bowtie_index)
    missing = [float(k) for k in range(1, len(kmerIndex) + 1) if float(k) not in present]
    if not missing:  # R only defines the results when some k-mers are unmapped
        r_stop("Error: object 'final_kmer_pos_index' not found (every k-mer mapped)")
    final_kmer_pos_index = bowtie_index + missing
    final_kmer_pos = bowtie_pos + [ref_length + 100000 + k * 0.01 for k in range(len(missing))]
    r_cat("Got final kmer positions", "\n")
    return {"final_kmer_pos_index": final_kmer_pos_index, "final_kmer_pos": final_kmer_pos}


def main():
    rcompat.script_setup(__file__)
    parser = argparse.ArgumentParser(description="plotManhattanbowtie.py plot QQ and Manhattan plots", allow_abbrev=False)
    for name in ("output-prefix", "analysis-dir", "kmerfile-prefix", "ref-gb", "ref-fa", "id-file", "kmer-type",
                 "kmer-length", "minor-allele-threshold", "samtools-filter", "software-file", "blastident", "ngenes"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--annotate-gene-file", default=None, help="genes/IRs to annotate (instead of the top genes)")
    parser.add_argument("--override-signif", default="FALSE", help="TRUE/FALSE: plot alignments whatever the significance")
    plotManhattan.run(parser.parse_args(), bowtie=True)


if __name__ == "__main__":
    main()
