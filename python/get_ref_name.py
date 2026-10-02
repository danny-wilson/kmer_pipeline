#!/usr/bin/env python3
"""get_ref_name.py: read the name of the supplied reference genome.
Port of get_ref_name.Rscript. Prints the first word of the FASTA header, without
the '>', to stdout with no newline. kmer_pipeline.nf requires stderr to be empty."""
import argparse
import sys

import rcompat


def main():
    rcompat.script_setup(__file__, announce=False)
    parser = argparse.ArgumentParser(
        description="Read the name of the supplied reference genome. Daniel Wilson (2022)",
        allow_abbrev=False)
    parser.add_argument("--fasta-file", required=True, help="reference genome in FASTA format")
    args = parser.parse_args()
    REF_FA = args.fasta_file

    with rcompat.r_open(REF_FA) as f:
        text = f.read()
    lines = text.split("\n")
    if lines[-1] == "":
        lines.pop()
    # scan(nlines = 1) gives nothing when the first line is blank
    ref_name = lines[0] if lines and lines[0] != "" else None
    if sum(line.startswith(">") for line in lines) > 1:
        raise RuntimeError(f"Error: reference fasta file {REF_FA} contains more than one record; "
                           "only single-record references are supported")
    if ref_name is not None:
        ref_name = ref_name.split(" ")[0]
        ref_name.encode("utf-8")  # fails on bytes that are not UTF-8, as R's substr() does
    if ref_name is None or ref_name[:1] != ">":
        raise RuntimeError("Error: reference fasta file does not begin with a name starting with >")
    ref_name = ref_name[1:1000000]
    sys.stdout.write(ref_name)


if __name__ == "__main__":
    main()
