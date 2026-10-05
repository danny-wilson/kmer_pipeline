#!/usr/bin/env python3
"""plotManhattan.py: QQ and Manhattan plots, top genes, k-mer tables, and the
close-up alignment figures of the top genes. Port of plotManhattan.Rscript.
Writes <prefix>_<type><k>.summary.json for the reports (PLAN 5.3)."""
import argparse
import math
import os
import sys
import time

import numpy as np

import rcompat
from rcompat import r_cat, r_paste0, r_stop

###################################################################################################
## Functions and software paths
###################################################################################################


def read_ref_length(ref_gb):
    """The genome length: the third word of the GenBank file's first line."""
    with rcompat.r_open(ref_gb) as f:
        first = f.readline().rstrip("\n")
    toks = [t for t in first.split(" ") if t != ""]
    ref_length = rcompat.r_as_numeric(toks[2]) if len(toks) >= 3 else None
    if ref_length is None:
        r_stop("Error retrieving the reference genome length from the genbank file", "\n")
    return ref_length


def first_gene_ids(ref, lookup_rows, ref_length):
    """ref.pos.gene.id[[pos]][1] from create_gene_lookup: for each reference
    position, the first gene (in reference order) covering it, else the
    intergenic region; 0 for none. As an array indexed by position."""
    import numpy as np
    n = int(ref_length)
    first = np.zeros(n + 2, dtype=np.int64)
    nref = len(ref)
    for gid, (name, _, s, e, _) in enumerate(lookup_rows, start=1):
        lo, hi = int(min(s, e)), int(max(s, e))
        if gid > nref and gid == len(lookup_rows):  # the final region wraps round to before the first gene
            ranges = [(lo, n), (1, int(lookup_rows[0][2]) - 1)]
        else:
            ranges = [(lo, hi)]
        for a, b in ranges:
            a, b = max(a, 1), min(b, n)
            if a <= b:
                seg = first[a:b + 1]
                seg[seg == 0] = gid
    return first


def get_gene_xpos(gene_lookup_file, ref_gb):
    import sequence_functions
    ref_length = read_ref_length(ref_gb)
    r_cat("Reference genome length:", ref_length, "\n")

    # Read in gene lookup
    gene_lookup = rcompat.r_read_table(gene_lookup_file, header=False, sep="\t", quote="")
    names_lookup = [rcompat.r_as_character(v) for v in gene_lookup.iloc[:, 0]]

    # Read in reference genbank file
    ref = sequence_functions.reorder_reference_gbk(ref_gb=ref_gb)
    ref_names = list(ref["name"])
    ref_start = [float(v) for v in ref["start"]]
    ref_end = [float(v) for v in ref["end"]]

    def where(name):
        return [k for k, n in enumerate(ref_names) if n == name]

    # Get the midpoint for all genes/intergenic regions to plot them at on the Manhattan
    gene_lookup_pos = []
    n = len(names_lookup)
    for i in range(n):
        name = names_lookup[i]
        if i != n - 1:
            if ":" in name:
                genes = name.split(":")
                w1, w2 = where(genes[0]), where(genes[1] if len(genes) > 1 else None)
                start = [ref_end[k] + 1 for k in w1]
                end = [ref_start[k] - 1 for k in w2]
            else:
                w = where(name)
                start = [ref_start[k] for k in w]
                end = [ref_end[k] for k in w]
        else:
            gene = [g for g in name.split(":") if g != ""]
            w = [k for g in gene for k in where(g)]
            start = [ref_start[k] for k in w]
            end = [ref_length]
        if len(start) != 1 or len(end) != 1:
            r_stop("Error in gene_lookup_pos[i] = ...: no single reference position for ", name)
        gene_lookup_pos.append(((end[0] - start[0]) / 2) + start[0])
    r_cat("Read in reference genbank file", "\n")
    return {"gene_lookup": names_lookup, "ref": ref, "gene_lookup_pos": gene_lookup_pos, "ref_length": ref_length}


def read_kmer_alignment(alignPosFile, alignCountFile, min_count, gene_lookup_pos, gene_lookup, kmerIndex, ref_length):
    lines = rcompat.r_scan_lines(alignPosFile, quiet=True)
    pairs = [l.split(",") for l in lines]
    if any(len(p) != 2 for p in pairs):
        r_stop("Error: Align pos matrix is not two rows", "\n")
    count = []
    for tok in rcompat.r_scan_lines(alignCountFile, quiet=True):
        v = rcompat.r_as_numeric(tok)
        count.append(math.nan if v is None else v)
    keep = [k for k in range(len(count)) if count[k] >= min_count]
    pairs = [pairs[k] for k in keep if k < len(pairs)]
    count = [count[k] for k in keep]
    # R: rep(1, length(alignPos)) where alignPos is the 2 x n matrix, so twice as long
    alignPosPCH = [1.0] * (2 * len(pairs))
    for k, c in enumerate(count):
        if c == 1:
            alignPosPCH[k] = 2.0
        elif 1 < c <= 5:
            alignPosPCH[k] = 0.0

    final_kmer_pos_index = [float(p[0]) for p in pairs]
    final_kmer_pos = [gene_lookup_pos[int(float(p[1])) - 1] for p in pairs]
    final_kmer_genes = [gene_lookup[int(float(p[1])) - 1] for p in pairs]
    present = set(final_kmer_pos_index)
    missing = [k for k in range(1, len(kmerIndex) + 1) if float(k) not in present]
    if missing:
        final_kmer_pos_index += [float(k) for k in missing]
        final_kmer_pos += [ref_length + 100000 + k * 0.01 for k in range(len(missing))]
        final_kmer_genes += [None] * len(missing)
        alignPosPCH += [1.0] * len(missing)
    r_cat("Got final kmer positions", "\n")
    return {"alignPosPCH": alignPosPCH, "final_kmer_pos_index": final_kmer_pos_index, "final_kmer_pos": final_kmer_pos,
            "final_kmer_genes": final_kmer_genes}


###################################################################################################


def main():
    rcompat.script_setup(__file__)
    parser = argparse.ArgumentParser(description="plotManhattan.py plot QQ and Manhattan plots", allow_abbrev=False)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--analysis-dir", required=True)
    parser.add_argument("--kmerfile-prefix", required=True)
    parser.add_argument("--ref-gb", required=True)
    parser.add_argument("--ref-fa", required=True)
    parser.add_argument("--gene-lookup-file", required=True)
    parser.add_argument("--id-file", required=True)
    parser.add_argument("--nucmerident", required=True)
    parser.add_argument("--plot-min-genomes", "--min-count", dest="min_count", required=True,
                        help="genomes a k-mer/gene combination must be seen in to be plotted (--min-count "
                             "is the old name)")
    parser.add_argument("--kmer-type", required=True)
    parser.add_argument("--kmer-length", required=True)
    parser.add_argument("--minor-allele-threshold", required=True)
    parser.add_argument("--software-file", required=True)
    parser.add_argument("--blastident", required=True)
    parser.add_argument("--ngenes", required=True)
    parser.add_argument("--annotate-gene-file", default=None, help="genes/IRs to annotate (instead of the top genes)")
    parser.add_argument("--override-signif", default="FALSE", help="TRUE/FALSE: plot alignments whatever the significance")
    run(parser.parse_args(), bowtie=False)


def run(args, bowtie):
    """The body of plotManhattan; plotManhattanbowtie.py runs it with bowtie=True
    (k-mer positions from the bowtie2 mapping instead of the contig alignment)."""
    start_time = time.monotonic()

    # Initialize variables
    output_prefix = args.output_prefix
    output_dir = args.analysis_dir
    kmerfilePrefix = args.kmerfile_prefix
    ref_gb = args.ref_gb
    ref_fa = args.ref_fa
    id_file = args.id_file
    if bowtie:
        gene_lookup_file = nucmerident = min_count = None
        samtools_filter = rcompat.r_as_integer(args.samtools_filter)
    else:
        gene_lookup_file = args.gene_lookup_file
        nucmerident = rcompat.r_as_integer(args.nucmerident)
        min_count = rcompat.r_as_integer(args.min_count)
    kmer_type = args.kmer_type.lower()
    kmer_length = rcompat.r_as_integer(args.kmer_length)
    minor_allele_threshold = rcompat.r_as_numeric(args.minor_allele_threshold)
    software_file = args.software_file
    blastident = rcompat.r_as_integer(args.blastident)
    ngenes = rcompat.r_as_integer(args.ngenes)
    annotateGeneFile = args.annotate_gene_file
    logical = {"TRUE": True, "T": True, "true": True, "True": True, "FALSE": False, "F": False, "false": False,
               "False": False}
    override_signif = logical.get(args.override_signif) if annotateGeneFile is not None else False

    # Check file inputs
    if not os.path.exists(output_dir):
        r_stop("Error: output directory doesn't exist", "\n")
    if not output_dir.endswith("/"):
        output_dir = output_dir + "/"
    if not os.path.exists(ref_gb):
        r_stop("Error: reference genbank file doesn't exist", "\n")
    if not os.path.exists(ref_fa):
        r_stop("Error: reference fasta file doesn't exist", "\n")
    if not bowtie and not os.path.exists(gene_lookup_file):
        r_stop("Error: reference gene ID file doesn't exist", "\n")
    if not os.path.exists(id_file):
        r_stop("Error: sample ID file doesn't exist", "\n")
    if not bowtie and (nucmerident is None or nucmerident > 100 or nucmerident < 0):
        r_stop("Error: nucmer identity threshold must be between 0-100", "\n")
    if not bowtie and min_count is None:
        r_stop("Error: min count must be an integer", "\n")
    if kmer_type != "protein" and kmer_type != "nucleotide":
        r_stop("Error: kmer type must be either 'protein' or 'nucleotide'", "\n")
    if kmer_length is None:
        r_stop("Error: kmer length must be an integer", "\n")
    if minor_allele_threshold is None:
        r_stop("Error: minor allele threshold must be a number", "\n")
    if not bowtie and 0.5 < minor_allele_threshold < 1:
        r_stop("Error: minor allele threshold must be <=0.5 or >=1")
    if not os.path.exists(software_file):
        r_stop("Error: software file doesn't exist", "\n")
    if blastident is None or blastident > 100 or blastident < 0:
        r_stop("Error: blast identity threshold must be between 0-100", "\n")
    if annotateGeneFile is not None and not os.path.exists(annotateGeneFile):
        r_stop("Error: annotate gene file doesn't exist", "\n")
    if override_signif is None:
        r_stop("Error: override_signif must be a logical", "\n")

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import kmercontigalignonly
    ref_name = kmercontigalignonly.read_reference_name(ref_fa)

    kmerKeySizeFile = kmerfilePrefix + ".patternmerge.patternKeySize.txt"
    kmerIndexFile = kmerfilePrefix + ".patternmerge.patternIndex.txt.gz"
    kmerPresenceCountFile = kmerfilePrefix + ".patternmerge.presenceCount.txt.gz"
    kmerSeqFile = kmerfilePrefix + ".kmermerge.txt.gz"
    if bowtie:
        mappingFile = r_paste0(kmerfilePrefix, ".", ref_name, ".SAMq", samtools_filter, ".bowtie2map.txt.gz")
        inputs = ((kmerKeySizeFile, "kmer pattern key size"), (kmerIndexFile, "kmer pattern index"),
                  (kmerPresenceCountFile, "kmer presence count"), (mappingFile, "bowtie2 mapping"),
                  (kmerSeqFile, "kmer sequence"))
    else:
        alignPosFile = r_paste0(kmerfilePrefix, ".", ref_name, "_t", nucmerident, ".kmeralignmerge.txt.gz")
        alignCountFile = r_paste0(kmerfilePrefix, ".", ref_name, "_t", nucmerident, ".kmeralignmerge.count.txt.gz")
        inputs = ((kmerKeySizeFile, "kmer pattern key size"), (kmerIndexFile, "kmer pattern index"),
                  (kmerPresenceCountFile, "kmer presence count"), (alignPosFile, "align pos"),
                  (alignCountFile, "align count"), (kmerSeqFile, "kmer sequence"))
    for f, what in inputs:
        if not os.path.exists(f):
            r_stop("Error: " + what + " file doesn't exist" + (": " + f + " \n" if bowtie else ""), "" if bowtie else "\n")

    # Read in software file
    software_paths = rcompat.r_read_table(software_file, header=True, sep="\t", quote="")
    names = [rcompat.r_as_character(v) for v in software_paths["name"]]
    paths = [rcompat.r_as_character(v) for v in software_paths["path"]]
    # R draws the figures (plot_figures.R, run by Nextflow after this step), so its entry is required
    required_software = ["scriptpath", "genoPlotR", "blast", "R"]
    if any(r not in names for r in required_software):
        r_stop("Error: missing required software path in the software file - requires " + ", ".join(required_software), "\n")

    def software(name):
        return [pth for nm, pth in zip(names, paths) if nm.lower() == name][0]
    script_location = software("scriptpath")
    if not os.path.isdir(script_location):
        r_stop("Error: script location directory specified in the software paths file doesn't exist", "\n")
    for mod in ("sequence_functions.py", "Manhattan_functions.py", "alignmentfunctions.py"):
        if not os.path.exists(script_location + "/" + mod):
            r_stop("Error: " + mod + " path doesn't exist - check pipeline script location in the software file", "\n")
    sys.path.insert(0, script_location)
    import Manhattan_functions as mf
    import sequence_functions

    Rscriptpath = software("r") + "script"
    if not os.path.exists(Rscriptpath):
        r_stop("Error: Rscript path ", Rscriptpath, " doesn't exist (R draws the figures)", "\n")
    if not os.path.exists(script_location + "/plot_figures.R"):
        r_stop("Error: plot_figures.R path doesn't exist - check pipeline script location in the software file", "\n")

    blast_dir = software("blast")
    if not os.path.isdir(blast_dir):
        r_stop("Error: blast installation directory specified in the software paths file doesn't exist", "\n")
    blastname = "blastp" if kmer_type == "protein" else "blastn"
    blastPath = blast_dir + "/" + blastname
    if not os.path.exists(blastPath):
        r_stop("Error: blast path", blastPath, "doesn't exist", "\n")

    # Report variables
    r_cat("#############################################", "\n")
    r_cat("Running on host: ", rcompat.r_system_intern("hostname"), "\n")
    r_cat("Command line arguments", "\n")
    r_cat(sys.argv[1:], "\n\n")
    r_cat("Parameters:", "\n")
    for label, v in (("Output prefix:", output_prefix), ("Analysis directory:", output_dir),
                     ("Kmer file prefix:", kmerfilePrefix), ("Kmer pattern key size file:", kmerKeySizeFile),
                     ("Kmer pattern index file:", kmerIndexFile), ("Kmer pattern presence count file:", kmerPresenceCountFile),
                     ("Kmer list file:", kmerSeqFile), ("Reference genbank file:", ref_gb),
                     ("Reference fasta file:", ref_fa), ("ID file path:", id_file), ("Kmer type:", kmer_type),
                     ("Kmer length:", kmer_length), ("Minor allele threshold:", minor_allele_threshold),
                     ("BLAST alignment minimum % identity:", blastident), ("Number of top genes to output:", ngenes),
                     ("Software file:", software_file), ("Script location:", script_location)):
        r_cat(label, v, "\n")
    if annotateGeneFile is not None:
        r_cat("Annotate gene file:", annotateGeneFile, "\n")
        r_cat("Override significance threshold:", override_signif, "\n")
    r_cat("#############################################", "\n\n")

    # Create an output directory
    alignmenttype = "bowtie2mapping" if bowtie else "kmergenealign"
    figures_dir = mf.create_figures_dir(dir=output_dir, kmer_type=kmer_type, kmer_length=kmer_length,
                                        alignmenttype=alignmenttype)

    # Get GEMMA input directory
    gemma_dir = output_dir + "/" + r_paste0(kmer_type, "kmer", kmer_length, "_gemma") + "/" + "output/"
    if not os.path.isdir(gemma_dir):
        r_stop("Error: gemma directory", gemma_dir, " doesn't exist", "\n")

    if minor_allele_threshold == 0:
        r_cat("Assuming no minor allele threshold - plotting results for all kmers", "\n")
        macormaf = "maf"
    elif minor_allele_threshold < 1:
        r_cat("Minor allele threshold below 1 - reading as a minor allele frequency (MAF) threshold of:",
              minor_allele_threshold, "\n")
        macormaf = "maf"
    else:
        r_cat("Minor allele threshold above 1 - reading as a minor allele count (MAC) threshold of:",
              minor_allele_threshold, "\n")
        macormaf = "mac"

    # Read in ID file
    id_table = rcompat.r_read_table(id_file, header=True, sep="\t")
    ids = [rcompat.r_as_character(v) for v in id_table["id"]]
    pheno = [rcompat.r_as_numeric_value(v) for v in id_table["pheno"]]
    nsamples = sum(1 for p in pheno if p is not None)

    # Check count threshold variable
    if not bowtie and (min_count < 1 or min_count > len(ids)):
        r_stop("Error: minimum count must be at least 1 and less than the total number of samples", "\n")

    # Read in total number of kmer patterns and the index
    nPatterns = float(open(kmerKeySizeFile).read().split()[0])
    kmerIndex = [int(float(t)) + 1 for t in rcompat.r_open(kmerIndexFile).read().split()]

    r_cat("Read in number of patterns and kmer index", "\n")
    r_cat("Number of kmers:", len(kmerIndex), "\n")
    r_cat("Number of patterns:", nPatterns, "\n")
    if len(set(kmerIndex)) != nPatterns:
        r_stop("Error: number of unique kmer indices does not equal the number of patterns", "\n")

    # Read in reference and gene look up
    if bowtie:
        ref_length = read_ref_length(ref_gb)
        print("Read 1 item", file=sys.stderr, flush=True)
        r_cat("Reference genome length:", ref_length, "\n")
        ref = sequence_functions.reorder_reference_gbk(ref_gb=ref_gb)
        r_cat("Read in reference genbank file", "\n")
        import alignmentfunctions
        lookup_rows = alignmentfunctions.create_gene_lookup(ref, ref_length)
        gene_lookup = [row[0] for row in lookup_rows]
        first_id = first_gene_ids(ref, lookup_rows, ref_length)
    else:
        gene_xpos = get_gene_xpos(gene_lookup_file=gene_lookup_file, ref_gb=ref_gb)
        ref_length = gene_xpos["ref_length"]
        ref = gene_xpos["ref"]
        gene_lookup = gene_xpos["gene_lookup"]
        gene_lookup_pos = gene_xpos["gene_lookup_pos"]

    # Read in gemma files
    assoc = mf.read_gemma_files(input_dir=gemma_dir, prefix=output_prefix, kmer_type=kmer_type, kmer_length=kmer_length,
                                nPatterns=nPatterns)
    neglog10 = mf.assoc_column(assoc, 6)
    beta_patterns = mf.assoc_column(assoc, 2)

    ## Read in MAF
    macpatterns = []
    for tok in rcompat.r_scan_lines(kmerPresenceCountFile, quiet=True):
        v = rcompat.r_as_numeric(tok)
        macpatterns.append(math.nan if v is None else v)
    macpatterns = np.array(macpatterns, dtype=float)
    length_nonNApheno = nsamples
    if any(p is None for p in pheno):
        r_cat(r_paste0("Calculating MAFs using number of samples with non NA phenotypes (", length_nonNApheno,
                       ") as denominator"), "\n")
    with np.errstate(invalid="ignore"):
        hi = macpatterns > (length_nonNApheno / 2)
    macpatterns[hi] = length_nonNApheno - macpatterns[hi]
    ki = np.array(kmerIndex) - 1
    mac = macpatterns[ki]
    mafpatterns = macpatterns / length_nonNApheno
    maf = mafpatterns[ki]

    r_cat("Read in pattern counts and converted into MAC and MAF", "\n")
    if macormaf == "maf":
        mapatterns, ma = mafpatterns, maf
    else:
        mapatterns, ma = macpatterns, mac

    ## Get Bonferroni threshold
    r_cat("Bonferroni threshold calculated using", macormaf, "threshold", minor_allele_threshold, "\n")
    with np.errstate(invalid="ignore"):
        tested = [kmerIndex[k] for k in range(len(kmerIndex)) if ma[k] >= minor_allele_threshold and assoc[kmerIndex[k] - 1] is not None]
    n_tests = len(set(tested))
    bonferroni = -math.log10(0.05 / n_tests)
    r_cat("Bonferroni threshold:", bonferroni, "\n")
    mf.write_summary_json(summary_file=r_paste0(output_dir, output_prefix, "_", kmer_type, kmer_length,
                                                ".bowtie2mapping.summary.json" if bowtie else ".summary.json"),
                          n_kmers=len(kmerIndex), n_patterns=nPatterns,
                          n_untested_patterns=sum(1 for r in assoc if r is None),
                          max_neglog10p=float(np.nanmax(neglog10)), minor_allele_threshold=minor_allele_threshold,
                          macormaf=macormaf, n_tests=n_tests, bonferroni=bonferroni)

    ## Figure data for plot_figures.R (PLAN 5.5): QQ plots
    fd = mf.FigureData(figures_dir)
    fd.table("patterns", [("neglog10p", "numeric"), ("beta", "numeric"), ("maf", "numeric"), ("ma", "numeric")],
             zip(neglog10, beta_patterns, mafpatterns, mapatterns))
    fd.table("kmers", [("kmer_index", "integer"), ("ma", "numeric")], zip(kmerIndex, ma))
    stem_qq = figures_dir + output_prefix + "_" + kmer_type + rcompat.r_as_character(kmer_length)
    for thr in (0.0, minor_allele_threshold):
        fd.expect(stem_qq + ("_QQplot_allkmers.png" if thr == 0 else r_paste0("_QQplot_", macormaf, thr, ".png")))

    ## Read in alignment results
    if bowtie:
        import plotManhattanbowtie
        km = plotManhattanbowtie.read_bowtie_pos(mappingFile, kmerIndex, ref_length)
        final_kmer_pos_index = km["final_kmer_pos_index"]
        final_kmer_pos = km["final_kmer_pos"]
        lookup_start1 = lookup_rows[0][2]
        final_kmer_genes = [None if pos > ref_length or pos < lookup_start1 else
                            (gene_lookup[first_id[int(pos)] - 1] if first_id[int(pos)] > 0 else None)
                            for pos in final_kmer_pos]
        r_cat("Assigned genes/IRs to each kmer", "\n")
    else:
        ka = read_kmer_alignment(alignPosFile, alignCountFile, min_count, gene_lookup_pos, gene_lookup, kmerIndex,
                                 ref_length)
        final_kmer_pos_index = ka["final_kmer_pos_index"]
        final_kmer_pos = ka["final_kmer_pos"]
        final_kmer_genes = ka["final_kmer_genes"]

    ## Get y position
    fki = np.array([int(v) for v in final_kmer_pos_index]) - 1
    ypos = neglog10[ki[fki]]
    r_cat("Got ypos", "\n")

    pheno_type = mf.get_pheno_type(pheno)

    # Top genes are chosen from all kmers, before any subsampling for plotting
    ma_full = ma[fki]
    gene_conversion_full = {g: g for g in final_kmer_genes if g is not None}
    mf.top20genes(final_kmer_genes, ma_full, minor_allele_threshold, ypos, macormaf, figures_dir, output_prefix, min_count,
                  nucmerident, kmer_type, kmer_length, ref_name)
    which_genes_full = [k for k in range(len(final_kmer_genes))
                        if final_kmer_genes[k] is not None and ma_full[k] >= minor_allele_threshold]
    rows = mf.get_genes_to_plot([final_kmer_genes[k] for k in which_genes_full], [ypos[k] for k in which_genes_full],
                                gene_conversion_full, [10, 0], [], ref, ngenes=ngenes)
    topngenes = [r[0] for r in rows]

    ## Figure data: genome-wide Manhattan plots (drawn, coloured and subsampled by plot_figures.R)
    fd.table("positions", [("kmer", "numeric"), ("position", "numeric"), ("gene", "character")],
             zip(final_kmer_pos_index, final_kmer_pos, final_kmer_genes))
    fd.table("reference_cds", [("name", "character"), ("start", "numeric"), ("end", "numeric"), ("strand", "numeric")],
             zip(ref["name"], ref["start"], ref["end"], ref["strand"]))
    if annotateGeneFile is not None:
        with open(fd.dir + "annotate_genes.txt", "w") as f:
            f.write("".join(g + "\n" for g in rcompat.r_scan_lines(annotateGeneFile, quiet=True)))
    if bowtie:
        manhattan_stem = r_paste0(figures_dir, output_prefix, "_", kmer_type, kmer_length, "_", ref_name,
                                  "_LMM_bowtie2mapping")
    else:
        manhattan_stem = r_paste0(figures_dir, output_prefix, "_", kmer_type, kmer_length, "_", ref_name,
                                  "_LMM_kmergenealign_ct", min_count)
    for key, value in (("figures_dir", figures_dir), ("output_prefix", output_prefix), ("kmer_type", kmer_type),
                       ("kmer_length", kmer_length), ("ref_name", ref_name), ("ref_length", float(ref_length)),
                       ("macormaf", macormaf), ("minor_allele_threshold", float(minor_allele_threshold)),
                       ("bonferroni", float(bonferroni)), ("pheno_type", pheno_type), ("nsamples", nsamples),
                       ("override_signif", bool(override_signif)), ("manhattan_stem", manhattan_stem)):
        fd.param(key, value)
    # The Manhattan plots R draws: four colourings, and a ylim 50 version when the maximum is at least 75
    finite = ypos[~np.isnan(ypos)]
    ylim50 = len(finite) > 0 and not (finite.max() < 50 + 50 / 2)
    for col, thr in (("alignCOL", minor_allele_threshold), ("betaCOL", minor_allele_threshold), ("mafCOL", 0.0),
                     ("mafCOL", minor_allele_threshold)):
        prefix_i = r_paste0(manhattan_stem, "_Manhattan_", col, "_", macormaf, thr)
        fd.expect(prefix_i + ".png")
        if ylim50:
            fd.expect(prefix_i + "_ylim50.png")

    final_kmer_list = rcompat.r_scan_lines(kmerSeqFile, quiet=True)

    r_cat("Writing unaligned significant kmers to file", "\n")
    if bowtie:
        output_file_prefix = r_paste0(figures_dir, output_prefix, "_", kmer_type, kmer_length, "_", ref_name, "_",
                                      macormaf, "_", minor_allele_threshold, "_bowtie2mapping")
    else:
        output_file_prefix = r_paste0(figures_dir, output_prefix, "_", kmer_type, kmer_length, "_", ref_name, "_",
                                      macormaf, "_", minor_allele_threshold, "_alignIdent_", nucmerident,
                                      "_alignPosMinCount_", min_count)
    wh_i = [k for k in range(len(final_kmer_pos)) if final_kmer_pos[k] > ref_length]
    mf.write_top_gene_kmers_to_file(wh_i, final_kmer_list, final_kmer_pos_index, assoc, kmerIndex, mac,
                                    output_file_prefix + "_unaligned_kmersandpvals.txt")

    if annotateGeneFile is not None:
        r_cat("Writing kmers for gene names in", annotateGeneFile, "to file", "\n")
        annotateGene = rcompat.r_scan_lines(annotateGeneFile, quiet=True)
        outputfiles = [output_file_prefix + "_namedgene_" + str(k + 1) + "_" + g + "_kmersandpvals.txt"
                       for k, g in enumerate(annotateGene)]
        gene_list = annotateGene
    else:
        r_cat(rcompat.r_paste("Writing kmers for top", ngenes, "genes to file"), "\n")
        outputfiles = [output_file_prefix + "_topgene_" + str(k + 1) + "_" + g + "_kmersandpvals.txt"
                       for k, g in enumerate(topngenes)]
        gene_list = topngenes
    for g, out in zip(gene_list, outputfiles):
        wh_i = [k for k in range(len(final_kmer_genes)) if final_kmer_genes[k] == g]
        mf.write_top_gene_kmers_to_file(wh_i, final_kmer_list, final_kmer_pos_index, assoc, kmerIndex, mac, out)
    genes_all = {"genes": outputfiles, "genes_names": gene_list}
    del final_kmer_list

    ## Plot close up alignments
    import alignmentfunctions
    alignmentfunctions.plot_closeup_alignments(
        ref=ref, ref_length=ref_length, ref_gb=ref_gb, ref_fa=ref_fa, figures_dir=figures_dir,
        output_prefix=output_prefix, ngenes=ngenes, nsamples=nsamples, bonferroni=bonferroni, gene_lookup=gene_lookup,
        oneLetterCodes=alignmentfunctions.oneLetterCodes, kmer_type=kmer_type, kmer_length=kmer_length,
        blastPath=blastPath, perident=blastident, ref_name=ref_name, alignmenttype=alignmenttype,
        override_signif=override_signif, genes_all=genes_all, minor_allele_threshold=minor_allele_threshold,
        macormaf=macormaf, figure_data=fd)
    fd.close()
    r_cat("Written figure data for plot_figures.R:", fd.dir, "\n")

    r_cat("Finished in", (time.monotonic() - start_time) / 60, "minutes\n")


if __name__ == "__main__":
    main()
