#!/usr/bin/env python3
"""gen-protein-report.py: generate the k-mer GWAS HTML report for one of the top
genes (protein k-mers). Port of gen-protein-report.Rscript (Daniel Wilson, 2022).
The HTML is built line by line as the R script builds it; the opening of the
report is shared with gen-gene-report.py, as the two R scripts share it."""
import argparse
import math
import os
import sys

import numpy as np

import rcompat
from rcompat import r_s3 as s3

NL = "\n"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
gr = __import__("gen-report")
ggr = __import__("gen-gene-report")


def match_first(values, vec):
    """match(values, vec) for a 1-based vec given as a dict value -> first index."""
    return [None if v is None else vec.get(v) for v in values]


def coord_frame(ref_length, lo, hi, forward):
    """gene.coords.in.ref[, fm]: coordinates 1, 2, ... every third base from lo
    (forward) or hi (reverse), recycled as R recycles 1:((hi-lo+1)/3) over the
    longer seq(). Returned as value -> first reference position."""
    positions = list(range(lo, hi + 1, 3)) if forward else list(range(hi, lo - 1, -3))
    n = int(math.floor((hi - lo + 1) / 3))
    vals = list(range(1, n + 1)) if n >= 1 else [1, 0]
    first = {}
    for k, pos in enumerate(positions):
        v = vals[k % len(vals)]
        first[v] = min(first.get(v, pos), pos)
    return first


def main():
    rcompat.script_setup(__file__)
    parser = argparse.ArgumentParser(description="gen-protein-report.py Generate a kmer GWAS report for a specific "
                                                 "protein. Daniel Wilson (2022)", allow_abbrev=False)
    for name in ("hit-num", "prefix", "anatype", "k", "refname", "ref-gb", "maf", "alignident", "mincount", "srcdir",
                 "outdir", "logdir"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    out = ggr.run(args, protein=True)
    if out is None:
        return
    outfile_html, html_head, html_body, html_foot, v = out
    FIGDIR, stem, gene, gene_html, MACORMAF, MAF = v["FIGDIR"], v["stem"], v["gene"], v["gene_html"], v["MACORMAF"], v["MAF"]
    is_maf_text3, thr_signif, is_intergenic, genes = v["is_maf_text3"], v["thr_signif"], v["is_intergenic"], v["genes"]
    gbk, ref_length, first = v["gbk"], v["ref_length"], v["first"]

    def rdir(pattern):
        d = os.path.dirname(pattern)
        return rcompat.r_dir(d, glob=os.path.basename(pattern), full_names=True) if os.path.isdir(d) else []

    cf_maf = rdir(FIGDIR + stem + "_" + gene + "_correct_frame_*_Manhattan_" + MACORMAF + MAF + ".png")
    cf_maf0 = rdir(FIGDIR + stem + "_" + gene + "_correct_frame_*_Manhattan_allkmers.png")
    filenames_Manhattan = [cf_maf[0] if cf_maf else "NA", cf_maf0[0] if cf_maf0 else "NA",
                           FIGDIR + stem + "_" + gene + "_allframes_Manhattan_" + MACORMAF + MAF + ".png",
                           FIGDIR + stem + "_" + gene + "_allframes_Manhattan_allkmers.png"]
    filename_kmer_maf = (FIGDIR + stem + "_" + gene + "_correct_frame_*_plot_*_aminoacids_*_to_*_" + MACORMAF + MAF
                         + "_alignment.png")
    filename_kmer_maf0 = FIGDIR + stem + "_" + gene + "_correct_frame_*_plot_*_aminoacids_*_to_*_alignment.png"
    filename_mapped = FIGDIR + stem + "_top_gene_*_" + gene + "_*_frame_blast_results.txt"
    filename_unmapped = FIGDIR + stem + "_top_gene_*_" + gene + "_no_blast_result_or_poor_alignment.txt"

    # In R the statement ends after "or by frame." (no %/%): the paragraph on frames that follows
    # it in the source is evaluated separately and not added to the report
    html_body = NL.join([
        html_body,
        "  <h2>Manhattan plot for " + gene_html + " </h2>",
        "  <p>The Figure displays the significance of each kmer against the position in the",
        "  reference genome to which it mapped, with a focus on " + gene_html + " .",
        "  The Bonferroni-corrected significance threshold is shown as a horizontal black dashed line.",
        "  Annotated features are plotted below. Points are shaded by direction of association",
        "  (in-frame kmers only) or by frame."])
    t3 = is_maf_text3
    html_body = gr.slideshow(html_body, filenames_Manhattan,
                             ["Kmers mapping to the region, in the correct frame, filtered by " + t3 + ".",
                              "Kmers mapping to the region, in the correct frame. No " + t3 + " filter.",
                              "Kmers mapping to the region, in any frame, filtered by " + t3 + ".",
                              "Kmers mapping to the region, in any frame. No " + t3 + " filter."], "")
    html_body = NL.join([html_body, NL])

    filenames_mapped = rdir(filename_mapped)
    frames = []
    for s in filenames_mapped:
        s2 = s.split("_")
        frames.append(s2[len(s2) - 5] + " " + s2[len(s2) - 4])
    filenames_unmapped = rdir(filename_unmapped)
    blast_nul_exists = len(filenames_unmapped) > 0
    if not filenames_mapped:
        rcompat.r_stop("Error in file(file, \"rt\"): invalid 'description' argument (no BLAST results for ", gene, ")")
    blast_map = []
    for fn, fr in zip(filenames_mapped, frames):
        for r in ggr.read_delim(fn)[1]:
            r["frame"] = fr
            blast_map.append(r)
    if blast_nul_exists:
        blast_nul = []
        for fn in filenames_unmapped:
            blast_nul += ggr.read_delim(fn)[1]
        mapped = {r["kmer"] for r in blast_map}
        blast_nul_gd = [r["kmer"] not in mapped and r["negLog10"] >= thr_signif for r in blast_nul]

    filenames_kmer_maf = rdir(filename_kmer_maf)
    maf_set = set(filenames_kmer_maf)
    filenames_kmer_maf0 = [f for f in rdir(filename_kmer_maf0) if f not in maf_set]

    if filenames_kmer_maf:
        pstem = FIGDIR + stem + "_" + gene + "_correct_frame_"  # a regular expression, as in R

        def f(s):
            import re
            parts = [rcompat.r_as_numeric(p) for p in re.sub(pstem, "", s).split("_")]
            get = lambda k: parts[k] if k < len(parts) else None  # noqa: E731
            return get(2), get(4), get(6)
        win_maf = [f(s) for s in filenames_kmer_maf]
        win_maf0 = [f(s) for s in filenames_kmer_maf0]
        if win_maf != win_maf0:
            rcompat.r_stop("Error: all(unname(win.maf) == unname(win.maf0)) is not TRUE")

        gstart = [float(x) for x in gbk["start"]]
        gend = [float(x) for x in gbk["end"]]
        coords = {}
        if is_intergenic:
            idx = [first.get(g) for g in genes]
            lo = int(max(1, min(gend[k] for k in idx if k is not None) + 1 - 999))
            hi = int(min(ref_length, max(gstart[k] for k in idx if k is not None) - 1 + 999))
            for r in blast_map:
                for c in ("sstart", "send"):
                    x = r[c]
                    r[c + ".ref"] = lo + int(x) - 1 if x is not None and 1 <= x <= hi - lo + 1 else None
            win_beg = [w[1] + 333 for w in win_maf]
            win_end = [w[2] + 333 for w in win_maf]
        else:
            if not all(a <= b for a, b in zip(gstart, gend)):
                rcompat.r_stop("Error: all(gbk$start <= gbk$end) is not TRUE")
            if not all(r["sstart"] <= r["send"] for r in blast_map):
                rcompat.r_stop("Error: all(blast.map$sstart <= blast.map$send) is not TRUE")
            k = first.get(gene)
            for r in blast_map:
                r["fm"] = int(r["frame"].split(" ")[0])
            for fm in range(1, 7):
                if fm <= 3:
                    lo = int(max(1, gstart[k] - 999 + (fm - 1)))
                    hi = int(min(ref_length, gend[k] + 999 + (fm - 1)))
                    coords[fm] = coord_frame(ref_length, lo, hi, True)
                else:
                    lo = int(max(1, gstart[k] - 999 - (fm - 4)))
                    hi = int(min(ref_length, gend[k] + 999 - (fm - 4)))
                    coords[fm] = coord_frame(ref_length, lo, hi, False)
            for r in blast_map:
                fm = r["fm"]
                a = coords[fm].get(int(r["sstart"]))
                b = coords[fm].get(int(r["send"]))
                if fm <= 3:
                    r["sstart.ref"], r["send.ref"] = a, None if b is None else b + 2
                else:
                    r["sstart.ref"], r["send.ref"] = None if a is None else a - 2, b
            beg_local = [w[1] + 333 for w in win_maf]
            end_local = [w[2] + 333 for w in win_maf]
            if not all(a < b for a, b in zip(beg_local, end_local)):
                rcompat.r_stop("Error: win.beg.local < win.end.local are not all TRUE")
            win_fm = [r["fm"] for r in blast_map if "correct" in r["frame"]][0]
            cf = coords[win_fm]
            if win_fm <= 3:
                win_beg = [cf.get(int(x)) for x in beg_local]
                win_end = [None if cf.get(int(x)) is None else cf.get(int(x)) + 2 for x in end_local]
            else:
                win_beg = [None if cf.get(int(x)) is None else cf.get(int(x)) - 2 for x in beg_local]
                win_end = [cf.get(int(x)) for x in end_local]

        cols = ["kmer", "Signif", "beta", "MAC", "qstart", "qend", "sstart", "send", "frame", "pident", "length", "mism",
                "gapo", "eval"]
        rows_all = []
        for r in blast_map:
            ev = r["evalue"]
            rows_all.append({"kmer": r["kmer"], "Signif": s3(r["negLog10"]), "beta": s3(r["beta"]), "MAC": r["mac"],
                             "qstart": r["qstart"], "qend": r["qend"], "sstart": r["sstart.ref"], "send": r["send.ref"],
                             "frame": r["frame"], "pident": s3(r["pident"]), "length": r["length"], "mism": r["mismatch"],
                             "gapo": r["gapopen"],
                             "eval": None if ev is None else (math.inf if ev == 0 else float(np.rint(-math.log10(ev))))})

        def le(a, b):
            return None if a is None or b is None else a <= b
        blast_html = []
        for b, e in zip(win_beg, win_end):
            wlo = None if b is None or e is None else min(b, e)
            whi = None if b is None or e is None else max(b, e)
            rows = []
            for r, tr in zip(blast_map, rows_all):
                a, c = r["sstart.ref"], r["send.ref"]
                slo = None if a is None or c is None else min(a, c)
                shi = None if a is None or c is None else max(a, c)
                x, y = le(slo, whi), le(wlo, shi)
                if x is False or y is False:
                    continue
                if x is None or y is None:  # an NA row, whose Signif makes R's if() stop
                    rcompat.r_stop("Error in if (as.numeric(blast.tb[[i]]$Signif[j]) >= thr.signif): missing value "
                                   "where TRUE/FALSE needed")
                rows.append(tr)
            if not rows:
                rcompat.r_stop("Error in if (as.numeric(blast.tb[[i]]$Signif[j]) >= thr.signif): argument is of length zero")
            blast_html.append(ggr.kmer_table_html(cols, rows, thr_signif))

        html_body = NL.join([
            html_body,
            "  <h2>High-resolution Earle plots for " + gene_html + " </h2>",
            "  <p>The series of Figures below are used to identify the underlying",
            "  variants tagged by significant kmers. Resembling Manhattan plots,",
            "  these are high-resolution figures plotting individual kmers",
            "  against the position to which they mapped in the reference genome,",
            "  in the region of " + gene_html + " .",
            "  The kmers are sorted vertically in order of significance, with",
            "  the most significant kmers at the top. The horizontal black",
            "  dashed line demarcates kmers above and below the Bonferroni-corrected",
            "  significance threshold. Only kmers in the correct reading frame are plotted.</p>",
            "",
            "  <p>The kmers are shaded light (<i>&beta;</i>&nbsp;<&nbsp;0)",
            "  or dark (<i>&beta;</i>&nbsp;>&nbsp;0) to indicate direction of association.",
            "  Where there is sequence variation relative to the reference genome,",
            "  individual sites are colour-coded by allele according to the key.",
            "  The reference allele is indicated at the bottom. Only invariant sites",
            "  are coloured grey.",
            "",
            "  Use the arrows to scroll through and jump between windows of",
            "  significance within the region. By default, low-" + t3 + " kmers are filtered",
            "  out. Use the checkbox to remove this filter, which can sometimes assist",
            "  in interpretation of the signal of association. For instance, in the case",
            "  of antimicrobial resistance, there are often multiple very low-" + t3,
            "  mutants associated with increased resistance (darker kmers) which can fall below",
            "  the " + t3 + " threshold. These mutants might have evolved independently, and",
            "  show lower significance than wild types associated with reduced",
            "  resistance (lighter kmers) because their low frequency reduces statistical",
            "  power.</p>",
            "",
            NL])
        html_body = ggr.earle_slideshow(html_body, filenames_kmer_maf, filenames_kmer_maf0, blast_html, t3)
        note = list(ggr.TABLE_NOTE)
        k = note.index('    to least significant (bottom). In the table, <code>beta</code>')
        note[k:k + 1] = ['    to least significant (bottom). All kmers are shown, irrespective of',
                         '    reading frame. In the table, <code>beta</code>']
        html_body = NL.join([html_body] + note + [NL])

    if blast_nul_exists and sum(blast_nul_gd) > 0:
        html_body = NL.join([
            html_body,
            "    <h2>Poorly mapped kmers",
            "    <p>In some cases, kmers were localized to " + gene_html + " by",
            "    the genome aligner (nucmer or bowtie2), but could not be mapped",
            "    with accuracy by BLAST. This discrepancy arises because the former",
            "    used the flanking sequence for context, but the BLAST mapping did",
            "    not. Significant kmers localized to " + gene_html + " but not mapped by",
            "    BLAST are listed in the Table below.</p>",
            NL])
        rows = [{"kmer": r["kmer"], "Signif": s3(r["negLog10"]), "beta": s3(r["beta"]), "MAC": r["mac"]}
                for r, g in zip(blast_nul, blast_nul_gd) if g]
        html_body = NL.join([html_body, ggr.kmer_table_html(["kmer", "Signif", "beta", "MAC"], rows, thr_signif,
                                                              all_bold=True)])

    with open(outfile_html, "w") as fh:
        fh.write(html_head + " " + html_body + " " + html_foot)


if __name__ == "__main__":
    main()
