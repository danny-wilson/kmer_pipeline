#!/usr/bin/env python3
"""gen-gene-report.py: generate the k-mer GWAS HTML report for one of the top
genes (nucleotide k-mers). Port of gen-gene-report.Rscript (Daniel Wilson, 2022).
The HTML is built line by line as the R script builds it."""
import argparse
import math
import os
import re
import sys
import time

import numpy as np

import rcompat
import sequence_functions as sf
from rcompat import r_s3 as s3

NL = "\n"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
gr = __import__("gen-report")


def nth(n):
    if 10 < n % 100 < 20:
        suffix = "th"
    else:
        suffix = (["th", "st", "nd", "rd"] + ["th"] * 96)[n % 10]
    return str(n) + suffix


def rstr(v):
    """paste() of one data frame cell."""
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "NA"
    if isinstance(v, str):
        return v
    return rcompat.r_as_character(v)


def read_delim(path):
    """read.delim(path, stringsAsFactors = FALSE) as a list of row dicts (None for NA)."""
    df = rcompat.r_read_table(path, header=True, sep="\t", quote="\"", comment_char="")
    rows = []
    for rec in df.astype(object).itertuples(index=False, name=None):
        row = {}
        for c, v in zip(df.columns, rec):
            if rcompat._is_na(v):
                v = None
            elif isinstance(v, (int, np.integer)) and not isinstance(v, (bool, np.bool_)):
                v = int(v)
            elif isinstance(v, (float, np.floating)):
                v = float(v)
            row[c] = v
        rows.append(row)
    return list(df.columns), rows


def f_window(s, stem):
    """The window number, start and end parsed from an alignment figure name."""
    res = re.sub(stem, "", s)
    parts = [rcompat.r_as_numeric(p) for p in res.split("_")]
    get = lambda k: parts[k] if k < len(parts) else None  # noqa: E731
    return get(0), get(2), get(4)


def kmer_table_html(columns, rows, thr_signif, signif_col="Signif", all_bold=False):
    lines = ["  <div class='divkmertab'>", "  <table class='kmertab'>", "    <tr>",
             "      <th>" + "</th><th>".join(columns) + "</th>", "    </tr>"]
    for r in rows:
        cells = [sf.kmer_qc_html(r[c]) if c == "kmer" else rstr(r[c]) for c in columns]
        sig = rcompat.r_as_numeric(r[signif_col]) if r[signif_col] is not None else None
        if all_bold or (sig is not None and sig >= thr_signif):
            lines += ["    <tr>", "      <td><b>" + "</b></td><td><b>".join(cells) + "</b></tr>", "    </tr>"]
        else:
            lines += ["    <tr>", "      <td>" + "</td><td>".join(cells) + "</tr>", "    </tr>"]
    legend = sf.kmer_qc_legend([r["kmer"] for r in rows if "kmer" in r])
    lines += ["  </table>", "  </div>"] + ([legend] if legend else []) + [NL]
    return NL.join(lines)


def main():
    rcompat.script_setup(__file__)
    parser = argparse.ArgumentParser(description="gen-gene-report.py Generate a kmer GWAS report for a specific gene. "
                                                 "Daniel Wilson (2022)", allow_abbrev=False)
    for name in ("hit-num", "prefix", "anatype", "k", "refname", "ref-gb", "maf", "alignident", "mincount", "srcdir",
                 "outdir", "logdir"):
        if name == "mincount":  # D4: --plot-min-genomes; --mincount kept as an alias
            parser.add_argument("--plot-min-genomes", "--mincount", dest="mincount", required=True,
                                help="genomes a k-mer/gene combination must be seen in to be plotted (as step 6)")
        else:
            parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    run(args, protein=False)


def run(args, protein):
    import sequence_functions
    HIT_NUM = rcompat.r_as_integer(args.hit_num)
    PREFIX, ANATYPE, K, REFNAME, REF_GB = args.prefix, args.anatype, args.k, args.refname, args.ref_gb
    MAF, ALIGNIDENT, MINCOUNT = args.maf, args.alignident, args.mincount
    SRC, PWD, LOGDIR = args.srcdir, args.outdir, args.logdir
    MACORMAF = "maf" if gr.r_lt(MAF, "1") else "mac"
    FIGDIR = ANATYPE + "kmer" + K + "_kmergenealign_figures/"

    is_maf = True if MAF == "0" else gr.r_lt(MAF, "1")
    is_maf_text3 = "MAF" if is_maf else "MAC"

    stem = PREFIX + "_" + ANATYPE + K + "_" + REFNAME
    filename_topgenes = (FIGDIR + stem + "_top20genes_toppvals_" + MACORMAF + "_" + MAF + "_nucmerAlign_alignIdent_"
                         + ALIGNIDENT + "_alignPosMinCount_" + MINCOUNT + ".txt")
    os.chdir(PWD)
    table_topgenes = rcompat.r_read_table(filename_topgenes)
    if HIT_NUM > len(table_topgenes):
        rcompat.r_cat("Hit", HIT_NUM, "not required: only", len(table_topgenes), "regions in", filename_topgenes, "\n")
        return
    gene = rcompat.r_as_character(table_topgenes.iloc[HIT_NUM - 1, 0])
    gene_html = "<i>" + gene + "</i>"
    gene_filename = gene.replace(":", "_")
    gene_max_signif = rcompat.r_as_numeric_value(table_topgenes.iloc[HIT_NUM - 1, 1])
    is_intergenic = ":" in gene
    genes = [gene]
    genes_html = [gene_html]
    if is_intergenic:
        genes = sequence_functions.r_strsplit(gene, ":")
        genes_html = ["<i>" + g + "</i>" for g in genes]

    html_head = NL.join(["<!DOCTYPE html>", "<html>", "<head>", "  <title>Kmer GWAS report:  " + gene + " </title>",
                         "  <link rel='stylesheet' href='report.css'>", NL])
    html_body = NL.join([
        "</head>", "<body>", "  <h1>Kmer GWAS report: " + gene_html + " </h1>",
        "  <div><p class='timestamp'><code>Prefix: " + PREFIX + "; KmerType: " + ANATYPE + "; K:",
        "  " + K + "; ReferenceGenome: " + REFNAME + "; " + is_maf_text3 + ": " + MAF + "; MinCount:",
        "  " + MINCOUNT + "; AlignIdent: " + ALIGNIDENT + "; ReportTimeStamp:",
        "  " + time.ctime() + ".</code></p></div>", NL])
    html_foot = NL.join(["<script src='report.js'></script>", "</body>", "</html>", ""])

    gbk = sequence_functions.read_dna_seg_from_file(REF_GB, tagsToParse=("CDS",))
    with rcompat.r_open(REF_GB) as f:
        first = f.readline().rstrip("\n")
    print("Read 1 item", file=sys.stderr, flush=True)
    toks = [t for t in first.split(" ") if t != ""]
    ref_length = rcompat.r_as_numeric(toks[2]) if len(toks) >= 3 else None
    if ref_length is None:
        rcompat.r_stop("Error retrieving the reference genome length from the genbank file", "\n")
    ref_length = int(ref_length)
    import reference
    if reference.n_records(REF_GB) > 1:  # D6: the records laid end to end
        ref_length = reference.total_length(REF_GB)

    summary = gr.read_summary_json(PREFIX + "_" + ANATYPE + K + ".summary.json")
    thr_signif = float(summary["bonferroni_threshold"])

    outfile_prefix = PREFIX + "_" + ANATYPE + K + "."
    outfile_html = outfile_prefix + "report_" + gene_filename + ".html"

    gbk_names = list(gbk["name"])
    first = {}
    records = list(gbk["record"]) if "record" in gbk.columns else None
    for k, n in enumerate(gbk_names):
        first.setdefault(n, k)
        if records is not None:  # D6: names used in several records appear as name@record
            first.setdefault(n + "@" + records[k], k)

    def gb(name, col):
        k = first.get(name)
        return None if k is None else gbk[col].iloc[k]

    ggbk = [{c: gb(g, c) for c in ("name", "start", "end", "synonym", "product", "proteinid", "strand")} for g in genes]
    while len(ggbk) < 2:
        ggbk.append({c: None for c in ("name", "start", "end", "synonym", "product", "proteinid", "strand")})

    html_body = NL.join([html_body,
                         "  <p>" + gene_html + " was the " + nth(HIT_NUM) + " most significant",
                         "  region, with a minimum <i>p</i>-value of 10<sup>-" + s3(gene_max_signif) + "</sup>.</p>"])

    def syn(g):  # ifelse(synonym == name, "", " (synonym)"), NA when either is NA
        if g["synonym"] is None or g["name"] is None:
            return "NA"
        return "" if g["synonym"] == g["name"] else " (" + g["synonym"] + ")"

    def length(g):
        if g["start"] is None or g["end"] is None:
            return "NA"
        return str(int(abs(g["end"] - g["start"]) + 1))
    gh = genes_html + ["NA"] * (2 - len(genes_html))
    if is_intergenic:
        html_body = NL.join([
            html_body,
            "  <p>This is an intergenic region.",
            "  The user-provided Genbank file lists " + gh[0] + syn(ggbk[0]),
            "  as " + length(ggbk[0]) + " nucleotides long and " + gh[1] + syn(ggbk[1]),
            "  as " + length(ggbk[1]) + " nucleotides long.",
            "  They encode the " + rstr(ggbk[0]["product"]) + " (protein ID " + rstr(ggbk[0]["proteinid"]) + ")",
            "  and the " + rstr(ggbk[1]["product"]) + " (protein ID " + rstr(ggbk[1]["proteinid"]) + ").</p>", NL])
    else:
        html_body = NL.join([
            html_body,
            "  <p>The user-provided Genbank file lists " + gh[0] + syn(ggbk[0]),
            "  as " + length(ggbk[0]) + " nucleotides long.",
            "  It encodes the " + rstr(ggbk[0]["product"]) + " (protein ID " + rstr(ggbk[0]["proteinid"]) + ").</p>", NL])

    if protein:
        return outfile_html, html_head, html_body, html_foot, locals()

    filename_Manhattan_maf = FIGDIR + stem + "_" + gene + "_Manhattan_" + MACORMAF + MAF + ".png"
    filename_Manhattan_maf0 = FIGDIR + stem + "_" + gene + "_Manhattan_allkmers.png"
    filename_kmer_maf = FIGDIR + stem + "_" + gene + "_plot_*_pos_*_to_*_" + MACORMAF + MAF + "_alignment.png"
    filename_kmer_maf0 = FIGDIR + stem + "_" + gene + "_plot_*_pos_*_to_*_alignment.png"
    filename_mapped = FIGDIR + stem + "_top_gene_*_" + gene + "_blast_results.txt"
    filename_unmapped = FIGDIR + stem + "_top_gene_*_" + gene + "_no_blast_result_or_poor_alignment.txt"

    html_body = NL.join([
        html_body,
        "  <h2>Manhattan plot for " + gene_html + " </h2>",
        "  <p>The Figure displays the significance of each kmer against the position in the",
        "  reference genome to which it mapped, with a focus on " + gene_html + " .",
        "  The Bonferroni-corrected significance threshold is shown as a horizontal black dashed line.",
        "  Annotated features are plotted below. Points are shaded light (<i>&beta;</i>&nbsp;<&nbsp;0)",
        "  or dark (<i>&beta;</i>&nbsp;>&nbsp;0) to indicate direction of association, and colour-coded",
        "  grey (unique) or orange (non-unique) to indicate the quality of mapping.",
        "  When <i>&beta;</i>&nbsp;>&nbsp;0, the presence of the kmer is associated with larger values of the phenotype.",
        "  The figure can be displayed with or without filtering of kmers below",
        "  the " + is_maf_text3 + " threshold (although the significance",
        "  threshold is not updated since we do not recommend reporting low-" + is_maf_text3 + " kmers",
        "  as significant).</p>", NL])
    # The filtered plot is not drawn when no kmer passes the threshold
    plots = [(f, c) for f, c in ((filename_Manhattan_maf, "Kmers mapping to the region, filtered by " + is_maf_text3 + "."),
                                 (filename_Manhattan_maf0, "Kmers mapping to the region. No " + is_maf_text3 + " filter."))
             if os.path.exists(f)]
    if plots:
        html_body = gr.slideshow(html_body, [f for f, _ in plots], [c for _, c in plots], "")
    html_body = NL.join([html_body, NL])

    def rdir(pattern):
        d = os.path.dirname(pattern)
        return rcompat.r_dir(d, glob=os.path.basename(pattern), full_names=True) if os.path.isdir(d) else []
    filenames_mapped = rdir(filename_mapped)
    filenames_unmapped = rdir(filename_unmapped)
    blast_nul_exists = len(filenames_unmapped) > 0

    if not filenames_mapped:
        rcompat.r_stop("Error in file(file, \"rt\"): invalid 'description' argument (no BLAST results for ", gene, ")")
    blast_map = []
    for fn in filenames_mapped:
        blast_map += read_delim(fn)[1]
    if blast_nul_exists:
        blast_nul = []
        for fn in filenames_unmapped:
            blast_nul += read_delim(fn)[1]
        mapped_kmers = {r["kmer"] for r in blast_map}
        blast_nul_gd = [r["kmer"] not in mapped_kmers and r["negLog10"] is not None and r["negLog10"] >= thr_signif
                        for r in blast_nul]

    filenames_kmer_maf = rdir(filename_kmer_maf)
    maf_set = set(filenames_kmer_maf)
    filenames_kmer_maf0 = [f for f in rdir(filename_kmer_maf0) if f not in maf_set]

    if filenames_kmer_maf:
        pstem = FIGDIR + stem + "_" + gene + "_plot_"  # used as a regular expression, as in R
        win_maf = [f_window(f, pstem) for f in filenames_kmer_maf]
        win_maf0 = [f_window(f, pstem) for f in filenames_kmer_maf0]
        if win_maf != win_maf0:
            rcompat.r_stop("Error: all(unname(win.maf) == unname(win.maf0)) is not TRUE")
        win_beg = [w[1] for w in win_maf]
        win_end = [w[2] for w in win_maf]

        gb_names = gbk_names
        gstart = [float(v) for v in gbk["start"]]
        gend = [float(v) for v in gbk["end"]]
        gstrand = [float(v) for v in gbk["strand"]]
        # The window stays within the reference (the region's record, D6)
        rec_lo, rec_hi = 1, ref_length
        if records is not None:
            k0 = first.get(genes[0] if is_intergenic else gene)
            if k0 is not None:
                rec = reference.record_of(reference.records(REF_GB), gstart[k0])
                rec_lo, rec_hi = rec.start, rec.end
        if is_intergenic:
            idx = [first.get(g) for g in genes]
            lo = max(rec_lo, min(gend[k] for k in idx if k is not None) + 1 - 999)
            hi = min(rec_hi, max(gstart[k] for k in idx if k is not None) - 1 + 999)
            forward = True
        else:
            k = first.get(gene)
            lo = max(rec_lo, gstart[k] - 999)
            hi = min(rec_hi, gend[k] + 999)
            forward = gstrand[k] == 1
        lo, hi = int(lo), int(hi)

        def to_ref(x):  # match(x, gene.coords.in.ref)
            if x is None:
                return None
            x = int(x)
            if not 1 <= x <= hi - lo + 1:
                return None
            return lo + x - 1 if forward else hi - x + 1
        for r in blast_map:
            r["sstart.ref"] = to_ref(r["sstart"])
            r["send.ref"] = to_ref(r["send"])

        cols = ["kmer", "Signif", "beta", "MAC", "qstart", "qend", "sstart", "send", "pident", "length", "mism", "gapo",
                "eval"]
        table_rows = []
        for r in blast_map:
            ev = r["evalue"]
            table_rows.append({"kmer": r["kmer"], "Signif": s3(r["negLog10"]), "beta": s3(r["beta"]), "MAC": r["mac"],
                               "qstart": r["qstart"], "qend": r["qend"], "sstart": r["sstart.ref"], "send": r["send.ref"],
                               "pident": s3(r["pident"]), "length": r["length"], "mism": r["mismatch"],
                               "gapo": r["gapopen"],
                               "eval": None if ev is None else (math.inf if ev == 0 else float(np.rint(-math.log10(ev))))})
        blast_html = []
        for b, e in zip(win_beg, win_end):
            wlo, whi = min(b, e), max(b, e)
            rows = []
            for r, tr in zip(blast_map, table_rows):
                a, c = r["sstart.ref"], r["send.ref"]
                if a is None or c is None:
                    continue  # NA in a logical row index gives an NA row in R; not reached when mapped
                if min(a, c) <= whi and max(a, c) >= wlo:
                    rows.append(tr)
            if not rows:
                rcompat.r_stop("Error in if (as.numeric(blast.tb[[i]]$Signif[j]) >= thr.signif): argument is of length zero")
            blast_html.append(kmer_table_html(cols, rows, thr_signif))

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
            "  significance threshold.</p>",
            "",
            "  <p>The kmers are shaded light (<i>&beta;</i>&nbsp;<&nbsp;0)",
            "  or dark (<i>&beta;</i>&nbsp;>&nbsp;0) to indicate direction of association.",
            "  Where there is sequence variation relative to the reference genome,",
            "  individual sites are colour-coded by allele according to the key.",
            "  The reference allele is indicated at the bottom. Only invariant sites",
            "  are coloured grey.",
            "",
            "  Use the arrows to scroll through and jump between windows of",
            "  significance within the region. By default, low-" + is_maf_text3 + " kmers are filtered",
            "  out. Use the checkbox to remove this filter, which can sometimes assist",
            "  in interpretation of the signal of association. For instance, in the case",
            "  of antimicrobial resistance, there are often multiple very low-" + is_maf_text3,
            "  mutants associated with increased resistance (darker kmers) which can fall below",
            "  the " + is_maf_text3 + " threshold. These mutants might have evolved independently, and",
            "  show lower significance than wild types associated with reduced",
            "  resistance (lighter kmers) because their low frequency reduces statistical",
            "  power.</p>",
            "",
            NL])
        html_body = earle_slideshow(html_body, filenames_kmer_maf, filenames_kmer_maf0, blast_html, is_maf_text3)
        html_body = NL.join([html_body] + TABLE_NOTE + [NL])

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
        html_body = NL.join([html_body, kmer_table_html(["kmer", "Signif", "beta", "MAC"], rows, thr_signif,
                                                         all_bold=True)])

    with open(outfile_html, "w") as f:
        f.write(html_head + " " + html_body + " " + html_foot)


TABLE_NOTE = [
    '    <p>The Table above provides detailed information on the',
    '    kmers plotted in the Figure, ordered from most significant (top)',
    '    to least significant (bottom). In the table, <code>beta</code>',
    '    provides the direction and magnitude of the association between the',
    '    phenotype and the presence of the kmer, and <code>MAC</code> provides the',
    '    minor allele count (no filter was applied to the Table).',
    '    The remaining columns were produced by BLAST: <code>qstart</code>, ',
    '    <code>qend</code>, <code>sstart</code> and <code>send</code> provide',
    '    the start and end coordinates of the BLAST match for the query (kmer)',
    '    and subject (reference genome). The match is further summarized by',
    '    the <code>pident</code> (percent identity), <code>length</code>,',
    '    number of <code>mism[atches]</code>, <code>gapo[pen]</code> events,',
    '    and the log<sub>10</sub> of the <code>eval[ue]</code>.</p>']


def earle_slideshow(html_body, files_maf, files_maf0, blast_html, is_maf_text3):
    n = len(files_maf)
    lines = [html_body, '  <div class="slideshow-container">',
             '      <div class="toggler"><label>Filter by ' + is_maf_text3 + ' <input type="checkbox" id="toggler" '
             'value="yes" onclick="showToggled()" checked></label></div>']
    for i in range(n):
        lines += ['    <div class="mySlides fade">',
                  '      <div class="numbertext">' + str(i + 1) + ' / ' + str(n) + '</div>',
                  '      <img src="' + files_maf[i] + '" class="center toggled" style="width:80%">',
                  '      <img src="' + (files_maf0[i] if i < len(files_maf0) else "NA") + '" class="center untoggled" '
                  'style="width:80%">',
                  blast_html[i],
                  '      <div class="text toggled"></div>',
                  '      <div class="text untoggled"></div>',
                  '    </div>']
    lines += ['    <a class="prev" onclick="plusSlides(-1,this)">&#10094;</a>',
              '    <a class="next" onclick="plusSlides(1,this)">&#10095;</a>',
              '    <br>',
              '    <div style="text-align:center">']
    lines += ['    <span class="dot" onclick="currentSlide(' + str(i) + ',this)"></span>' for i in range(1, n + 1)]
    lines += ['    </div>', '  </div>', NL]
    return NL.join(lines)


if __name__ == "__main__":
    main()
