"""Functions for the Manhattan and QQ plots and the top-gene tables.
Port of Manhattan_functions.R.

R keeps the GEMMA results as a character matrix, so every number read from it
(beta, -log10 p) is the value R printed with 15 significant digits; the port
does the same (as_r_text_number). Figures are drawn with matplotlib: they show
the same data as R's but are not pixel copies (PLAN 5.4)."""
import math
import os
import sys

import numpy as np

import rcompat
from rcompat import r_cat, r_paste0, r_stop

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

### Get colours
colour_selection = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00"]

CM = 1 / 2.54


def as_r_text_number(v):
    """A double as R stores it in a character matrix and reads it back:
    as.numeric(as.character(v)), i.e. rounded to 15 significant digits."""
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return math.nan
    return float(rcompat.r_as_character(float(v)))


def create_figures_dir(dir, kmer_type, kmer_length, alignmenttype):
    figures_dir = dir + "/" + r_paste0(kmer_type, "kmer", kmer_length, "_", alignmenttype, "_figures/")  # file.path
    if not os.path.isdir(figures_dir):
        rcompat.r_dir_create(figures_dir)
    return figures_dir


def match_first(values, table):
    """match(values, table): 1-based index of the first match, or None."""
    first = {}
    for k, t in enumerate(table):
        first.setdefault(t, k + 1)
    return [first.get(v) for v in values]


def get_genes_to_plot(gene_names, y, gene_conversion, ymax, gene_panel, ref, xadjust=None, ngenes=20):
    """Rows (genes, ytop, xadjust, replace_gene_name, gene_col) for the genes to
    label: the ngenes genes with the highest y, in reference order."""
    if gene_names is None:
        return None
    if y is None:
        top_genes = list(gene_names)
    else:
        o = rcompat.r_order([math.nan if v is None else v for v in y], decreasing=True)
        top_genes = rcompat.r_unique([gene_names[k] for k in o])[:ngenes]
    ref_names = list(ref["name"])
    m = match_first(top_genes, ref_names)
    o2 = rcompat.r_order([math.nan if v is None else float(v) for v in m])
    if xadjust is not None:
        xadjust = [xadjust[k] for k in o2]
    top_genes = [top_genes[k] for k in o2]
    gene_name_conversion = [gene_conversion.get(g) if gene_conversion is not None else None for g in top_genes]
    gene_name_conversion = [t if c is None else c for c, t in zip(gene_name_conversion, top_genes)]
    gene_col = ["black"] * len(top_genes)
    panel = set(gene_panel or [])
    for k, g in enumerate(gene_name_conversion):
        if g in panel:
            gene_col[k] = "#D55E00"
    wh = [sum(1 for part in g.split(":") if part in panel) for g in gene_name_conversion]
    if any(w > 0 for w in wh):
        for k in range(len(wh)):
            if wh[k] == 1 and gene_col[k] != "#D55E00":
                gene_col[k] = "#E69F00"
            if wh[k] == 2:
                gene_col[k] = "#D55E00"
    if xadjust is None:
        xadjust = [0.0] * len(top_genes)
    ytops = [ymax[0] + ymax[1] / 40, ymax[0] + ymax[1] / 12]
    return [(g, ytops[k % 2], xadjust[k], c, col)
            for k, (g, c, col) in enumerate(zip(top_genes, gene_name_conversion, gene_col))]


def gene_span(gene, ref):
    names = list(ref["name"])
    starts, ends = list(ref["start"]), list(ref["end"])
    if ":" not in gene:
        idx = [k for k, n in enumerate(names) if n == gene]
        return float(starts[idx[0]]), float(ends[idx[-1]])
    parts = gene.split(":")
    i1 = [k for k, n in enumerate(names) if n == parts[0]]
    i2 = [k for k, n in enumerate(names) if len(parts) > 1 and n == parts[1]]
    pos1 = float(ends[i1[-1]]) + 1 if i1 else math.nan
    pos2 = float(starts[i2[0]]) - 1 if i2 else math.nan
    return pos1, pos2


def plot_gene_lines(ax, rows, ref, col="#cecece", gene_name_cex=0.6):
    """Dotted lines from 0 to each label height, at the middle of each gene, with
    the gene name at 45 degrees (as plot_gene_lines(rect = FALSE))."""
    for genes, ytop, xadj, replace_name, name_col in rows:
        if genes is None:
            continue
        pos1, pos2 = gene_span(genes, ref)
        mid = pos1 + (pos2 - pos1) / 2
        ax.plot([mid, mid], [0, ytop], linestyle=":", color=col, linewidth=0.6, clip_on=False, zorder=1)
        label = replace_name if replace_name else genes
        ax.text(mid + xadj, ytop, label, rotation=45, ha="left", va="bottom", fontsize=7 * gene_name_cex,
                color=name_col, clip_on=False)


def extract_lambda_lognull(datafiles):
    a = rcompat.r_scan_lines(datafiles, quiet=True)
    b = a[12].split(" ")[-1]
    c = a[16].split(" ")[-1]
    return {"lambda": b, "lognull": c}


def get_loglik(LH1, lognull):
    return 2 * (LH1 - rcompat.r_as_numeric(lognull))


def read_gemma_files(input_dir, prefix, kmer_type, kmer_length, nPatterns):
    """The GEMMA results for every pattern, in pattern order: a list of rows
    [rs, beta, se, p_lrt, logl_H1, negLog10] of strings (None for untested
    patterns), as R's character matrix."""
    stem = r_paste0(input_dir, prefix, "_", kmer_type, kmer_length)
    files = rcompat.r_system_intern("ls " + stem + "*-*.assoc.txt.gz")
    file_range = [f.replace(stem + ".", "").replace(".assoc.txt.gz", "") for f in files]
    file_beg = [int(s.split("-")[0]) for s in file_range]
    file_end = [int(s.split("-")[1]) for s in file_range]
    if max(file_end) != nPatterns:
        r_stop("Error: max gemma pattern index does not equal total number of patterns", "\n")
    covered = set()
    for b, e in zip(file_beg, file_end):
        covered.update(rcompat.r_colon(b, e))
    if any(k not in covered for k in range(1, int(nPatterns) + 1)):
        r_stop("Error: not all patterns are present in gemma files", "\n")
    files = [files[k] for k in rcompat.r_order(file_beg)]

    assoc = []
    header = None
    for f in files:
        lines = rcompat.r_pipe("zcat " + f + " | cut -f2,5,6,10,12").split("\n")
        if lines and lines[-1] == "":
            lines.pop()
        if header is None:
            header = lines[0].split("\t")
        for line in lines[1:]:
            fields = line.split("\t")
            if len(fields) != 5:
                r_stop("Error in matrix(gemma.i, ncol = 5): GEMMA line with ", len(fields), " fields in ", f)
            assoc.append(fields)

    gemma_log_file = rcompat.r_system_intern("ls " + stem + ".1-*.log.txt.gz")
    l0 = extract_lambda_lognull(gemma_log_file[0])["lognull"]

    D = np.array([get_loglik(rcompat.r_as_numeric(r[4]), l0) for r in assoc], dtype=float)
    pvals = rcompat.neg_log10_pchisq1(D)
    for r, p in zip(assoc, pvals):
        r.append(rcompat.r_as_character(float(p)))
    tested = set()
    for r in assoc:
        tested.add(rcompat.r_as_numeric(r[0]))
    r_cat("Number of untested patterns:", sum(1 for k in range(1, int(nPatterns) + 1) if float(k) not in tested), "\n")
    pv = [rcompat.r_as_numeric(r[3]) for r in assoc]
    nl = [rcompat.r_as_numeric(r[5]) for r in assoc]
    r_cat("GEMMA range of pvalues:", min(pv), max(pv), "\n")
    r_cat("GEMMA range of -log10(pvalues):", min(nl), max(nl), "\n")

    # Match gemma results to patterns
    by_rs = {}
    for r in assoc:
        by_rs.setdefault(rcompat.r_as_numeric(r[0]), r)
    out = [by_rs.get(float(k)) for k in range(1, int(nPatterns) + 1)]
    r_cat("Matched gemma results to all patterns", "\n")
    return out


def assoc_column(assoc, j):
    """as.numeric(assoc[, j]) (1-based column) with NaN for untested patterns."""
    return np.array([math.nan if r is None else (rcompat.r_as_numeric(r[j - 1]) if rcompat.r_as_numeric(r[j - 1]) is not None
                                                  else math.nan) for r in assoc], dtype=float)


def plot_QQ(kmerIndex, assoc, output_dir, prefix, minor_allele_threshold, macormaf, mapatterns, kmer_type, kmer_length):
    uk = rcompat.r_unique(list(kmerIndex))
    ma_u = np.array([mapatterns[k - 1] for k in uk], dtype=float)
    with np.errstate(invalid="ignore"):
        which_kmers = np.flatnonzero(ma_u > 0) if minor_allele_threshold == 0 else \
            np.flatnonzero(ma_u >= minor_allele_threshold)
    n = len(which_kmers)
    qq_x = -np.log10(np.arange(1, n + 1) / n) if n else np.array([])
    neg = assoc_column(assoc, 6)
    qq_y = np.array([neg[uk[k] - 1] for k in which_kmers], dtype=float)
    qq_y = qq_y[rcompat.r_order(qq_y, decreasing=True)] if n else qq_y

    if minor_allele_threshold == 0:
        file_suffix = "_QQplot_allkmers.png"
    else:
        file_suffix = r_paste0("_QQplot_", macormaf, minor_allele_threshold, ".png")
    fig, ax = plt.subplots(figsize=(12 * CM, 12 * CM))
    ax.plot(qq_x, qq_y, color="black", linewidth=0.8)
    lim = [0, max(np.nanmax(qq_x) if n else 1, np.nanmax(qq_y) if n and np.isfinite(qq_y).any() else 1)]
    ax.plot(lim, lim, color="red", linestyle="--", linewidth=0.8)
    ax.set_xlabel(r"Null distribution of -log$_{10}$ $\it{p}$ values", fontsize=8)
    ax.set_ylabel(r"Empirical distribution of -log$_{10}$ $\it{p}$ values", fontsize=8)
    ax.tick_params(labelsize=7)
    fig.tight_layout()
    fig.savefig(output_dir + prefix + "_" + kmer_type + rcompat.r_as_character(kmer_length) + file_suffix, dpi=600)
    plt.close(fig)


def _ramp(c1, c2, t):
    """colorRamp(c(c1, c2))(t) then rgb(maxColorValue = 256)."""
    a = np.array(matplotlib.colors.to_rgb(c1)) * 255
    b = np.array(matplotlib.colors.to_rgb(c2)) * 255
    out = []
    for v in t:
        rgb = a + (b - a) * v
        out.append("#%02X%02X%02X" % tuple(int(round(x / 256 * 255)) for x in rgb))
    return out


def get_Manhattan_colours(final_kmer_pos_index, assoc_patterns, kmerIndex, colour_selection, ypos, bonferroni,
                          mafpatterns, pheno_type):
    n = len(final_kmer_pos_index)
    ## Colour by whether the kmer has mapped more than once
    counts = {}
    for v in final_kmer_pos_index:
        counts[v] = counts.get(v, 0) + 1
    multialignCOL = [colour_selection[5] if counts[v] > 1 else "grey50" for v in final_kmer_pos_index]
    r_cat("Created multialignCOL", "\n")

    beta_all = assoc_column(assoc_patterns, 2)
    r_cat("Range beta:", np.nanmin(beta_all), np.nanmax(beta_all), "\n")
    beta = np.array([beta_all[kmerIndex[int(v) - 1] - 1] for v in final_kmer_pos_index], dtype=float)
    betaCOL = ["grey50"] * len(ypos)
    if pheno_type == "binary":
        for k, b in enumerate(beta):
            if b > 0:
                betaCOL[k] = colour_selection[5]
            elif b < 0:
                betaCOL[k] = colour_selection[4]
    else:
        pos = np.flatnonzero(beta > 0)
        neg = np.flatnonzero(beta < 0)
        if len(pos):
            bp = beta[pos]
            bp = (bp - np.nanmin(bp)) / (np.nanmax(bp) - np.nanmin(bp)) if np.nanmax(bp) > np.nanmin(bp) else bp * np.nan
            for k, c in zip(pos, _ramp("grey50" if False else "#7F7F7F", colour_selection[5], np.nan_to_num(bp))):
                betaCOL[k] = c
        if len(neg):
            bn = beta[neg]
            bn = (bn - np.nanmin(bn)) / (np.nanmax(bn) - np.nanmin(bn)) if np.nanmax(bn) > np.nanmin(bn) else bn * np.nan
            for k, c in zip(neg, _ramp(colour_selection[4], "#7F7F7F", np.nan_to_num(bn))):
                betaCOL[k] = c
    with np.errstate(invalid="ignore"):
        for k in np.flatnonzero(np.asarray(ypos, dtype=float) < bonferroni):
            betaCOL[k] = "grey50"
    r_cat("Created betaCOL", "\n")

    maf = np.array([mafpatterns[kmerIndex[int(v) - 1] - 1] for v in final_kmer_pos_index], dtype=float)
    mafCOL = ["grey50"] * n
    with np.errstate(invalid="ignore"):
        for k, m in enumerate(maf):
            if m < 0.01:
                mafCOL[k] = colour_selection[5]
            elif 0.01 <= m < 0.05:
                mafCOL[k] = colour_selection[4]
            elif m >= 0.05:
                mafCOL[k] = colour_selection[2]
    r_cat("Created mafCOL", "\n")
    return {"multialignCOL": multialignCOL, "betaCOL": betaCOL, "mafCOL": mafCOL}


def write_top_gene_kmers_to_file(wh_i, final_kmer_list, final_kmer_pos_index, assoc, kmerIndex, mac, output_file):
    beta_all = assoc_column(assoc, 2)
    neg_all = assoc_column(assoc, 6)
    rows = []
    for w in wh_i:
        k = int(final_kmer_pos_index[w])
        pat = kmerIndex[k - 1]
        rows.append((final_kmer_list[k - 1], neg_all[pat - 1], beta_all[pat - 1], mac[k - 1]))
    # Remove kmers which have not been tested for this phenotype
    rows = [r for r in rows if not math.isnan(r[1])]
    # Order from most significant to least
    o = rcompat.r_order([r[1] for r in rows], decreasing=True)
    rows = [rows[k] for k in o]
    with open(output_file, "w") as f:
        f.write("kmer\tnegLog10\tbeta\tmac\n")
        for r in rows:
            f.write("\t".join(rcompat.r_str(v, 15) if not isinstance(v, str) else v for v in r) + "\n")


def top20genes(gene_names, ma, minor_allele_threshold, ypos, macormaf, output_dir, prefix, min_count, ident_threshold,
               kmer_type, kmer_length, ref_name):
    which = [k for k in range(len(gene_names)) if gene_names[k] is not None and ma[k] >= minor_allele_threshold]
    o = rcompat.r_order([ypos[k] for k in which], decreasing=True)
    top = rcompat.r_unique([gene_names[which[k]] for k in o])[:20]
    pvals = []
    for g in top:
        vals = [ypos[k] for k in which if gene_names[k] == g and not math.isnan(ypos[k])]
        pvals.append(max(vals) if vals else -math.inf)
    r_cat(r_paste0("Top 20 genes above the ", macormaf, " threshold ", minor_allele_threshold, ":"), "\n")
    r_cat("Gene -log10pvalue", "\n")
    for g, p in zip(top, pvals):  # apply() of a character matrix: the p-value as text
        r_cat(g, rcompat.r_as_character(p), "\n")
    outprefix = r_paste0(output_dir, prefix, "_", kmer_type, kmer_length, "_", ref_name, "_top20genes_toppvals_", macormaf,
                         "_", minor_allele_threshold)
    if min_count is not None:
        outfile = r_paste0(outprefix, "_nucmerAlign_alignIdent_", ident_threshold, "_alignPosMinCount_", min_count, ".txt")
    else:
        outfile = outprefix + "_bowtie2mapping.txt"
    with open(outfile, "w") as f:
        f.write("".join(g + "\t" + rcompat.r_as_character(p) + "\n" for g, p in zip(top, pvals)))


R_GREY50 = "#7F7F7F"


def _col(c):
    return R_GREY50 if c == "grey50" else c


def plot_manhattan(outfilename, xpos, ma_threshold_pass, ypos, ylims_i, annotateGeneFile, ref, which_genes_to_annotate_i,
                   allCOLS, allPCH, i, bonferroni, legendtext, legendcol, legendpch, legendlty, beta, gene_names,
                   gene_conversion, pheno_type, ref_length, filecol):
    fig = plt.figure(figsize=(22 * CM, 12 * CM))
    ax = fig.add_axes([0.08, 0.13, 0.62, 0.7])
    x = np.asarray(xpos, dtype=float)[ma_threshold_pass]
    y = np.asarray(ypos, dtype=float)[ma_threshold_pass]
    cols = [_col(allCOLS[i][k]) for k in ma_threshold_pass]
    if ylims_i is not None:
        ax.set_ylim(ylims_i)
    else:
        finite = y[np.isfinite(y)]
        top = finite.max() if len(finite) else 1
        ax.set_ylim(-0.04 * top, top * 1.04)
    if len(x):
        ax.set_xlim(np.nanmin(x) - 0.04 * (np.nanmax(x) - np.nanmin(x)), np.nanmax(x) + 0.04 * (np.nanmax(x) - np.nanmin(x)))
    lo, hi = ax.get_ylim()
    ymax = (hi, hi - lo)
    if annotateGeneFile is not None:
        annotateGene = rcompat.r_scan_lines(annotateGeneFile, quiet=True)
        r_cat("Genes/IRs to annotate on the Manhattan plot:", " ".join(annotateGene), "\n")
        rows = get_genes_to_plot(annotateGene, None, {g: g for g in annotateGene}, ymax, [], ref,
                                 xadjust=[0.0] * len(annotateGene))
    else:
        rows = get_genes_to_plot([gene_names[k] for k in which_genes_to_annotate_i],
                                 [ypos[k] for k in which_genes_to_annotate_i],
                                 gene_conversion, ymax, [], ref)
    plot_gene_lines(ax, rows, ref)
    ax.scatter(x, y, s=4, facecolors="none", edgecolors=cols, linewidths=0.4, zorder=2)
    ax.set_xlabel("Position in reference genome (Mb)", fontsize=8)
    ax.set_ylabel(r"Significance (-log$_{10}$ $\it{p}$) LMM", fontsize=8)
    ticks = [k * 1e6 for k in range(int(math.floor(ref_length / 1e6)) + 1)]
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(k) for k in range(len(ticks))], fontsize=7)
    ax.tick_params(axis="y", labelsize=7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.axhline(bonferroni, color="black", linestyle="--", linewidth=0.8)
    handles, labels = [], []
    for t, c, p, l in zip(legendtext[i], legendcol[i], legendpch, legendlty):
        if t == "" or (filecol[i] == "betaCOL" and pheno_type == "continuous" and t in legendtext[i][2:5]):
            continue
        if l == 2:
            handles.append(plt.Line2D([], [], color=c, linestyle="--"))
        else:
            handles.append(plt.Line2D([], [], color=_col(c), marker="o" if p == 16 else None, linestyle=""))
        labels.append(t)
    fig.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.72, 0.97), fontsize=6, frameon=True)
    fig.savefig(outfilename, dpi=600)
    plt.close(fig)


def write_summary_json(summary_file, n_kmers, n_patterns, n_untested_patterns, max_neglog10p, minor_allele_threshold,
                       macormaf, n_tests, bonferroni):
    lines = ["{",
             '  "n_kmers": %d,' % int(n_kmers),
             '  "n_patterns": %d,' % int(n_patterns),
             '  "n_untested_patterns": %d,' % int(n_untested_patterns),
             '  "max_neglog10p": %.17g,' % max_neglog10p,
             '  "minor_allele_threshold": %.17g,' % minor_allele_threshold,
             '  "macormaf": "%s",' % macormaf,
             '  "n_tests": %d,' % int(n_tests),
             '  "bonferroni_threshold": %.17g' % bonferroni,
             "}"]
    rcompat.r_cat_lines(lines, summary_file)
    r_cat("Written summary for reports:", summary_file, "\n")


def get_pheno_type(pheno):
    values = {p for p in pheno if p is not None and not (isinstance(p, float) and math.isnan(p))}
    return "binary" if len(values) == 2 else "continuous"
