"""Functions for the Manhattan and QQ plots and the top-gene tables.
Port of Manhattan_functions.R.

R keeps the GEMMA results as a character matrix, so every number read from it
(beta, -log10 p) is the value R printed with 15 significant digits; the port
does the same (as_r_text_number). Figures are drawn by plot_figures.R from the
figure-data files written here (FigureData; PLAN 5.5)."""
import math
import os
import sys

import numpy as np

import rcompat
from rcompat import r_cat, r_paste0, r_stop

### Get colours
colour_selection = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00"]


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


def extract_lambda_lognull(datafiles):
    a = rcompat.r_scan_lines(datafiles, quiet=True)
    b = a[12].split(" ")[-1]
    c = a[16].split(" ")[-1]
    return {"lambda": b, "lognull": c}


def get_loglik(LH1, lognull):
    return 2 * (LH1 - rcompat.r_as_numeric(lognull))


def r_range_text(values):
    """min and max as R's cat(min(x), max(x)) prints them: NA NA if any value is NA (None)."""
    if not values or any(v is None or v != v for v in values):
        return ["NA", "NA"]
    return [min(values), max(values)]


def read_gemma_files(input_dir, prefix, kmer_type, kmer_length, nPatterns):
    """The GEMMA results for every pattern, in pattern order: a list of rows
    [rs, beta, se, p_lrt, logl_H1, negLog10] of strings (None for untested
    patterns), as R's character matrix."""
    stem = r_paste0(input_dir, prefix, "_", kmer_type, kmer_length)
    files = rcompat.r_system_intern("ls " + stem + "*-*.assoc.txt.gz")
    file_range = [f.replace(stem + ".", "").replace(".assoc.txt.gz", "") for f in files]
    file_beg = [rcompat.parse_index(s.split("-")[0]) for s in file_range]  # also "1e+05" (D1a)
    file_end = [rcompat.parse_index(s.split("-")[1]) for s in file_range]
    if max(file_end) != nPatterns:
        r_stop("Error: max gemma pattern index does not equal total number of patterns", "\n")
    # N1: the batches must cover patterns 1..n exactly once (results of another run with a
    # different number of tasks would overlap them)
    expected = 1
    for b, e in sorted(zip(file_beg, file_end)):
        if b != expected:
            r_stop("Error: the GEMMA result files in ", input_dir, " do not cover the patterns exactly once (",
                   ", ".join(str(x) + "-" + str(y) for x, y in sorted(zip(file_beg, file_end))),
                   "): remove the results of earlier runs (rerun step 4 with overwrite = true)", "\n")
        expected = e + 1
    if expected != int(nPatterns) + 1:
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
    pvals = rcompat.neg_log10_pchisq1_r(D)
    for r, p in zip(assoc, pvals):
        r.append(rcompat.r_as_character(float(p)))
    tested = set()
    for r in assoc:
        tested.add(rcompat.r_as_numeric(r[0]))
    r_cat("Number of untested patterns:", sum(1 for k in range(1, int(nPatterns) + 1) if float(k) not in tested), "\n")
    pv = [rcompat.r_as_numeric(r[3]) for r in assoc]
    nl = [rcompat.r_as_numeric(r[5]) for r in assoc]
    # As R's min() and max(): NA if any value is NA (N10: GEMMA writes -nan for a pattern it
    # cannot fit, e.g. one collinear with the covariates among the analysed genomes)
    r_cat("GEMMA range of pvalues:", r_range_text(pv), "\n")
    r_cat("GEMMA range of -log10(pvalues):", r_range_text(nl), "\n")

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


def assoc_value(assoc, pat, j):
    """as.numeric(assoc[, j])[pat] keeping R's distinction: None for NA (untested
    pattern), nan for NaN (GEMMA's "-nan")."""
    row = assoc[pat - 1]
    return None if row is None else rcompat.r_as_numeric(row[j - 1])


def write_top_gene_kmers_to_file(wh_i, final_kmer_list, final_kmer_pos_index, assoc, kmerIndex, mac, output_file):
    rows = []
    for w in wh_i:
        k = int(final_kmer_pos_index[w])
        pat = kmerIndex[k - 1]
        rows.append((final_kmer_list[k - 1], assoc_value(assoc, pat, 6), assoc_value(assoc, pat, 2), mac[k - 1]))
    # Remove kmers which have not been tested for this phenotype (is.na: NA or NaN)
    rows = [r for r in rows if r[1] is not None and not math.isnan(r[1])]
    # Order from most significant to least
    o = rcompat.r_order([r[1] for r in rows], decreasing=True)
    rows = [rows[k] for k in o]

    def text(v):  # as.character() in cbind(): NaN stays "NaN"
        if isinstance(v, str):
            return v
        if isinstance(v, float) and math.isnan(v):
            return "NaN"
        return rcompat.r_str(v, 15)
    with open(output_file, "w") as f:
        f.write("kmer\tnegLog10\tbeta\tmac\n")
        for r in rows:
            f.write("\t".join(text(v) for v in r) + "\n")


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


# --------------------------------------------------------------------------
# Figure data for plot_figures.R (PLAN 5.5)
# --------------------------------------------------------------------------

FIGURE_DATA_VERSION = "kmer_pipeline figure data v1"
NA_SENTINEL = "__NA__"


def _fd_value(v, typ):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return NA_SENTINEL
    if typ == "numeric":
        v = float(v)
        if math.isinf(v):
            return "Inf" if v > 0 else "-Inf"
        return "%.17g" % v  # reads back exactly
    if typ == "integer":
        return str(int(v))
    if typ == "logical":
        return "TRUE" if v else "FALSE"
    v = str(v)
    if "\t" in v or "\n" in v or "\r" in v or v == NA_SENTINEL:
        r_stop("Error: value ", repr(v), " can't be written to a figure data file")
    return v


class FigureData:
    """The data behind every figure, written for plot_figures.R into
    <figures_dir>figure_data/: params.tsv (key, value), tables as gzipped
    tab-separated files whose first line names the format and whose second line
    lists the columns as name:type (character, numeric, integer, logical), and
    expected_figures.txt, the figures R must draw (and no others)."""

    def __init__(self, figures_dir):
        self.dir = figures_dir + "figure_data/"
        if not os.path.isdir(self.dir):
            rcompat.r_dir_create(self.dir)
        self.params = {}
        self.expected = []

    def param(self, key, value):
        if isinstance(value, bool):
            value = "TRUE" if value else "FALSE"
        elif isinstance(value, float):
            value = _fd_value(value, "numeric")
        elif value is None:
            value = NA_SENTINEL
        self.params[key] = str(value)

    def table(self, name, columns, rows):
        """columns: list of (name, type); rows: iterable of tuples."""
        import gzip
        types = [t for _, t in columns]
        with gzip.open(self.dir + name + ".tsv.gz", "wt", encoding="utf-8") as f:
            f.write("# " + FIGURE_DATA_VERSION + "\n")
            f.write("\t".join(n + ":" + t for n, t in columns) + "\n")
            f.write("".join("\t".join(_fd_value(v, t) for v, t in zip(r, types)) + "\n" for r in rows))

    def expect(self, path):
        if path not in self.expected:
            self.expected.append(path)

    def close(self):
        with open(self.dir + "params.tsv", "w") as f:
            f.write("# " + FIGURE_DATA_VERSION + "\n")
            f.write("".join(k + "\t" + v + "\n" for k, v in self.params.items()))
        with open(self.dir + "expected_figures.txt", "w") as f:
            f.write("".join(p + "\n" for p in self.expected))
