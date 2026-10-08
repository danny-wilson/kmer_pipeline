"""Close-up alignment figures and tables for the top genes: BLAST of each gene's
k-mers against the gene region, the BLAST result tables, alignment figures in
sliding windows, per-gene Manhattan plots, and the table of significant k-mers
per alignment figure. Port of alignmentfunctions.R.

The tables are exact ports. The figures are drawn by plot_figures.R (PLAN 5.5)
from the figure data written here: each gene's processed BLAST tables, the k-mers
without a result, the reference region and its translations, and the GenBank
features. The decisions about which figures exist are kept here, as R makes them,
to list the figures R must draw."""
import math
import os
import sys

import numpy as np

import rcompat
from rcompat import r_cat, r_paste0, r_stop

# Figure data for plot_figures.R (set by plot_closeup_alignments)
FIGURES = None

# oneLetterCodes as alignmentfunctions.R defines it: unknown codons are "X"
oneLetterCodes = {"Gly": "G", "Ala": "A", "Leu": "L", "Met": "M", "Phe": "F", "Trp": "W", "Lys": "K", "Gln": "Q",
                  "Glu": "E", "Ser": "S", "Pro": "P", "Val": "V", "Ile": "I", "Cys": "C", "Tyr": "Y", "His": "H",
                  "Arg": "R", "Asn": "N", "Asp": "D", "Thr": "T", "STO": "*", "---": "X"}

col_lib_nuc = {"A": "#009E73", "C": "#0072B2", "G": "black", "T": "#E69F00"}
col_lib_pro = {"*": "#000000", "A": "#009E73", "C": "#0072B2", "D": "#E69F00", "E": "#A01FF0", "F": "#50FF00",
               "G": "#FAC0CB", "H": "#F8A503", "I": "#ADD8E6", "K": "#0C008B", "L": "#8B0000", "M": "#1A6400",
               "N": "#a52a2a", "P": "#ffbbff", "Q": "#F78C02", "R": "#df9797", "S": "#90EE90", "T": "#FDFF00",
               "V": "#9d9d00", "W": "#ff0f39", "Y": "#A52A29", "-": "#ffffff"}
colour_selection = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00"]

# --------------------------------------------------------------------------
# Small R helpers
# --------------------------------------------------------------------------


def pstr(x):
    """A value as paste() writes it: NA (None) becomes "NA"."""
    return "NA" if x is None else x


def r_substr(x, start, stop):
    """substr(x, start, stop): 1-based inclusive; start below 1 counts from 1;
    NA (None) stays NA."""
    if x is None:
        return None
    start = max(int(start), 1)
    stop = int(stop)
    if stop < start:
        return ""
    return x[start - 1:stop]


def rc_full_str(x):
    """paste(rc_full(strsplit(x, "")), collapse = ""): unknown bases become "NA"."""
    tab = {"A": "T", "C": "G", "G": "C", "T": "A", "-": "-", "N": "N"}
    return "".join(tab.get(b, "NA") for b in reversed(x))


def nmatch(a, b):
    """length(which(strsplit(a) == strsplit(b))) with R's recycling; NA (None)
    compares as NA, so never matches."""
    if a is None or b is None or not a or not b:
        return 0
    n = max(len(a), len(b))
    return sum(1 for k in range(n) if a[k % len(a)] == b[k % len(b)])


class Table:
    """A small data frame: ordered columns of Python values, as R's read.table
    types them (int, float, str)."""

    def __init__(self, columns, rows):
        self.columns = list(columns)
        self.rows = [list(r) for r in rows]

    def col(self, name):
        k = self.columns.index(name)
        return [r[k] for r in self.rows]

    def set_col(self, name, values):
        if name not in self.columns:
            self.columns.append(name)
            for r, v in zip(self.rows, values):
                r.append(v)
        else:
            k = self.columns.index(name)
            for r, v in zip(self.rows, values):
                r[k] = v

    def subset(self, idx):
        return Table(self.columns, [self.rows[k] for k in idx])

    def __len__(self):
        return len(self.rows)

    def write(self, path, col_names=True):
        with open(path, "w") as f:
            if col_names:
                f.write("\t".join(self.columns) + "\n")
            for r in self.rows:
                f.write("\t".join("NA" if v is None else (v if isinstance(v, str) else rcompat.r_str(v, 15))
                                  for v in r) + "\n")


def read_table_header(path):
    """read.table(path, h = T, sep = "\t") as a Table."""
    df = rcompat.r_read_table(path, header=True, sep="\t")
    rows = []
    for rec in df.astype(object).itertuples(index=False, name=None):
        rows.append([None if rcompat._is_na(v) else (int(v) if isinstance(v, (int, np.integer)) and not isinstance(v, bool)
                                                      else (float(v) if isinstance(v, (float, np.floating)) else v))
                     for v in rec])
    return Table(df.columns, rows)


def num(v):
    return math.nan if v is None else float(v)


# --------------------------------------------------------------------------
# Gene look-up and reference region
# --------------------------------------------------------------------------


def create_gene_lookup(ref, ref_length, recs=None):
    """sequence_functions.R's create_gene_lookup()$gene_lookup: rows (name, id,
    start, end, strand) for genes then intergenic regions. With several reference
    records (recs, D6), the regions of reference.regions (a wrap-round region's row
    gives its first range, after the record's last gene, as for one record)."""
    if recs is not None and len(recs) > 1:
        import reference
        regs = reference.regions(ref, recs)
        return [(r.name, str(k + 1), r.ranges[0][0], r.ranges[0][1], r.strand) for k, r in enumerate(regs)]
    names = list(ref["name"])
    starts = [float(v) for v in ref["start"]]
    ends = [float(v) for v in ref["end"]]
    strands = [float(v) for v in ref["strand"]]
    intergenic, istart, iend = [], [], []
    n = len(names)
    for i in rcompat.r_colon(2, n):
        if i > n:
            r_stop("Error in if (inter_start <= inter_end): missing value where TRUE/FALSE needed")
        inter_start = ends[i - 2] + 1
        inter_end = starts[i - 1] - 1
        if inter_start <= inter_end:
            intergenic.append(names[i - 2] + ":" + names[i - 1])
            istart.append(inter_start)
            iend.append(inter_end)
        if i == n:
            intergenic.append(names[-1] + ":")
            istart.append(ends[n - 1] + 1)
            iend.append(float(ref_length))
    allnames = names + intergenic
    return [(nm, str(k + 1), s, e, st) for k, (nm, s, e, st) in
            enumerate(zip(allnames, starts + istart, ends + iend, strands + [1.0] * len(intergenic)))]


def get_ref_gene_i(gene_lookup, genename_i, ref_length, ref_fa, oneLetterCodes, bounds=None):
    """The gene's region of the reference, with 999 bases either side, within bounds (the
    gene's record, D6; the whole reference by default)."""
    import sequence_functions as sf
    wh = [k for k, row in enumerate(gene_lookup) if row[0] == genename_i]
    if not wh:
        r_stop("Error: no match for gene name", genename_i, "\n")
    row = gene_lookup[wh[0]]
    lo, hi = bounds if bounds is not None else (1.0, float(ref_length))
    ref_start_i = row[2] - 999
    if ref_start_i < lo:
        ref_start_i = float(lo)
    ref_end_i = row[3] + 999
    if ref_end_i > hi:
        ref_end_i = float(hi)
    ref_gene_i = rcompat.r_paste_collapse(rcompat.r_index(ref_fa, rcompat.r_colon(int(ref_start_i), int(ref_end_i))))
    length_protein = len(rcompat.r_colon(int(row[2]), int(row[3]))) / 3
    all_translations = translate_6_frames_alignment(ref_gene_i, oneLetterCodes)
    correct_frame = 1 if row[4] == 1 else 4
    return {"ref_start_i": ref_start_i, "ref_end_i": ref_end_i, "ref_gene_i": ref_gene_i,
            "length_protein": length_protein, "all_translations": all_translations, "correct_frame": correct_frame,
            "length_correct": len(all_translations[correct_frame - 1])}


def translate_function_alignment(contig, frame, oneLetterCodes):
    import sequence_functions as sf
    seq = list(contig)
    if frame > 3:
        seq = [sf.revcompl_full.get(b, "NA") for b in reversed(seq)]
        frame -= 3
    part = ["NA" if b is None else b for b in rcompat.r_index(seq, rcompat.r_colon(frame, len(seq)))]
    return sf.one_letter(sf.translate(sf.totriplet(part)), oneLetterCodes)


def translate_6_frames_alignment(contig, oneLetterCodes):
    return [translate_function_alignment(contig, f, oneLetterCodes) for f in range(1, 7)]


# --------------------------------------------------------------------------
# k-mers, BLAST and result tables
# --------------------------------------------------------------------------


def read_kmer_file(kmerfile, i, prefix, kmer_type, kmer_length, output_dir):
    t = read_table_header(kmerfile)
    seen, rows = set(), []
    for r in t.rows:  # unique() of the data frame rows
        key = tuple(r)
        if key not in seen:
            seen.add(key)
            rows.append(r)
    t = Table(t.columns, rows)
    o = rcompat.r_order([num(v) for v in t.col("negLog10")], decreasing=True)
    t = t.subset(o)
    lines = []
    for k, km in enumerate(t.col("kmer")):
        lines += [">kmer" + str(k + 1), km]
    rcompat.r_cat_lines(lines, r_paste0(output_dir, prefix, "_", kmer_type, kmer_length, "_gene_", i, "_tmp_kmer_file.txt"))
    return t


BLAST_COLS = ["qseqid", "sseqid", "sacc", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "evalue",
              "sstart", "send", "sseq", "qseq"]


def _blast_table(blast_output_file, ncol, perident, kmers_gene_i, nsamples, extra_cols):
    with open(blast_output_file) as f:
        tokens = f.read().split()
    if not tokens:
        return None
    if len(tokens) % ncol:
        r_stop("Error in matrix(blast.search, ncol = ", ncol, "): BLAST output is not a multiple of ", ncol, " fields")
    rows = [tokens[k:k + ncol] for k in range(0, len(tokens), ncol)]
    rows = [r for r in rows if rcompat.r_as_numeric(r[3]) >= perident]
    out = []
    for r in rows:
        idx = int(rcompat.r_as_numeric(r[0][4:]))
        out.append(list(kmers_gene_i.rows[idx - 1]) + r)
    t = Table(["kmer", "negLog10", "beta", "mac"] + BLAST_COLS + extra_cols, out)
    t.set_col("maf", [num(m) / nsamples for m in t.col("mac")])
    return t


def get_top_genes(input_dir, prefix, ngenes, nsamples, bonferroni, ref_name, kmer_type=None, kmer_length=None):
    """The top-gene k-mer files written by plotManhattan, in gene order, and the
    gene names. Not called (plot_closeup_alignments is given them), as in R."""
    genes = rcompat.ls(r_paste0(input_dir, prefix, "_", kmer_type, kmer_length, "_", ref_name,
                                                     "*topgene_*kmersandpvals.txt"))
    order = [float(g.split("/")[-1].split("topgene_")[1].split("_")[0]) for g in genes]
    genes = [genes[k] for k in rcompat.r_order(order)]
    genes_names = [g.split("_topgene_" + str(x + 1) + "_")[1].split("_kmersandpvals.txt")[0] for x, g in enumerate(genes)]
    r_cat("Read in top gene names", "\n")
    return {"genes": genes, "genes_names": genes_names, "genes_all": list(genes)}


def process_blast_protein(prefix, i, kmers_gene_i, genes_names, j, blastPath, correct_frame, genename_i, all_translations,
                          kmer_type, kmer_length, perident, nsamples, output_dir, ref_name):
    correct_or_wrong = "correct_frame" if j == correct_frame else "wrong_frame"
    header = r_paste0(">", genename_i, "_", correct_or_wrong, "_", j)
    ref_trans_file = r_paste0(output_dir, prefix, "_", kmer_type, kmer_length, "_", genename_i, "_", correct_or_wrong, "_", j,
                              "_pro.fa")
    rcompat.r_cat_lines([header, all_translations[j - 1]], ref_trans_file)
    kmer_sequence_file = r_paste0(output_dir, prefix, "_", kmer_type, kmer_length, "_gene_", i, "_tmp_kmer_file.txt")
    blast_output_file = r_paste0(output_dir, prefix, "_", kmer_type, kmer_length, "_blast_out_", i, ".txt")
    rcompat.r_system(blastPath + " -query " + kmer_sequence_file + " -subject " + ref_trans_file
                     + " -max_hsps 1 -outfmt '6 qseqid sseqid sacc pident length mismatch gapopen qstart qend evalue "
                       "sstart send sseq qseq' -evalue 200000 -word_size 2 -gapopen 9 -gapextend 1 -matrix PAM30 -out "
                     + blast_output_file)
    t = _blast_table(blast_output_file, 14, perident, kmers_gene_i, nsamples, [])
    if t is None:
        r_stop("Error: object 'blast_output_file_processed' not found (no BLAST hits for ", genename_i, " frame ", j, ")")
    processed = r_paste0(output_dir, prefix, "_", kmer_type, kmer_length, "_", ref_name, "_top_gene_", i, "_",
                         genes_names[i - 1], "_", j, "_", correct_or_wrong, "_blast_results.txt")
    t.write(processed)
    rcompat.remove(ref_trans_file)
    if j == 6:
        rcompat.remove(kmer_sequence_file)
    rcompat.remove(blast_output_file)
    return processed


def process_blast_nucleotide(blastPath, prefix, kmer_type, kmer_length, i, perident, kmers_gene_i, genes_names, ref_gene_i,
                             genename_i, nsamples, output_dir, ref_name):
    ref_file = r_paste0(output_dir, prefix, "_", kmer_type, kmer_length, "_", genename_i, "_nuc.fa")
    rcompat.r_cat_lines([">" + genename_i, ref_gene_i], ref_file)
    kmer_sequence_file = r_paste0(output_dir, prefix, "_", kmer_type, kmer_length, "_gene_", i, "_tmp_kmer_file.txt")
    blast_output_file = r_paste0(output_dir, prefix, "_", kmer_type, kmer_length, "_blast_out_", i, ".txt")
    rcompat.r_system(blastPath + " -query " + kmer_sequence_file + " -subject " + ref_file
                     + " -max_hsps 1 -outfmt '6 qseqid sseqid sacc pident length mismatch gapopen qstart qend evalue "
                       "sstart send sseq qseq sstrand' -evalue 0.05 -gapopen 5 -gapextend 2 -penalty -3 -reward 2 "
                       "-word_size 4 -out " + blast_output_file)
    t = _blast_table(blast_output_file, 15, perident, kmers_gene_i, nsamples, ["sstrand"])
    if t is None:
        r_stop("Error: object 'blast_output_file_processed' not found (no BLAST hits for ", genename_i, ")")
    processed = r_paste0(output_dir, prefix, "_", kmer_type, kmer_length, "_", ref_name, "_top_gene_", i, "_",
                         genes_names[i - 1], "_blast_results.txt")
    t.write(processed)
    rcompat.remove(ref_file)
    rcompat.remove(kmer_sequence_file)
    rcompat.remove(blast_output_file)
    return processed


def _extend(res, refseq):
    """The shared part of fix_res_table_protein/_nucleotide: extend each BLAST
    hit so the whole k-mer is aligned."""
    kmer = res.col("kmer")
    klen = [len(k) for k in kmer]
    sstart = [num(v) for v in res.col("sstart")]
    send = [num(v) for v in res.col("send")]
    qstart = [num(v) for v in res.col("qstart")]
    qend = [num(v) for v in res.col("qend")]
    sseq = list(res.col("sseq"))
    qseq = list(res.col("qseq"))

    w = [k for k in range(len(res)) if qstart[k] > 1]
    new_start = {k: sstart[k] for k in w}
    for k in w:
        sstart[k] = new_start[k] - qstart[k] + 1
    for k in range(len(res)):
        if sstart[k] < 1:
            sstart[k] = 1.0
    for k in w:
        if sstart[k] != new_start[k]:
            o, s, q = new_start[k], sstart[k], qstart[k]
            qseq[k] = pstr(r_substr(kmer[k], q - len(rcompat.r_colon(int(s), int(o - 1))), q - 1)) + pstr(qseq[k])
            lo, hi = sorted((o - q + 1, o - 1))
            sseq[k] = pstr(r_substr(refseq, lo, hi)) + pstr(sseq[k])

    w = [k for k in range(len(res)) if qend[k] < klen[k]]
    new_end = {k: send[k] for k in w}
    for k in w:
        send[k] = new_end[k] - qend[k] + klen[k]
    for k in range(len(res)):
        if send[k] > len(refseq):
            send[k] = float(len(refseq))
    for k in w:
        if send[k] != new_end[k]:
            o, s, q, kl = new_end[k], send[k], qend[k], klen[k]
            sseq[k] = pstr(sseq[k]) + pstr(r_substr(refseq, o + 1, o - q + kl))
            qseq[k] = pstr(qseq[k]) + pstr(r_substr(kmer[k], q + 1, q + len(rcompat.r_colon(int(o + 1), int(s)))))
    res.set_col("sstart", sstart)
    res.set_col("send", send)
    res.set_col("sseq", sseq)
    res.set_col("qseq", qseq)


def fix_res_table_protein(res, refseq):
    _extend(res, refseq)
    res.set_col("origkmer", list(res.col("kmer")))
    return res


def fix_res_table_nucleotide(res, ref_fa, kmer_len=None):
    res.set_col("origkmer", list(res.col("kmer")))
    # For those which aligned to the reverse strand, reverse complement the sequences
    strand = res.col("sstrand")
    kmer = res.col("kmer")
    for k in range(len(res)):
        if strand[k] == "minus":
            r = res.rows[k]
            c = res.columns.index
            klen = len(kmer[k])
            qs, qe = num(r[c("qstart")]), num(r[c("qend")])
            ss, se = r[c("sstart")], r[c("send")]
            r[c("sseq")] = rc_full_str(r[c("sseq")])
            r[c("qseq")] = rc_full_str(r[c("qseq")])
            r[c("kmer")] = rc_full_str(r[c("kmer")])
            r[c("sstart")], r[c("send")] = num(se), num(ss)
            r[c("qstart")], r[c("qend")] = klen - qe + 1, klen - qs + 1
    _extend(res, ref_fa)
    return res


def read_res_table(resfile, refseq, i, genes_names, perident=70, kmer_type=None, kmer_length=None):
    res = read_table_header(resfile)
    if len(res) == 0:
        r_stop("Error in if (nchar(refseq) > nchar(res$kmer[1])): no BLAST results in ", resfile)
    res = res.subset(rcompat.r_order([num(v) for v in res.col("negLog10")], decreasing=True))
    if kmer_type == "protein":
        res = fix_res_table_protein(res, refseq)
    else:
        res = fix_res_table_nucleotide(res, refseq, kmer_length)
    nm = [nmatch(a, b) for a, b in zip(res.col("sseq"), res.col("qseq"))]
    k1 = len(res.col("kmer")[0])
    thr = (k1 if len(refseq) > k1 else len(refseq)) * (perident / 100)
    low = res.subset([k for k in range(len(res)) if nm[k] < thr])
    res = res.subset([k for k in range(len(res)) if nm[k] >= thr])
    return {"res": res, "low_match_results_j": low}


def get_kmers_noresult(gene_i_results_list, kmers_gene_i, prefix, i, genes_names, kmer_type, kmer_length, output_dir,
                       ref_name, nsamples):
    matches = set()
    for res in gene_i_results_list:
        matches.update(res.col("origkmer"))
    w = [k for k, km in enumerate(kmers_gene_i.col("kmer")) if km not in matches]
    if not w:
        return None
    t = kmers_gene_i.subset(w)
    t.columns = ["kmer", "negLog10", "beta", "mac"]
    # MAC is always present, so MAF = MAC / n (as for the BLAST-matched k-mers); without it the MAF
    # threshold plots silently drop these k-mers and can be left with nothing to plot
    t.set_col("maf", [num(m) / nsamples for m in t.col("mac")])
    t = t.subset(rcompat.r_order([num(v) for v in t.col("negLog10")], decreasing=True))
    t.write(r_paste0(output_dir, prefix, "_", kmer_type, kmer_length, "_", ref_name, "_top_gene_", i, "_", genes_names[i - 1],
                     "_no_blast_result_or_poor_alignment.txt"))
    return t


# --------------------------------------------------------------------------
# Alignment figures (drawn by plot_figures.R; here, which ones exist)
# --------------------------------------------------------------------------


def _which_to_align(res, genestart, geneend, macormaf, minor_allele_threshold):
    send, sstart, ma = res.col("send"), res.col("sstart"), res.col(macormaf)
    return [k for k in range(len(res)) if num(send[k]) >= genestart and num(sstart[k]) <= geneend
            and num(ma[k]) >= minor_allele_threshold]


def _out_rows(res, wta, bonferroni, override_signif):
    neg = [num(v) for v in res.col("negLog10")]
    n_sig = sum(1 for k in wta if neg[k] >= bonferroni)
    take = wta if override_signif else [x for x in rcompat.r_index(wta, rcompat.r_colon(1, n_sig)) if x is not None]
    cols = ["kmer", "negLog10", "beta", "mac", "maf", "sstart"]
    return [[res.rows[k][res.columns.index(c)] for c in cols] for k in take]


def plot_alignment_function_nucleotide(genestart, geneend, minor_allele_threshold, res, nsamples, maname, bonferroni,
                                       ref_fa, prefix, gene_name, p, plot_ref, main, reverse_xaxis_start,
                                       forward_xaxis_start, override_signif, alignment_range, macormaf, output_dir,
                                       kmer_type, kmer_length, ref_name):
    alignment_start, alignment_end = int(alignment_range[0]), int(alignment_range[1])
    nlen = len(rcompat.r_colon(int(genestart), int(geneend)))
    if reverse_xaxis_start is None:
        plot_start_position = forward_xaxis_start
        plot_end_position = float(forward_xaxis_start + nlen - 1)
    else:
        plot_start_position = reverse_xaxis_start
        plot_end_position = float(reverse_xaxis_start - nlen + 1)
    wta = _which_to_align(res, genestart, geneend, macormaf, minor_allele_threshold)
    if not wta:
        return None
    neg = [num(v) for v in res.col("negLog10")]
    bonferroni_lim = sum(1 for k in wta if neg[k] < bonferroni)
    if not (bonferroni_lim != len(wta) or override_signif):
        return None
    FIGURES.expect(r_paste0(output_dir, prefix, "_", kmer_type, kmer_length, "_", ref_name, "_", gene_name, "_plot_", p,
                            "_pos_", plot_start_position, "_to_", plot_end_position, maname, "_alignment.png"))
    rows = _out_rows(res, wta, bonferroni, override_signif)
    span = rcompat.r_colon(alignment_start, alignment_end) if reverse_xaxis_start is None else \
        rcompat.r_colon(alignment_end, alignment_start)
    for r in rows:
        x = int(num(r[5]))
        r[5] = span[x - 1] if 1 <= x <= len(span) else None
    plot = r_paste0(gene_name, "_plot_", p, "_ps_", plot_start_position, "_to_", plot_end_position)
    return [[gene_name, plot] + r for r in rows]


def run_alignment_nplots_nucleotide(ref_gene_i, res, nsamples, bonferroni, prefix, gene_name, override_signif, gene_lookup,
                                    wh_genelookup, minor_allele_threshold, macormaf, output_dir, kmer_type, kmer_length,
                                    ref_name):
    ref_start_i, ref_end_i = int(ref_gene_i["ref_start_i"]), int(ref_gene_i["ref_end_i"])
    seq = ref_gene_i["ref_gene_i"]
    nplots = _windows(len(seq), 40, 99, max(len(k) for k in res.col("kmer")) if len(res) else 0)
    out = []
    for p, (gs, ge) in enumerate(nplots, start=1):
        if gene_lookup[wh_genelookup][4] != 1:
            rev = rcompat.r_index(rcompat.r_colon(ref_end_i, ref_start_i), [int(gs)])[0]
            fwd = None
        else:
            rev = None
            fwd = rcompat.r_index(rcompat.r_colon(ref_start_i, ref_end_i), [int(gs)])[0]
        o = plot_alignment_function_nucleotide(gs, ge, 0.0, res, nsamples, "", bonferroni, seq, prefix, gene_name, p, False,
                                               "All nucleotide kmers", rev, fwd, override_signif, (ref_start_i, ref_end_i),
                                               macormaf, output_dir, kmer_type, kmer_length, ref_name)
        plot_alignment_function_nucleotide(gs, ge, minor_allele_threshold, res, nsamples,
                                           r_paste0("_", macormaf, minor_allele_threshold), bonferroni, seq, prefix,
                                           gene_name, p, False,
                                           r_paste0("Nucleotide kmers ≥ ", macormaf, " ", minor_allele_threshold),
                                           rev, fwd, override_signif, (ref_start_i, ref_end_i), macormaf, output_dir,
                                           kmer_type, kmer_length, ref_name)
        if o is not None:
            out += o
    return out


def _windows(n, step, width, max_klen):
    """nplots: windows (start, end) every `step` positions, `width`+1 long, the
    last cut at n; the last dropped if shorter than the longest k-mer and the one
    before already reaches n."""
    starts = [1 + step * k for k in range(int((n - 1) // step) + 1)] if n >= 1 else [1]
    w = [(float(s), float(min(s + width, n))) for s in starts]
    if len(w) > 1:
        last = w[-1]
        if len(rcompat.r_colon(int(last[0]), int(last[1]))) < max_klen and w[-2][1] == n:
            w = w[:-1]
    return w


def plot_alignment_function_protein(genestart, geneend, minor_allele_threshold, res, nsamples, maname, bonferroni,
                                    translation, prefix, gene_name, correct_or_wrong, j, p, plot_ref, main, x_adjust,
                                    override_signif, macormaf, output_dir, kmer_type, kmer_length, ref_name):
    wta = _which_to_align(res, genestart, geneend, macormaf, minor_allele_threshold)
    if not wta:
        return None
    neg = [num(v) for v in res.col("negLog10")]
    bonferroni_lim = sum(1 for k in wta if neg[k] < bonferroni)
    if not (bonferroni_lim != len(wta) or override_signif):
        return None
    FIGURES.expect(r_paste0(output_dir, prefix, "_", kmer_type, kmer_length, "_", ref_name, "_", gene_name, "_",
                            correct_or_wrong, "_", j, "_plot_", p, "_aminoacids_", genestart - x_adjust, "_to_",
                            geneend - x_adjust, maname, "_alignment.png"))
    if correct_or_wrong != "correct_frame":
        return None
    rows = _out_rows(res, wta, bonferroni, override_signif)
    for r in rows:
        r[5] = num(r[5]) - x_adjust
    plot = r_paste0(gene_name, "_plot_", p, "_aminoacids_", genestart - x_adjust, "_to_", geneend - x_adjust)
    return [[gene_name, plot] + r for r in rows]


def run_alignment_nplots_protein(ref_gene_i, res, nsamples, bonferroni, prefix, gene_name, j, override_signif,
                                 minor_allele_threshold, macormaf, output_dir, kmer_type, kmer_length, ref_name):
    correct_or_wrong = "correct_frame" if j == ref_gene_i["correct_frame"] else "wrong_frame"
    tr = ref_gene_i["all_translations"][j - 1]
    nplots = _windows(len(tr), 20, 39, max(len(k) for k in res.col("kmer")) if len(res) else 0)
    out = []
    for p, (gs, ge) in enumerate(nplots, start=1):
        o = plot_alignment_function_protein(gs, ge, 0.0, res, nsamples, "", bonferroni, tr, prefix, gene_name,
                                            correct_or_wrong, j, p, False, "All protein kmers", 333, override_signif,
                                            macormaf, output_dir, kmer_type, kmer_length, ref_name)
        plot_alignment_function_protein(gs, ge, minor_allele_threshold, res, nsamples,
                                         r_paste0("_", macormaf, minor_allele_threshold), bonferroni, tr, prefix,
                                         gene_name, correct_or_wrong, j, p, False,
                                         r_paste0("Protein kmers ≥ ", macormaf, " ", minor_allele_threshold), 333,
                                         override_signif, macormaf, output_dir, kmer_type, kmer_length, ref_name)
        if o is not None:
            out += o
    return out if j == ref_gene_i["correct_frame"] else None


# --------------------------------------------------------------------------
# Per-gene Manhattan plots
# --------------------------------------------------------------------------


def _manhattan_values(res, which_kmers_no_result, macormaf, minor_allele_threshold):
    """The -log10 p values plotted in a gene's Manhattan plot, and which pass the MAF/MAC
    threshold. The k-mers without a BLAST result carry both mac and maf, so they are
    thresholded like the rest (the original R dropped them under a MAF threshold)."""
    ypos = [num(v) for v in res.col("negLog10")]
    ma = [num(v) for v in res.col(macormaf)]
    if which_kmers_no_result is not None:
        ypos += [num(v) for v in which_kmers_no_result.col("negLog10")]
        ma += [num(v) for v in which_kmers_no_result.col(macormaf)]
    which = [k for k in range(len(ma)) if ma[k] >= minor_allele_threshold]
    return ypos, which


def run_manhattan_single(res, which_kmers_no_result, prefix_path, gene_name, bonferroni, macormaf, minor_allele_threshold,
                         name_part):
    """The files R's run_manhattan_single_protein/_nucleotide write: all k-mers and those
    passing the threshold, each with a ylim 50 version when the maximum exceeds 100."""
    ypos, which = _manhattan_values(res, which_kmers_no_result, macormaf, minor_allele_threshold)
    for subset, maname in ((list(range(len(ypos))), "_allkmers"), (which, r_paste0("_", macormaf, minor_allele_threshold))):
        y = [ypos[k] for k in subset]
        if not y:
            if name_part:  # protein: R plots only if there are points, and png() then writes no file
                continue
            r_cat("Warning: no k-mers to plot for ", gene_name, maname, "; figure skipped\n")
            continue
        FIGURES.expect(prefix_path + "_" + gene_name + name_part + "_Manhattan" + maname + ".png")
        if max(y) > 100:
            FIGURES.expect(prefix_path + "_" + gene_name + name_part + "_Manhattan_ylim50" + maname + ".png")


def run_manhattan_allframes(gene_i_results_list, prefix_path, gene_name, which_kmers_no_result, minor_allele_threshold,
                            macormaf):
    """The files R's run_manhattan_allframes writes (always drawn; ylim 50 versions when
    plot_allframes_manhattan's maximum exceeds 100)."""
    yall = [num(v) for res in gene_i_results_list for v in res.col("negLog10")]
    ythr = [num(v) for res in gene_i_results_list for v, m in zip(res.col("negLog10"), res.col(macormaf))
            if num(m) >= minor_allele_threshold]
    if which_kmers_no_result is not None:
        yall += [num(v) for v in which_kmers_no_result.col("negLog10")]
        ythr += [num(v) for v, m in zip(which_kmers_no_result.col("negLog10"), which_kmers_no_result.col(macormaf))
                 if num(m) >= minor_allele_threshold]
    for y, maname in ((yall, "_allkmers"), (ythr, r_paste0("_", macormaf, minor_allele_threshold))):
        FIGURES.expect(prefix_path + "_" + gene_name + "_allframes_Manhattan" + maname + ".png")
        if y and max(y) > 100:
            FIGURES.expect(prefix_path + "_" + gene_name + "_allframes_Manhattan_ylim50" + maname + ".png")


STRING_COLUMNS = {"kmer", "qseqid", "sseqid", "sacc", "sseq", "qseq", "sstrand", "origkmer"}


def write_table_for_r(fd, name, table):
    """A Table (as R's data frame) for plot_figures.R: strings as character, all else numeric."""
    cols = [(c, "character" if c in STRING_COLUMNS else "numeric") for c in table.columns]
    fd.table(name, cols, table.rows)


def write_features(fd, ref_gb):
    """The GenBank features R's gene-arrow plots draw (alignmentfunctions.R:1724)."""
    import sequence_functions as sf
    t = sf.read_dna_seg_from_file(ref_gb, tagsToParse=("CDS", "repeat_region", "tRNA", "rRNA", "ncRNA"))
    numeric = {"start", "end", "strand", "length", "lty", "lwd", "pch", "cex"}
    cols = [(c, "numeric" if c in numeric else "character") for c in t.columns]
    fd.table("features", cols, (tuple(None if (isinstance(v, float) and v != v) else v for v in r)
                                for r in t.astype(object).itertuples(index=False, name=None)))


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------


def plot_closeup_alignments(ref, ref_length, ref_gb, ref_fa, figures_dir, output_prefix, ngenes, nsamples, bonferroni,
                            gene_lookup, oneLetterCodes, kmer_type, kmer_length, blastPath, perident, ref_name,
                            alignmenttype, genes_all, override_signif, minor_allele_threshold, macormaf, figure_data,
                            correct_only=True):
    """BLAST each top gene's k-mers against its region, write the BLAST tables and the
    table of significant k-mers per alignment plot, and write the figure data R needs
    to draw the alignment and gene Manhattan figures."""
    global FIGURES
    import sequence_functions as sf
    FIGURES = figure_data
    write_features(figure_data, ref_gb)
    import reference
    recs = reference.records(ref_gb)
    gene_lookup = create_gene_lookup(ref, ref_length, recs)
    ref_fa_seq = sf.read_reference(ref_fa)

    genes = genes_all["genes"]
    genes_names = genes_all["genes_names"]
    out_rows = []
    out_cols = ["gene", "plot", "kmer", "negLog10", "beta", "mac", "maf", "ps"]
    gene_rows = []

    for i in range(1, len(genes_names) + 1):
        genename_i = genes_names[i - 1]
        wh_genelookup = [k for k, row in enumerate(gene_lookup) if row[0] == genename_i]
        if not wh_genelookup:
            r_stop("Error: no match for gene name", genename_i, "\n")
        wh_genelookup = wh_genelookup[0]
        bounds = None
        if len(recs) > 1:
            rec = reference.record_of(recs, gene_lookup[wh_genelookup][2])
            bounds = (float(rec.start), float(rec.end))
        ref_gene_i = get_ref_gene_i(gene_lookup, genename_i, ref_length, ref_fa_seq, oneLetterCodes, bounds)
        if kmer_type == "nucleotide" and ref_gene_i["correct_frame"] != 1:
            ref_gene_i["ref_gene_i"] = rc_full_str(ref_gene_i["ref_gene_i"])

        kmers_gene_i = read_kmer_file(genes[i - 1], i, output_prefix, kmer_type, kmer_length, figures_dir)

        if kmer_type == "protein":
            gene_i_results_list = []
            for j in range(1, 7):
                blast_search = process_blast_protein(output_prefix, i, kmers_gene_i, genes_names, j, blastPath,
                                                     ref_gene_i["correct_frame"], genename_i, ref_gene_i["all_translations"],
                                                     kmer_type, kmer_length, perident, nsamples, figures_dir, ref_name)
                res = read_res_table(blast_search, ref_gene_i["all_translations"][j - 1], i, genes_names, perident,
                                     kmer_type)["res"]
                gene_i_results_list.append(res)
                if (correct_only and j == ref_gene_i["correct_frame"]) or not correct_only:
                    o = run_alignment_nplots_protein(ref_gene_i, res, nsamples, bonferroni, output_prefix, genes_names[i - 1],
                                                     j, override_signif, minor_allele_threshold, macormaf, figures_dir,
                                                     kmer_type, kmer_length, ref_name)
                    if j == ref_gene_i["correct_frame"] and o:
                        out_rows += o
        else:
            blast_search = process_blast_nucleotide(blastPath, output_prefix, kmer_type, kmer_length, i, perident,
                                                    kmers_gene_i, genes_names, ref_gene_i["ref_gene_i"], genename_i,
                                                    nsamples, figures_dir, ref_name)
            res = read_res_table(blast_search, ref_gene_i["ref_gene_i"], i, genes_names, perident, kmer_type,
                                 kmer_length)["res"]
            gene_i_results_list = [res]
            out_rows += run_alignment_nplots_nucleotide(ref_gene_i, res, nsamples, bonferroni, output_prefix, genename_i,
                                                        override_signif, gene_lookup, wh_genelookup,
                                                        minor_allele_threshold, macormaf, figures_dir, kmer_type,
                                                        kmer_length, ref_name)

        which_kmers_no_result = get_kmers_noresult(gene_i_results_list, kmers_gene_i, output_prefix, i, genes_names,
                                                   kmer_type, kmer_length, figures_dir, ref_name, nsamples)
        prefix_path = r_paste0(figures_dir, output_prefix, "_", kmer_type, kmer_length, "_", ref_name)
        if kmer_type == "protein":
            for j in range(1, 7):
                if (correct_only and j == ref_gene_i["correct_frame"]) or not correct_only:
                    cw = "correct_frame" if j == ref_gene_i["correct_frame"] else "wrong_frame"
                    run_manhattan_single(gene_i_results_list[j - 1], which_kmers_no_result, prefix_path, genes_names[i - 1],
                                         bonferroni, macormaf, minor_allele_threshold, "_" + cw + "_" + str(j))
            run_manhattan_allframes(gene_i_results_list, prefix_path, genes_names[i - 1], which_kmers_no_result,
                                    minor_allele_threshold, macormaf)
        else:
            run_manhattan_single(gene_i_results_list[0], which_kmers_no_result, prefix_path, genes_names[i - 1], bonferroni,
                                 macormaf, minor_allele_threshold, "")

        # Figure data for this gene
        gene_rows.append((i, genename_i))
        fd = figure_data
        with open(fd.dir + "gene_" + str(i) + "_params.tsv", "w") as f:
            import Manhattan_functions as mf
            f.write("# " + mf.FIGURE_DATA_VERSION + "\n")
            for key in ("ref_start_i", "ref_end_i", "length_protein", "correct_frame", "length_correct"):
                f.write(key + "\t" + mf._fd_value(ref_gene_i[key], "numeric") + "\n")
            f.write("strand\t" + mf._fd_value(gene_lookup[wh_genelookup][4], "numeric") + "\n")
        fd.table("gene_" + str(i) + "_sequences", [("name", "character"), ("sequence", "character")],
                 [("region", ref_gene_i["ref_gene_i"])]
                 + [("frame" + str(f), t) for f, t in enumerate(ref_gene_i["all_translations"], start=1)])
        for j, res in enumerate(gene_i_results_list, start=1):
            write_table_for_r(fd, "gene_" + str(i) + "_res_" + str(j), res)
        if which_kmers_no_result is not None:
            write_table_for_r(fd, "gene_" + str(i) + "_no_result", which_kmers_no_result)

    figure_data.table("genes", [("index", "integer"), ("gene", "character")], gene_rows)

    Table(out_cols, out_rows).write(r_paste0(figures_dir, output_prefix, "_", kmer_type, kmer_length, "_", ref_name, "_",
                                             alignmenttype, "_all_top_genes_significant_kmers_per_alignment_plot.txt"))
