"""Sequence functions: reference reading, translation and contig reading.
Port of sequence_functions.R.

R works with character vectors of single bases; here a sequence is a str and a
character vector is a list of str. reorder_reference_gbk and create_gene_lookup
(which need the GenBank parser) are ported with the alignment step. transcribe()
is not called anywhere and is not ported.
"""
import sys

import rcompat
from rcompat import r_colon, r_index, r_stop


def read_reference(ref_file):
    """The reference sequence (every line after the first, concatenated), as a
    str. scan() is not quiet here, so "Read N items" goes to stderr as in R."""
    r = rcompat.r_scan_lines(ref_file, quiet=False)
    if sum(line.startswith(">") for line in r) > 1:
        r_stop("Error: reference fasta file ", ref_file, " contains more than one record; "
               "only single-record references are supported", "\n")
    # r[2:length(r)]: with one line, R's 2:1 gives c(NA, r[1])
    return rcompat.r_paste_collapse(r_index(r, r_colon(2, len(r))))


rev_base = {"A": "T", "T": "A", "C": "G", "G": "C", "-": "-"}

revcompl = {"A": "T", "C": "G", "G": "C", "T": "A"}
revcompl_full = {"A": "T", "C": "G", "G": "C", "T": "A", "-": "-", "N": "N"}


def rc(x):
    """revcompl[rev(x)] for a str of bases. Bases outside the table become None
    (R's NA), so a list is returned."""
    return [revcompl.get(b) for b in reversed(x)]


def rc_full(x):
    return [revcompl_full.get(b) for b in reversed(x)]


def totriplet(x):
    """Codons of a sequence (list or str of bases), dropping any partial codon at
    the end. As R's seq(1, 0, by = 3), fails when there is no whole codon."""
    L = (len(x) // 3) * 3
    if L == 0:
        raise ValueError("wrong sign in 'by' argument (totriplet: fewer than 3 bases)")
    return ["".join(x[i:i + 3]) for i in range(0, L, 3)]


geneticCode = {
    "TTT": "Phe", "TTC": "Phe", "TTA": "Leu", "TTG": "Leu",
    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser",
    "TAT": "Tyr", "TAC": "Tyr", "TAA": "STO", "TAG": "STO",
    "TGT": "Cys", "TGC": "Cys", "TGA": "STO", "TGG": "Trp",
    "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
    "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "CAT": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "ATT": "Ile", "ATC": "Ile", "ATA": "Ile", "ATG": "Met",
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "AAT": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "AGT": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "GAT": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"}

oneLetterCodes = {"Gly": "G", "Ala": "A", "Leu": "L", "Met": "M", "Phe": "F", "Trp": "W", "Lys": "K",
                  "Gln": "Q", "Glu": "E", "Ser": "S", "Pro": "P", "Val": "V", "Ile": "I", "Cys": "C",
                  "Tyr": "Y", "His": "H", "Arg": "R", "Asn": "N", "Asp": "D", "Thr": "T", "STO": "*"}


def translate(codons):
    """translate() for one row of codons: three-letter amino acids, "---" for a
    codon not in the genetic code (after upper-casing). R's oneLetter = TRUE
    option is not used by the pipeline and not ported."""
    return [geneticCode.get(c.upper(), "---") for c in codons]


def one_letter(amino_acids, codes):
    """paste(unlist(codes[amino_acids]), collapse = ""): amino acids missing from
    `codes` are NULL in R's list lookup, so unlist drops them."""
    return "".join(codes[a] for a in amino_acids if a in codes)


def translate_function(contig, frame, oneLetterCodes, revcompl):
    """Translate a contig (str of A, C, G, T) in reading frame 1-3 (forward) or
    4-6 (reverse complement)."""
    contig = list(contig)
    if any(b not in ("A", "C", "G", "T") for b in contig):
        r_stop("Error: contig contains non base characters", "\n")
    if frame > 3:
        contig = [revcompl.get(b) for b in reversed(contig)]
        if any(b is None for b in contig):
            r_stop("Error: reverse complemented contig contains NAs", "\n")
        frame -= 3
    # contig[frame:length(contig)], with R's NA for positions beyond the end
    part = ["NA" if b is None else b for b in r_index(contig, r_colon(frame, len(contig)))]
    return one_letter(translate(totriplet(part)), oneLetterCodes)


def get_output_file(outDir, id):
    return outDir + id + "_translated_all_reading_frames.fa"


def _recycle(*vectors):
    """Lengths after R's recycling of vectors in paste()."""
    n = max(len(v) for v in vectors)
    return 0 if any(len(v) == 0 for v in vectors) else n


def write_proteins_to_file(contig_names, frame, contigs, id, append=False, outDir=None):
    names = [n + "_rf" + str(frame) for n in contig_names]
    n = _recycle(names, contigs)
    out = [names[i % len(names)] + "\n" + contigs[i % len(contigs)] for i in range(n)]
    rcompat.r_cat_lines(out, get_output_file(outDir, id), append=append)


def remove_Ns(contig):
    """Split a contig at runs of N, dropping empty pieces."""
    return [c for c in contig.split("N") if c != ""]


def read_contigs_removeNs(contig_file):
    """Read a (possibly gzipped) multi-FASTA of contigs. Each contig is upper-cased
    and split at Ns; piece j of contig "name" is called name_j. Returns
    (contig_seq, contig_names_final)."""
    contigs = rcompat.r_scan_lines(contig_file, quiet=True)
    contig_name_positions = [i + 1 for i, x in enumerate(contigs) if ">" in x]  # 1-based, as R
    rcompat.r_cat("Number of contigs:", len(contig_name_positions), "\n")
    contig_names = [contigs[p - 1].split(" ")[0] for p in contig_name_positions]
    contig_names_final, contig_seq = [], []

    def lines(a, b):  # contigs[a:b]
        return rcompat.r_paste_collapse(r_index(contigs, r_colon(a, b)))

    n = len(contig_names)
    if n == 0:
        r_stop("Error in contig_names[i]: subscript out of bounds (no contigs in ", contig_file, ")")
    for i in range(1, n + 1):
        if i != n:
            contig_i = lines(contig_name_positions[i - 1] + 1, contig_name_positions[i] - 1)
        else:
            contig_i = lines(contig_name_positions[i - 1] + 1, len(contigs))
        contig_i = contig_i.upper()
        contig_i = remove_Ns(contig_i) if "N" in contig_i else [contig_i]
        contig_seq += contig_i
        contig_names_final += [contig_names[i - 1] + "_" + str(j) for j in r_colon(1, len(contig_i))]
    return contig_seq, contig_names_final


def translate_6_frames(contig_path, id, outDir, oneLetterCodes, revcompl):
    """Translate every contig in all six reading frames and write them, frame by
    frame, to <outDir><id>_translated_all_reading_frames.fa, then gzip it.
    Returns the path of the gzipped file."""
    contig_seq, contig_names_final = read_contigs_removeNs(contig_path)
    rf = [[translate_function(c, frame, oneLetterCodes, revcompl) for c in contig_seq]
          for frame in range(1, 7)]
    for frame in range(1, 7):
        write_proteins_to_file(contig_names_final, frame, rf[frame - 1], id, append=frame > 1, outDir=outDir)
    final_file = get_output_file(outDir, id)
    rcompat.r_system("gzip " + final_file)
    return final_file + ".gz"


# --------------------------------------------------------------------------
# GenBank reading: a port of genoPlotR 0.8.11's read_dna_seg_from_file (the
# GenBank branch) with its helpers, written against its regular expressions so
# that every quirk is kept (PLAN 5.4): CDS only; joins split into exons with
# "intron" rows between; partial (<, >) and single-base locations dropped;
# colons and double quotes removed from qualifiers; name from /gene, else
# /locus_tag; "length" as (end - start + 1)/3 - 1.
# --------------------------------------------------------------------------

import re as _re

_PRINT = r"[^\x00-\x1f\x7f]"        # [[:print:]] (includes space)
_GRAPH = r"[^\x00-\x20\x7f]"        # [[:graph:]]
_ALPHA = r"[^\W\d_]"                # [[:alpha:]]
_SPACE = r"[ \t\n\r\f\v]"           # [[:space:]]


def r_strsplit(s, pattern):
    """strsplit(s, pattern)[[1]] for a regular expression: R drops one empty
    piece at the end, and an empty string gives character(0)."""
    if s == "":
        return []
    parts = _re.split(pattern, s)
    if parts and parts[-1] == "":
        parts.pop()
    return parts


def _grepl(pattern, x):
    return _re.search(pattern, x) is not None


def extract_data(extract, cF):
    vals = [_re.sub(extract, "", el) for el in cF if _re.search("^" + extract, el)]
    return vals[0] if vals else "NA"


def _as_numeric_list(strings):
    return [rcompat.r_as_numeric(s) for s in strings]


def get_start(line):
    if _grepl("complement", line):
        hits = [line] if _re.search(r"^" + _GRAPH + r"+ complement\([0-9]+\.\.[0-9]+\)$", line) else []
    else:
        hits = [line] if _re.search(r"^" + _GRAPH + r"+ [0-9]+\.\.[0-9]+$", line) else []
    return _as_numeric_list(_re.sub(r"_|[ \t]|" + _ALPHA + r"|\(|\)|\.\..*", "", h) for h in hits)


def get_end(line):
    if _grepl("complement", line):
        hits = [line] if _re.search(r"^" + _GRAPH + r"+ complement\([0-9]+\.\.[0-9]+\)$", line) else []
    else:
        hits = [line] if _re.search(r"^" + _GRAPH + r"+ [0-9]+\.\.[0-9]+$", line) else []
    return _as_numeric_list(_re.sub(r"_|[ \t]|" + _ALPHA + r"|\(|\)|.*\.\.", "", h) for h in hits)


def _paste_num(v):
    """paste() of a numeric(0)-or-one-number result: "" when empty."""
    return rcompat.r_as_character(v[0]) if v else ""


_ARTEMIS = ["#FFFFFF", "#646464", "#FF0000", "#00FF00", "#0000FF", "#00FFFF", "#FF00FF", "#FFFF00", "#98FB98",
            "#87CEFA", "#FFA500", "#C89664", "#FFC8C8", "#AAAAAA", "#000000", "#FF3F3F", "#FF7F7F", "#FFBFBF"]


def read_dna_seg_from_file(file, tagsToParse=("CDS",), gene_type="auto"):
    """genoPlotR::read_dna_seg_from_file(file) for a GenBank file. Returns a
    pandas DataFrame with genoPlotR's columns (name, start, end, strand, length,
    pid, gene, synonym, product, proteinid, feature, gene_type, [col,] fill, lty,
    lwd, pch, cex), or None when no features were read. Row labels are R's row
    names (1, 2, ...)."""
    import warnings

    import pandas as pd
    with rcompat.r_open(file) as f:
        importedData = f.read().split("\n")
    if importedData and importedData[-1] == "":
        importedData.pop()
    TYPE = "Unknown"
    if importedData and ">" in importedData[0]:
        TYPE = "Fasta"
    if any(_re.search("^ID", l) for l in importedData):
        TYPE = "EMBL"
    if any(_re.search("^LOCUS", l) for l in importedData):
        TYPE = "Genbank"
    if TYPE == "Unknown":
        r_stop("fileType has to be either 'detect', 'embl', 'genbank' or 'ptt'. Note if file type is .ptt, "
               "please specify this rather than using 'detect'.")
    if TYPE != "Genbank":
        r_stop("Error: only GenBank reference annotation files are supported by the Python pipeline (", file,
               " looks like ", TYPE, ")")

    mainSegments = [k + 1 for k, l in enumerate(importedData) if l[:1].isalnum()]  # 1-based line numbers
    segNames = [_re.sub(r"\*| .*", "", importedData[k - 1]) for k in mainSegments]
    if sum(_grepl("FEATURES|DEFINITION", n) for n in segNames) < 2:
        r_stop("FEATURES or DEFINITION segment missing in GBK File.")
    if sum(_grepl("LOCUS", n) for n in segNames) != 1:
        r_stop("Number of LOCUS should be 1.")
    seg_name = _re.sub("DEFINITION {1,}", "", importedData[mainSegments[segNames.index("DEFINITION")] - 1]) \
        if "DEFINITION" in segNames else None
    kF = segNames.index("FEATURES") if "FEATURES" in segNames else None
    if kF is None:
        r_stop("Error in importedData[mainSegments[\"FEATURES\"]:...]: NA/NaN argument")
    if kF == len(mainSegments) - 1:
        dataFeatures = rcompat.r_index(importedData, r_colon(mainSegments[kF], len(importedData) - 1))
    else:
        dataFeatures = rcompat.r_index(importedData, r_colon(mainSegments[kF], mainSegments[kF + 1] - 1))
    dataFeatures = ["NA" if l is None else l for l in dataFeatures]
    if len(dataFeatures) < 2:
        r_stop("No FEATURES in GBK file.")
    indented = [k + 1 for k, l in enumerate(dataFeatures) if _re.search("^ {6,}", l)]
    if not indented:  # c(1:n)[-integer(0)] is empty in R
        startLineOfFeature = []
    else:
        skip = set(indented)
        startLineOfFeature = [k for k in range(1, len(dataFeatures) + 1) if k not in skip]
    startLineOfFeature = startLineOfFeature + [len(dataFeatures) + 1]
    nF = len(startLineOfFeature) - 1
    if nF < 1:
        r_stop("Error in dataFeatures[...]: NA/NaN argument (no feature lines)")

    cols = {k: [] for k in ("name", "start", "end", "strand", "length", "pid", "gene", "synonym", "product",
                            "proteinid", "feature", "gene_type", "color")}
    excluded = []
    for counter in range(1, nF + 1):
        lines = dataFeatures[startLineOfFeature[counter - 1] - 1:startLineOfFeature[counter] - 1]
        currentFeature = [_re.sub(r'^ |:|"| $', "", _re.sub(r"[ \t]+|" + _SPACE + "+", " ", piece))
                          for piece in r_strsplit("".join(lines), "   /")]
        if not currentFeature:
            continue
        key = _re.sub(" " + _PRINT + "+", "", currentFeature[0])
        if not any(_grepl(key, t) for t in tagsToParse):
            continue
        tag = _re.sub(" " + _GRAPH + "+", "", currentFeature[0])
        exons = r_strsplit(_re.sub(_ALPHA + r"|_| |\(|\)|", "", currentFeature[0]), ",")
        complement = _grepl("complement", currentFeature[0])
        if complement:
            exonVector = [tag + " complement(" + e + ")" for e in exons]
        else:
            exonVector = [tag + " " + e for e in exons]
        if len(exonVector) > 1:
            exonVector2 = []
            for i in range(len(exonVector) - 1):
                a = _paste_num([v + 1 for v in get_end(exonVector[i])])
                b = _paste_num([v - 1 for v in get_start(exonVector[i + 1])])
                if complement:
                    intron = tag + "_intron complement(" + a + ".." + b + ")"
                else:
                    intron = tag + "_intron " + a + ".." + b
                exonVector2 += [exonVector[i], intron]
            exonVector = exonVector2 + [exonVector[-1]]
        for currentExon in exonVector:
            currentFeature[0] = currentExon
            has_gene = any(_grepl("gene=", el) for el in currentFeature)
            nameTEMP = extract_data("gene=", currentFeature) if has_gene else extract_data("locus_tag=", currentFeature)
            startTEMP = get_start(currentFeature[0])
            endTEMP = get_end(currentFeature[0])
            if not startTEMP or not endTEMP:
                excluded.append(nameTEMP)
                continue
            cols["name"].append(nameTEMP)
            cols["start"].append(startTEMP[0])
            cols["end"].append(endTEMP[0])
            cols["length"].append((endTEMP[0] - startTEMP[0] + 1) / 3 - 1)
            cols["strand"].append(-1.0 if _grepl("complement", currentFeature[0]) else 1.0)
            cols["pid"].append(extract_data("db_xref=GI", currentFeature))
            cols["gene"].append(extract_data("gene=", currentFeature) if has_gene else "-")
            cols["synonym"].append(extract_data("locus_tag=", currentFeature))
            cols["proteinid"].append(extract_data("protein_id=", currentFeature))
            cols["product"].append(extract_data("product=", currentFeature))
            cols["color"].append(extract_data("(color|colour)=", currentFeature))
            cols["gene_type"].append("introns" if _grepl("intron", currentFeature[0]) else gene_type)
            feat = _re.sub(" " + _PRINT + "+", "", currentFeature[0])
            cols["feature"].append(feat + "_pseudo" if any(_re.search("^pseudo", el) for el in currentFeature) else feat)
    for nm in excluded:  # R gives these as warnings at the end
        print(f"Warning message:\nStart and stop position invalid for {nm}", file=sys.stderr, flush=True)

    color = cols.pop("color")
    artemis_numbers = [str(k) for k in range(len(_ARTEMIS))]
    if color and all(c in ["NA"] + artemis_numbers for c in color):
        color = [c if c == "NA" else _ARTEMIS[int(c)] for c in color]
    gt = cols["gene_type"]
    if gene_type == "auto":
        fill = "exons" if any("intron" in g for g in gt) else "bars"
        cols["gene_type"] = [fill if g == "auto" else g for g in gt]
    table = pd.DataFrame({k: pd.Series(v, dtype=float if k in ("start", "end", "strand", "length") else object)
                          for k, v in cols.items()})
    if not all(c == "NA" for c in color):
        table["col"] = ["blue" if c == "NA" else c for c in color]
    if len(table) == 0:
        return None
    # as.dna_seg()
    if "col" not in table:
        table["col"] = "blue"
    table["fill"] = "blue"
    table["lty"] = 1.0
    table["lwd"] = 1.0
    table["pch"] = 8.0
    table["cex"] = 1.0
    table.index = range(1, len(table) + 1)
    table.attrs["seg_name"] = seg_name
    return table


def reorder_reference_gbk(ref_gb):
    """The reference's CDS features, duplicate names suffixed _1, _2, ... (in
    file order), sorted by start position (stable)."""
    ref = read_dna_seg_from_file(ref_gb)
    if ref is None:
        r_stop("Error in ref[which(ref$feature == \"CDS\"), ]: incorrect number of dimensions")
    ref = ref[ref["feature"] == "CDS"].copy()
    # For each name, if there is more than one entry, label as 'gene_1', 'gene_2'
    counts = {}
    for nm in ref["name"]:
        counts[nm] = counts.get(nm, 0) + 1
    names = list(ref["name"])
    for nm in [n for n in counts if counts[n] > 1]:
        w = [k for k, n in enumerate(names) if n == nm]
        for j, k in enumerate(w):
            names[k] = nm + "_" + str(j + 1)
    ref["name"] = names
    return ref.iloc[rcompat.r_order(ref["start"].to_numpy(dtype=float))]
