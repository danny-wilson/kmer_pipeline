"""Sequence functions: reference reading, translation and contig reading.
Port of sequence_functions.R.

R works with character vectors of single bases; here a sequence is a str and a
character vector is a list of str. reorder_reference_gbk and create_gene_lookup
(which need the GenBank parser) are ported with the alignment step. transcribe()
is not called anywhere and is not ported.
"""
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
