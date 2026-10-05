#!/usr/bin/env python3
"""kmercontigalignonly.py: align contigs to the reference genome using nucmer and
assign k-mers to genes or intergenic regions. Port of kmercontigalignonly.Rscript.

Positions are 1-based throughout, as in R; R's NULL (no gene) is None."""
import argparse
import os
import sys
import time

import rcompat
from rcompat import r_cat, r_colon, r_index, r_paste0, r_stop

###################################################################################################
## Functions and software paths
###################################################################################################


def create_contigalign_dir(dir, kmer_type, kmer_length):
    contigalign_dir = dir + "/" + r_paste0(kmer_type, "kmer", kmer_length, "_kmergenealign/")  # file.path
    if not os.path.isdir(contigalign_dir):
        rcompat.r_dir_create(contigalign_dir)
    return contigalign_dir


def get_alignment_pos_indiv(x):
    """Reference start/end and contig start/end from a show-aligns BEGIN line."""
    t = x.split(" ")
    return [rcompat.r_as_numeric(t[k - 1]) if k <= len(t) else None for k in (6, 8, 11, 13)]


def get_alignment_pos(x):
    """The 8 fields of a show-coords line; the last is the contig name."""
    x = [t for t in x.split(" ") if t != "" and t != "|"]
    if not x:
        r_stop("Error in x[length(x)] = ...: replacement has length zero")
    last = x[-1].split("\t")
    x[-1] = last[1] if len(last) > 1 else None
    return x


def get_alignment(x):
    """The bases of a show-aligns sequence line: everything after the last space."""
    k = x.rfind(" ")
    if k < 0:
        r_stop("Error in (wh.space[length(wh.space)] + 1):length(align.split): argument of length 0")
    return x[k + 1:]


def kmerCount(protein, kmer_length):
    k = int(kmer_length)
    if len(protein) >= k:
        return [protein[s:s + k] for s in range(len(protein) - (k - 1))]
    return None


def count_protein_kmers(proteins, kmer_length=31):
    """Count kmers for one translated contig (R passes a single string)."""
    return kmerCount(proteins, kmer_length)


def count_varlength_kmers(seq, kmerLengths=(9, 100)):
    k1, k2 = kmerLengths
    if k2 > len(seq):
        k2 = len(seq)
    if k1 > len(seq):
        return None
    out = []
    for k in r_colon(int(k1), int(k2)):
        km = kmerCount(seq, k)
        if km is not None:
            out += km
    return out if out else None


def get_varlength_pos_kmer(pos, kmer_length):
    k = int(kmer_length)
    return [r_index(pos, r_colon(x, x + k - 1)) for x in r_colon(1, len(pos) - k + 1)]


class PosGenes:
    """R's ref.pos.gene.id: for each reference position (1-based) the gene and
    intergenic IDs assigned to it, in order of assignment. Stored as the first ID
    in an array plus a dict of any further IDs (overlapping genes), to keep a
    bacterial genome's worth of positions small."""

    def __init__(self, n):
        import numpy as np
        self.n = n
        self.first = np.zeros(n + 1, dtype=np.int64)
        self.extra = {}

    def add_range(self, a, b, gid):
        import numpy as np
        lo, hi = (a, b) if a <= b else (b, a)
        if lo < 1:
            r_stop("Error in ref.pos.gene.id[[j]]: subscript out of bounds (position ", lo, ")")
        if hi > self.n:  # R extends a list assigned beyond its end
            self.first = np.concatenate([self.first, np.zeros(hi - self.n, dtype=np.int64)])
            self.n = hi
        seg = self.first[lo:hi + 1]
        taken = np.flatnonzero(seg != 0)
        seg[seg == 0] = gid
        for k in taken:
            self.extra.setdefault(lo + int(k), []).append(gid)

    def genes(self, j):
        """The IDs at position j, or None (R's NULL)."""
        if j < 1 or j > self.n or self.first[j] == 0:
            return None
        g = (int(self.first[j]),)
        return g + tuple(self.extra[j]) if j in self.extra else g

    def n_distinct(self):
        import numpy as np
        ids = set(np.unique(self.first[1:]).tolist()) - {0}
        for v in self.extra.values():
            ids.update(v)
        return len(ids)


class WindowGenes:
    """unique(unlist(vals[window])) for windows of consecutive positions, via runs
    of positions with the same IDs. vals is 1-based (vals[0] unused); None is NULL."""

    def __init__(self, vals):
        self.vals = vals
        self.run_of = [0] * len(vals)
        self.run_start, self.run_val = [], []
        prev = object()
        for p in range(1, len(vals)):
            v = vals[p]
            if v != prev:
                self.run_start.append(p)
                self.run_val.append(v)
                prev = v
            self.run_of[p] = len(self.run_start) - 1

    def window(self, w):
        """w: list of positions (None for R's NA)."""
        n = len(self.vals)
        if w and None not in w and w[-1] - w[0] == len(w) - 1 and 1 <= w[0] and w[-1] < n:
            out, seen = [], set()
            r = self.run_of[w[0]]
            while r < len(self.run_start) and self.run_start[r] <= w[-1]:
                v = self.run_val[r]
                if v is not None:
                    for g in v:
                        if g not in seen:
                            seen.add(g)
                            out.append(g)
                r += 1
            return out if out else None
        return unique_unlist(self.vals[p] if p is not None and 1 <= p < n else None for p in w)


def read_reference_name(ref_fa):
    """ref.name as read in read_reference_files: the first word of the header, without '>'."""
    with rcompat.r_open(ref_fa) as f:
        lines = f.read().split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    ref_name = lines[0] if lines and lines[0] != "" else None
    if sum(l.startswith(">") for l in lines) > 1:
        r_stop("Error: reference fasta file ", ref_fa, " contains more than one record; only single-record "
               "references are supported", "\n")
    if ref_name is None:
        r_stop("Error in substr(ref.name, 1, 1): argument is of length zero")
    ref_name = ref_name.split(" ")[0]
    if ref_name[:1] != ">":
        r_stop("Error: reference fasta file does not begin with a name starting with '>'", "\n")
    return ref_name[1:1000000]


def read_reference_files(ref_gb, ref_fa, ref_length, process, output_dir, prefix, kmer_type, kmer_length):
    import sequence_functions

    # Get the reference name
    ref_name = read_reference_name(ref_fa)
    r_cat("Reference name:", ref_name, "\n")
    # scan(ref_gb, nlines = 1) is not quiet
    with rcompat.r_open(ref_gb) as f:
        first = f.readline().rstrip("\n")
    print("Read 1 item" if first != "" else "Read 0 items", file=sys.stderr, flush=True)
    toks = [t for t in first.split(" ") if t != ""]
    ref_length = rcompat.r_as_numeric(toks[2]) if len(toks) >= 3 else None
    if ref_length is None:
        r_stop("Error retrieving the reference genome length from the genbank file", "\n")
    r_cat("Reference genome length:", ref_length, "\n")

    # Read in reference genbank file
    ref = sequence_functions.reorder_reference_gbk(ref_gb=ref_gb)
    names = list(ref["name"])
    starts = [int(v) for v in ref["start"]]
    ends = [int(v) for v in ref["end"]]
    nref = len(ref)

    # Get a gene/IR ID for every position in the reference
    # A position can have multiple gene IDs as there are overlapping genes
    ref_pos_gene_id = PosGenes(int(ref_length))
    for i in range(1, nref + 1):
        ref_pos_gene_id.add_range(starts[i - 1], ends[i - 1], i)
    # Intergenic ID only assigned if the start of one gene is after the end of the previous gene
    intergenic_names = []
    for i in r_colon(2, nref):
        if i > nref:
            r_stop("Error in if (inter_start <= inter_end): missing value where TRUE/FALSE needed")
        inter_start = ends[i - 2] + 1
        inter_end = starts[i - 1] - 1
        if inter_start <= inter_end:
            intergenic_names.append(names[i - 2] + ":" + names[i - 1])
            ref_pos_gene_id.add_range(inter_start, inter_end, nref + len(intergenic_names))
        if i == nref:
            intergenic_names.append(names[-1] + ":")
            # The final intergenic region runs from after the last gene to the end of the
            # (circular) reference and wraps round to the base before the first gene
            gid = nref + len(intergenic_names)
            if ends[nref - 1] + 1 <= ref_pos_gene_id.n:
                ref_pos_gene_id.add_range(ends[nref - 1] + 1, ref_pos_gene_id.n, gid)
            if starts[0] - 1 >= 1:
                ref_pos_gene_id.add_range(1, starts[0] - 1, gid)

    r_cat("Read in reference and assigned a gene/intergenic region ID to every position", "\n")
    if process == 1:
        # Create a lookup table to get a gene name/IR name for every ID assigned
        lookup = names + intergenic_names
        n_ids = ref_pos_gene_id.n_distinct()
        if n_ids > len(lookup):
            r_stop("Error in names(gene_id_lookup) = ...: 'names' attribute must be the same length as the vector")
        ids = [str(k) for k in range(1, n_ids + 1)] + ["NA"] * (len(lookup) - n_ids)
        with open(r_paste0(output_dir, prefix, "_", kmer_type, kmer_length, "_", ref_name, "_gene_id_name_lookup.txt"),
                  "w") as f:
            f.write("".join(nm + "\t" + k + "\n" for nm, k in zip(lookup, ids)))
        r_cat("Written gene/IR ID lookup to file", "\n")

    return {"ref": ref, "ref.pos.gene.id": ref_pos_gene_id, "ref_length": ref_length, "ref.name": ref_name}


def read_ids(id_file):
    ids = rcompat.r_read_table(id_file, header=True, sep="\t", comment_char="")
    ids.columns = [c.lower() for c in ids.columns]
    if "id" not in ids.columns or "paths" not in ids.columns:
        r_stop("Column names for sample ID and assembly path file must be 'id' and 'paths' in any case", "\n")
    assemblies = [rcompat.r_as_character(v) for v in ids.iloc[:, list(ids.columns).index("paths")]]
    idv = [rcompat.r_as_character(v) for v in ids.iloc[:, list(ids.columns).index("id")]]
    if not all(os.path.exists(a) for a in assemblies):
        r_stop("Error: not all assembly files exist", "\n")
    return {"ids": idv, "assemblies": assemblies}


def read_sorted_kmers(kmerSeqFile):
    # Read in full list of sorted kmers
    return rcompat.r_scan_lines(kmerSeqFile, quiet=True)


def write_kmergenecomb_filepaths_to_file(ids, kmer_type, kmer_length, prefix, output_prefix, ident_threshold, ref_name,
                                         output_dir, final_file_prefix):
    # For the last process, create file containing paths to all output files
    # Won't check if they are all created - flag warning
    all_outfiles = [final_file_prefix + i + "_nucmeralign_kmer_list_gene_IDs.txt.gz" for i in ids]
    all_outfiles_path = r_paste0(output_dir, prefix, "_", kmer_type, kmer_length, "_", ref_name,
                                 "_kmergenecombination_filepaths.txt")
    r_cat(r_paste0("Writing file paths to all kmer gene combinations per sample for ", kmer_type, " kmer length ",
                   kmer_length),
          "to file (Warning: have not checked that all kmer contig alignment is completed and all files exist):",
          all_outfiles_path, "\n")
    rcompat.r_cat_lines(all_outfiles, all_outfiles_path)
    return all_outfiles_path


# Allowed characters in the alignments from nucmer - check all alignments against this to make sure no other characters are present
allowed_chars = set("acgt.n")


def unique_unlist(vectors):
    """unique(unlist(list_of_vectors)); None (R's NULL) when there is nothing."""
    out, seen = [], set()
    for v in vectors:
        if v is None:
            continue
        for g in v:
            if g not in seen:
                seen.add(g)
                out.append(g)
    return out if out else None


def gene_index_for_kmers(windows, out):
    """sapply(windows, function(x) list(unique(unlist(out[x])))): the genes of the
    positions of each k-mer (positions None are R's NA, giving NULL)."""
    wg = WindowGenes(out)
    return [wg.window(w) for w in windows]


def covered_contig_bases(alignfile, passing):
    """Contig positions covered by the passing alignments: each alignment's contig start to end
    (columns 3-4 of its show-aligns header). D1c: R's matrix(..., byrow = T) paired the starts
    with each other and the ends with each other when more than one alignment passed, so the
    unaligned-base count in the log was wrong."""
    covered = set()
    for j in passing:
        covered.update(r_colon(int(alignfile[j][2]), int(alignfile[j][3])))
    return covered


def kmer_gene_pairs(aggregate):
    """One (gene ID, k-mer) pair per gene of each k-mer, from {k-mer: unique gene IDs}, in the
    dict's order. D1d: R's sapply() made a matrix without names when every k-mer had the same
    number (> 1) of genes, so the genome's k-mer/gene pairs were written without their k-mers."""
    genes = [float(g) for v in aggregate.values() for g in v]
    names = [km for km, v in aggregate.items() for _ in v]
    return genes, names


def align_contig(c, contig, contig_id_c, kmer_type, kmer_length, kstart, kend, ident_threshold, alignment_pos,
                 delta, ref_name, mummer_path, ref_pos_gene_id, ids_i, oneLetterCodes, revcompl):
    """One pass of the R loop over contigs: (kmers, gene index of each k-mer, bases not aligned)."""
    import sequence_functions
    L = len(contig)
    if kmer_type == "protein":
        # Translate contig c
        contig_translate_c = [sequence_functions.translate_function(contig, frame, oneLetterCodes, revcompl)
                              for frame in range(1, 7)]
        # Get the kmers for each of the translated contigs
        kmers_contig_c = [count_protein_kmers(x, kmer_length) for x in contig_translate_c]
        r_cat("Translated contig", c, "\n")
        kmers_flat = [k for v in kmers_contig_c if v is not None for k in v]
        kmers_null = False
    else:
        # Pull out contig c and set to uppercase
        contigs_c = [contig.upper()]
        # Remove Ns, creating multiple 'contigs' from the contig if Ns present
        if "N" in contigs_c[0]:
            contigs_c = sequence_functions.remove_Ns(contigs_c[0])
        kmers_flat = []
        for x in contigs_c:
            km = kmerCount(x, kmer_length) if kmer_length > 0 else count_varlength_kmers(x, (kstart, kend))
            if km is not None:
                kmers_flat += km
        kmers_null = not kmers_flat

    # Run for every contig
    contig_alignment_c = rcompat.r_system_intern(mummer_path + "show-aligns " + delta + " " + ref_name + " " + contig_id_c)

    no_match = None
    if len(contig_alignment_c) > 0 and not kmers_null:
        # Find the lines in the alignment output containing the start and end of each alignment
        wh_begin = [k + 1 for k, l in enumerate(contig_alignment_c) if "BEGIN alignment" in l]
        wh_end = [k + 1 for k, l in enumerate(contig_alignment_c) if "END alignment" in l]
        if len(wh_begin) != len(wh_end):
            r_stop("Error: number of BEGIN alignment matches is not equal to the number of END alignment matches for ID",
                   ids_i, "and contig", c, contig_id_c, "\n")
        # Pull out from the alignment the start and end positions for contig c
        alignfile = [get_alignment_pos_indiv(contig_alignment_c[k - 1]) for k in wh_begin]
        # Pull out from the full alignment positions matrix the rows for contig c
        coordsfile = [row for row in alignment_pos if row[7] == contig_id_c]
        for j in range(len(alignfile)):
            if j >= len(coordsfile):
                r_stop("Error in alignment.pos.coordsfile[j, ]: subscript out of bounds")
            if any(a != rcompat.r_as_numeric(b) for a, b in zip(alignfile[j], coordsfile[j][:4])):
                r_stop("Error: alignment coordinates do not match between the alignment file and the coordinates file", "\n")
        passing = [j for j in range(len(coordsfile)) if rcompat.r_as_numeric(coordsfile[j][6]) >= ident_threshold]
        # Get all bases which are covered by an alignment for contig c that pass the identity threshold
        covered = covered_contig_bases(alignfile, passing)
        # Write the number of bases that are not part of any alignment for contig c
        no_match = sum(1 for p in range(1, L + 1) if p not in covered)

        contig_pos_all = []
        ref_pos_all = []
        # For every alignment for contig c
        for j in passing:
            if j >= len(wh_begin):
                r_stop("Error in (wh.begin[j] + 1):(wh.end[j] - 1): NA/NaN argument")
            # First pull out the lines of the alignment j
            lines = r_index(contig_alignment_c, r_colon(wh_begin[j] + 1, wh_end[j] - 1))
            lines = ["NA" if l is None else l for l in lines]
            # Remove all empty lines; keep those that don't start with a space (the alignments)
            lines = [l for l in lines if l != "" and l[0] != " "]
            lines = [get_alignment(l) for l in lines]
            # Concatenate the reference and query alignments, which alternate in lines
            n_half = -(-len(lines) // 2)  # seq(length.out = n/2) rounds up
            ref_align = "".join("NA" if v is None else v for v in r_index(lines, [1 + 2 * k for k in range(n_half)]))
            query_align = "".join("NA" if v is None else v for v in r_index(lines, [2 + 2 * k for k in range(n_half)]))
            # Check that there are no unknown characters in either alignment
            if not set(ref_align) <= allowed_chars:
                r_stop("Unknown characters in reference sequence", "\n")
            if not set(query_align) <= allowed_chars:
                r_stop("Unknown characters in reference sequence", "\n")
            # Get the length of the reference and query in the alignment that are not gaps
            ref_nchar = len(ref_align) - ref_align.count(".")
            query_nchar = len(query_align) - query_align.count(".")
            if len(ref_align) != len(query_align):
                r_stop("Error: reference alignment length not equal to query alignment length for ID", ids_i,
                       "and contig", c, contig_id_c, "j", j + 1, "\n")
            if str(ref_nchar) != coordsfile[j][4]:
                r_stop("Error: reference alignment extracted is not the correct length for ID", ids_i, "and contig", c,
                       contig_id_c, "j", j + 1, "\n")
            if str(query_nchar) != coordsfile[j][5]:
                r_stop("Error: reference alignment extracted is not the correct length for ID", ids_i, "and contig", c,
                       contig_id_c, "j", j + 1, "\n")
            # For every position in the reference alignment that is not a gap, assign it its position in the
            # reference; for every position that is a gap, the highest leftmost position that is not a gap
            positions = r_colon(int(alignfile[j][0]), int(alignfile[j][1]))
            ref_pos = []
            it = iter(positions)
            last = None
            for ch in ref_align:
                if ch != ".":
                    last = next(it)
                    ref_pos.append(last)
                else:
                    if last is None:
                        r_stop("Error in ref.pos.new[...] = ...: replacement has length zero")
                    ref_pos.append(last)
            contig_pos_all += r_colon(int(alignfile[j][2]), int(alignfile[j][3]))
            ref_pos_all += [rp for rp, q in zip(ref_pos, query_align) if q != "."]

        if contig_pos_all:
            # aggregate(ref_pos_all, by = list(contig_pos_all), FUN = "unique"): for each contig position,
            # the reference positions aligned to it (first appearance order)
            groups = {}
            for cp, rp in zip(contig_pos_all, ref_pos_all):
                g = groups.setdefault(cp, [])
                if rp not in g:
                    g.append(rp)
            # out[p]: the gene IDs (with repeats) of the reference positions aligned to contig position p
            out = [None] * (L + 1)
            for cp, rps in groups.items():
                if 1 <= cp <= L:
                    genes = tuple(g for rp in rps for g in (ref_pos_gene_id.genes(rp) or ()))
                    out[cp] = genes if genes else None

            if kmer_type == "protein":
                # Get a reference gene for every amino acid in the translated contigs (the middle base)
                frames = [(2, 3), (3, 3), (4, 3), (L - 1, -3), (L - 2, -3), (L - 3, -3)]
                kmers_genes_c = []
                for j in range(6):
                    start, by = frames[j]
                    n_aa = len(contig_translate_c[j])
                    idx = [start + by * k for k in range(n_aa)]
                    pos_frame = [None] + [out[p] if 1 <= p <= L else None for p in idx]  # 1-based
                    # positions of the amino acids that are not X (all of them now)
                    wh_notX = [k + 1 for k, a in enumerate(contig_translate_c[j]) if a != "X"]
                    nk = len(kmers_contig_c[j]) if kmers_contig_c[j] is not None else 0
                    k = int(kmer_length)
                    windows = [r_index(wh_notX, [m + x - 1 for m in range(1, k + 1)]) for x in r_colon(1, nk)]
                    kmers_genes_c += gene_index_for_kmers(windows, pos_frame)
            else:
                # Which positions in the contig are not Ns; split where the Ns were
                wh_notN = [k + 1 for k, b in enumerate(contig) if b != "N"]
                which_break = [k + 1 for k in range(len(wh_notN) - 1) if wh_notN[k + 1] - wh_notN[k] > 1]
                if which_break:
                    breaks = []
                    nb = len(which_break)
                    for kk in range(1, nb + 1):
                        if kk == 1:
                            breaks.append(wh_notN[0:which_break[0]])
                            if nb == 1:
                                breaks.append(wh_notN[which_break[0]:])
                        elif kk != nb:
                            breaks.append(wh_notN[which_break[kk - 2]:which_break[kk - 1]])
                        else:
                            breaks.append(wh_notN[which_break[kk - 2]:which_break[kk - 1]])
                            breaks.append(wh_notN[which_break[kk - 1]:])
                else:
                    breaks = [wh_notN]
                windows = []
                if kmer_length > 0:
                    k = int(kmer_length)
                    for seg in breaks:
                        xs = r_colon(1, len(seg) - k + 1)
                        if any(x < 0 for x in xs):
                            r_stop("Error in wh.notN.breaks[c(1:kmer_length) + x - 1]: can't mix positive and "
                                   "negative subscripts")
                        windows += [r_index(seg, [m + x - 1 for m in range(1, k + 1)]) for x in xs]
                else:
                    for seg in breaks:
                        kend_j = min(kend, len(seg))
                        for kk in r_colon(int(kstart), int(kend_j)):
                            windows += get_varlength_pos_kmer(seg, kk)
                kmers_genes_c = gene_index_for_kmers(windows, out)
        else:
            # If none of the alignments were above the identity threshold, set all kmers to NULL
            kmers_genes_c = [None] * len(kmers_flat)
    else:
        # If there were no alignments at all for contig c, set length of no match to the length of the contig
        no_match = L
        kmers_genes_c = [None] * len(kmers_flat)
    return kmers_flat, kmers_genes_c, no_match


###################################################################################################


def parse_args(prog_description):
    parser = argparse.ArgumentParser(description=prog_description, allow_abbrev=False)
    parser.add_argument("--task-id", required=True)
    parser.add_argument("--n", required=True, help="number of samples")
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--output-dir", required=True, help="analysis directory")
    parser.add_argument("--id-file", required=True)
    parser.add_argument("--ref-fa", required=True, help="reference genome (FASTA)")
    parser.add_argument("--ref-gb", required=True, help="reference annotation (GenBank)")
    parser.add_argument("--kmer-type", required=True, help="protein or nucleotide")
    parser.add_argument("--kmer-length", required=True, help="k-mer length (0: variable length nucleotide k-mers)")
    parser.add_argument("--nucmerident", required=True, help="minimum nucmer alignment identity (%%)")
    parser.add_argument("--kmerseqfile", required=True, help="merged list of all k-mers")
    parser.add_argument("--software-file", required=True)
    parser.add_argument("--kstart", default=None, help="shortest variable k-mer length (default 9)")
    parser.add_argument("--kend", default=None, help="longest variable k-mer length (default 100)")
    return parser.parse_args()


def run(script_path, description, merge_hook=None):
    """The body of kmercontigalignonly; kmercontigalign.py runs it with a
    merge_hook that launches kmercontigalignmerge.py."""
    rcompat.script_setup(script_path)
    start_time = time.monotonic()
    args = parse_args(description)

    # Initialise variables
    process = rcompat.r_as_integer(args.task_id)
    n = rcompat.r_as_integer(args.n)
    prefix = args.output_prefix
    output_dir = args.output_dir
    id_file = args.id_file
    ref_fa = args.ref_fa
    ref_gb = args.ref_gb
    kmer_type = args.kmer_type.lower()
    kmer_length = rcompat.r_as_numeric(args.kmer_length)
    ident_threshold = rcompat.r_as_numeric(args.nucmerident)
    kmerSeqFile = args.kmerseqfile
    software_file = args.software_file
    kparams_given = args.kstart is not None or args.kend is not None
    if kparams_given:
        kstart = rcompat.r_as_integer(args.kstart) if args.kstart is not None else None
        kend = rcompat.r_as_integer(args.kend) if args.kend is not None else None
    else:
        kstart = 9.0
        kend = 100.0
        if kmer_length == 0:
            r_cat("Assuming variable kmer lengths between 9-100 bases long", "\n")

    # Check file inputs
    if n is None:
        r_stop("Error: n must be an integer", "\n")
    if not os.path.exists(output_dir):
        r_stop("Error: output directory doesn't exist", "\n")
    if not output_dir.endswith("/"):
        output_dir = output_dir + "/"
    if not os.path.exists(id_file):
        r_stop("Error: sample ID file doesn't exist", "\n")
    if not os.path.exists(ref_fa):
        r_stop("Error: reference fasta file doesn't exist", "\n")
    if not os.path.exists(ref_gb):
        r_stop("Error: reference genbank file doesn't exist", "\n")
    if kmer_type != "protein" and kmer_type != "nucleotide":
        r_stop("Error: kmer type must be either protein or nucleotide", "\n")
    if kmer_length is None:
        r_stop("Error: kmer length must be an integer", "\n")
    if ident_threshold is None or ident_threshold > 100 or ident_threshold < 0:
        r_stop("Error: nucmer identity threshold must be between 0-100", "\n")
    if not os.path.exists(kmerSeqFile):
        r_stop("Error: kmer sequence file doesn't exist", "\n")
    if not os.path.exists(software_file):
        r_stop("Error: software file doesn't exist", "\n")
    # R checks kstart/kend whenever there are more than 10 arguments, i.e. always
    if kstart is None or kend is None:
        r_stop("Error: kstart and kend must be integers", "\n")

    # Read in software file
    software_paths = rcompat.r_read_table(software_file, header=True, sep="\t", quote="")
    names = [rcompat.r_as_character(v) for v in software_paths["name"]]
    paths = [rcompat.r_as_character(v) for v in software_paths["path"]]

    def software(name):
        return [pth for nm, pth in zip(names, paths) if nm.lower() == name][0]
    # Required software and script paths (R and genoPlotR are no longer used, but software
    # files keep their entries)
    required_software = ["scriptpath", "R", "mummer", "genoPlotR"]
    if any(r not in names for r in required_software):
        r_stop("Error: missing required software path in the software file - requires " + ", ".join(required_software), "\n")
    python_path = sys.executable
    script_location = software("scriptpath")
    if not os.path.isdir(script_location):
        r_stop("Error: script location directory specified in the software paths file doesn't exist", "\n")
    kmercontigalignmergepath = script_location + "/kmercontigalignmerge.py"
    if not os.path.exists(kmercontigalignmergepath):
        r_stop("Error: kmercontigalignmerge.py path doesn't exist - check pipeline script location in the software file", "\n")
    sequence_functions_file = script_location + "/sequence_functions.py"
    if not os.path.exists(sequence_functions_file):
        r_stop("Error: sequence_functions.py path doesn't exist - check pipeline script location in the software file", "\n")
    sys.path.insert(0, script_location)
    import sequence_functions

    mummer_path = software("mummer")
    if not os.path.exists(mummer_path):
        r_stop("Error: mummer software directory path doesn't exist", "\n")
    if not mummer_path.endswith("/"):
        mummer_path = mummer_path + "/"

    sortstringspath = script_location + "/sort_strings"
    if not os.path.exists(sortstringspath):
        r_stop("Error: sort_strings path doesn't exist - check pipeline script location in the software file", "\n")

    # Report variables
    r_cat("#############################################", "\n")
    r_cat("Running on host: ", rcompat.r_system_intern("hostname"), "\n")
    r_cat("Command line arguments", "\n")
    r_cat(sys.argv[1:], "\n\n")
    r_cat("Parameters:", "\n")
    r_cat("task_id:", process, "\n")
    r_cat("n:", n, "\n")
    r_cat("Output prefix:", prefix, "\n")
    r_cat("Analysis directory:", output_dir, "\n")
    r_cat("ID file path:", id_file, "\n")
    r_cat("Reference fasta file:", ref_fa, "\n")
    r_cat("Reference genbank file:", ref_gb, "\n")
    r_cat("Kmer type:", kmer_type, "\n")
    r_cat("Kmer length:", kmer_length, "\n")
    r_cat("Nucmer alignment minimum % identity:", ident_threshold, "\n")
    r_cat("Kmer list file:", kmerSeqFile, "\n")
    r_cat("Software file:", software_file, "\n")
    r_cat("Script location:", script_location, "\n")
    r_cat("Python path:", python_path, "\n")
    r_cat("mummer path:", mummer_path, "\n")
    r_cat("Kmer start length:", kstart, "\n")
    r_cat("Kmer end length:", kend, "\n")
    r_cat("#############################################", "\n\n")

    # Create an output directory
    contigalign_dir = create_contigalign_dir(dir=output_dir, kmer_type=kmer_type, kmer_length=kmer_length)

    # Read in sample IDs
    idinfo = read_ids(id_file=id_file)
    assemblies = idinfo["assemblies"]
    ids = idinfo["ids"]

    # Read in reference genbank file
    refinfo = read_reference_files(ref_gb=ref_gb, ref_fa=ref_fa, ref_length=None, process=process,
                                   output_dir=contigalign_dir, prefix=prefix, kmer_type=kmer_type,
                                   kmer_length=kmer_length)
    ref_pos_gene_id = refinfo["ref.pos.gene.id"]
    ref_name = refinfo["ref.name"]
    # Read in full list of sorted kmers
    final_kmer_list = read_sorted_kmers(kmerSeqFile=kmerSeqFile)
    first_index = {}
    for k, km in enumerate(final_kmer_list):
        first_index.setdefault(km, k + 1)

    final_file_prefix = None
    for i in [process]:
        if not 1 <= i <= len(ids):
            r_stop("Error in assemblies[i]: task_id ", i, " is not a row of the ID file")
        # Find the contig file for sample ID i
        contig_file = assemblies[i - 1]
        id_i = ids[i - 1]
        r_cat("\n")
        r_cat("Running nucmer", "\n")
        r_cat("\n")

        # Changed to remove Ns from contig before alignment
        contig_seq, contig_names_final = sequence_functions.read_contigs_removeNs(contig_file)
        temp_fa = r_paste0(contigalign_dir, id_i, "_", kmer_type, kmer_length, "_contigs_temp.fa")
        n_rec = max(len(contig_names_final), len(contig_seq)) if contig_names_final and contig_seq else 0
        rcompat.r_cat_lines([contig_names_final[k % len(contig_names_final)] + "\n" + contig_seq[k % len(contig_seq)]
                             for k in range(n_rec)], temp_fa)

        query = r_paste0(contigalign_dir, id_i, "_", kmer_type, kmer_length, "_query")
        rcompat.r_system(mummer_path + "nucmer --prefix=" + query + " " + ref_fa + " " + temp_fa)

        r_cat("\n")

        contig_id = [x[1:10000000] for x in contig_names_final]
        contigs = contig_seq
        contig_length = [len(x) for x in contigs]

        # Get positions for all alignments for contig file i
        coords_lines = rcompat.r_system_intern(mummer_path + "show-coords " + query + ".delta")
        coords_lines = r_index(coords_lines, r_colon(6, len(coords_lines)))
        if any(l is None for l in coords_lines):
            r_stop("Error in get_alignment_pos: no alignments in ", query, ".delta")
        alignment_pos = [get_alignment_pos(x) for x in coords_lines]
        if any(len(row) != 8 for row in alignment_pos):
            r_stop("Error in colnames(alignment.pos) = ...: show-coords lines do not have 8 fields")
        r_cat("Got all alignment positions", "\n")

        all_kmers_genome_i = []
        all_kmers_gene_index_i = []
        length_contig_no_match = []

        for c in r_colon(1, len(contig_id)):
            if c > len(contig_id) or c < 1:
                r_stop("Error in contigs[c]: no contigs in ", contig_file)
            kmers_c, genes_c, no_match = align_contig(
                c, contigs[c - 1], contig_id[c - 1], kmer_type, kmer_length, kstart, kend, ident_threshold,
                alignment_pos, query + ".delta", ref_name, mummer_path, ref_pos_gene_id, id_i,
                sequence_functions.oneLetterCodes, sequence_functions.revcompl)
            length_contig_no_match.append(no_match)
            # Add the kmers and gene pos for the kmers from contig c to the total
            all_kmers_genome_i += kmers_c
            all_kmers_gene_index_i += genes_c
            r_cat("Finished for contig", c, "of", len(contig_id), time.strftime("%Y-%m-%d %H:%M:%S"), "\n")
        if len(all_kmers_genome_i) != len(all_kmers_gene_index_i):
            r_stop("Length all_kmers_genome_i != length all_kmers_gene_index_i", "\n")
        # Remove all that do not have a gene
        keep = [k for k, g in enumerate(all_kmers_gene_index_i) if g is not None]
        all_kmers_genome_i = [all_kmers_genome_i[k] for k in keep]
        all_kmers_gene_index_i = [all_kmers_gene_index_i[k] for k in keep]

        # Get list of unique kmers in genome i
        # For every unique kmer sequence, get all gene IDs. (R orders the groups by
        # collating the k-mers; the order does not reach the output, which is sorted.)
        aggregate = {}
        for km, genes in zip(all_kmers_genome_i, all_kmers_gene_index_i):
            aggregate.setdefault(km, []).append(genes)
        if not aggregate:
            r_stop("Error in agg_kmers[[index]]: subscript out of bounds (no k-mers aligned for ID ", id_i, ")")
        aggregate = {km: unique_unlist(v) for km, v in sorted(aggregate.items())}
        unique_kmers_genes, unique_kmers_genes_names = kmer_gene_pairs(aggregate)
        ## Add in for reverse complement kmers
        if kmer_type == "nucleotide":
            rb = sequence_functions.rev_base
            unique_kmers_genes_names_revcompl = ["".join(rb[ch] for ch in reversed(x) if ch in rb)
                                                 for x in unique_kmers_genes_names]

        r_cat("Got list of all unique gene IDs for each kmer in genome", i, "\n")

        # Match to full list of sorted kmers
        unique_kmers_index_match = [first_index.get(x) for x in unique_kmers_genes_names]
        if kmer_type == "nucleotide":
            for k, x in enumerate(unique_kmers_genes_names_revcompl):
                if unique_kmers_index_match[k] is None:
                    unique_kmers_index_match[k] = first_index.get(x)

        # For any kmers not in the dsk output, just remove them with a warning
        if any(m is None for m in unique_kmers_index_match):
            r_cat("Warning: removing", sum(m is None for m in unique_kmers_index_match),
                  "kmers not in the provided kmer list", "\n")
            unique_kmers_genes = [g for g, m in zip(unique_kmers_genes, unique_kmers_index_match) if m is not None]
            unique_kmers_index_match = [m for m in unique_kmers_index_match if m is not None]

        # Write outfile - concatenate the index and the genes into a string for matching later
        if unique_kmers_index_match or not unique_kmers_genes:
            out = [str(m) + "," + rcompat.r_as_character(g) for m, g in zip(unique_kmers_index_match, unique_kmers_genes)]
        else:  # paste() recycles a zero-length vector as ""
            out = ["," + rcompat.r_as_character(g) for g in unique_kmers_genes]
        out = rcompat.r_unique(out)

        outfile = r_paste0(contigalign_dir, prefix, "_nucmeralign_kmer_list_gene_IDs_", id_i, "_t", ident_threshold,
                           "_unsorted.txt")
        with open(outfile, "w") as f:
            f.write("".join(o + "\t1\n" for o in out))
        # Sort the file
        final_file_prefix = r_paste0(contigalign_dir, prefix, "_", kmer_type, kmer_length, "_", ref_name, "_t",
                                     ident_threshold, "_")
        final_kmer_txt_gz = final_file_prefix + id_i + "_nucmeralign_kmer_list_gene_IDs.txt.gz"
        sortCommand = " ".join([sortstringspath, outfile, "| gzip -c >", final_kmer_txt_gz])
        r_cat("Sort command:", "\n")
        r_cat(sortCommand, "\n")
        r_cat("\n")
        rcompat.r_system(sortCommand)
        # Remove the temporary files
        rcompat.r_system("rm " + outfile)
        rcompat.r_system("rm " + query + ".delta")
        rcompat.r_system("rm " + temp_fa)
        if any(v is None for v in length_contig_no_match):
            unaligned = None
        else:
            unaligned = sum(length_contig_no_match)
        total = sum(contig_length)
        r_cat("Number of unaligned bases for ID", id_i, "NA" if unaligned is None else unaligned, "of", total,
              "NA" if unaligned is None else (unaligned / total) * 100, "%", "\n")

        # For the last process, create file containing paths to all kmer/gene combinations
        if i == len(ids):
            write_kmergenecomb_filepaths_to_file(ids=ids, kmer_type=kmer_type, kmer_length=kmer_length, prefix=prefix,
                                                 output_prefix=prefix, ident_threshold=ident_threshold,
                                                 ref_name=ref_name, output_dir=contigalign_dir,
                                                 final_file_prefix=final_file_prefix)

        # Create completed file
        outfile_completed = r_paste0(contigalign_dir, prefix, "_", kmer_type, kmer_length, "_", id_i,
                                     ".kmercontigalign.completed.txt")
        rcompat.r_system2("/bin/bash", "-c 'touch " + outfile_completed + "'")

        r_cat("Finished aligning contigs to the reference genome for process", i, "\n")

    state = {"process": process, "n": n, "prefix": prefix, "output_dir": output_dir, "contigalign_dir": contigalign_dir,
             "kmer_type": kmer_type, "kmer_length": kmer_length, "ref_name": ref_name, "ref_fa": ref_fa,
             "ident_threshold": ident_threshold, "ids": ids, "software_file": software_file,
             "python_path": python_path, "kmercontigalignmergepath": kmercontigalignmergepath,
             "start_time": start_time}
    if merge_hook is None:
        r_cat("Finished in", (time.monotonic() - start_time) / 60, "minutes\n")
    else:
        merge_hook(state)


def main():
    run(__file__, "kmercontigalignonly.py align contigs to the reference genome using nucmer and assign kmers to "
                  "genes or intergenic regions")


if __name__ == "__main__":
    main()
