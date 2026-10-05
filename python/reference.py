"""reference.py: references with several records (D6), such as a chromosome and its plasmids or
a draft assembly.

The records are laid end to end in file order into one coordinate system: the offset of a record
is the total length of the records before it. Every position in the pipeline is a position in
this system, so single-record references are handled exactly as before; the scripts take the
multi-record branches here only when there is more than one record.

- Records are named by their GenBank LOCUS names. FASTA and GenBank records are matched by order: the same number of records with the same
  lengths (an error otherwise); a name that differs is only a warning, as NCBI LOCUS names lack
  the version and Prokka shortens them.
- Gene names are made unique: duplicates within a record get _1, _2, ... as before; a name
  used in more than one record becomes name@record (":" already joins the genes either side
  of an intergenic region).
- Intergenic regions stay within a record: "a:b" between consecutive genes, and "z:" from the
  end of the record's last gene z round to the base before its first gene, as for the
  single-record reference (every record is treated as circular, whatever its LOCUS line says).
"""
import os

import rcompat
from rcompat import r_stop


class Record:
    def __init__(self, name, length, offset, locus=None, lines=None):
        self.name, self.length, self.offset, self.locus, self.lines = name, length, offset, locus, lines

    @property
    def start(self):
        return self.offset + 1

    @property
    def end(self):
        return self.offset + self.length


def fasta_records(ref_fa):
    """[(name, length)]: name is the first word of the header, without ">"."""
    out = []
    with rcompat.r_open(ref_fa) as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if line.startswith(">"):
                out.append([line[1:].split(" ")[0], 0])
            elif out:
                out[-1][1] += len(line.strip())
    return [(n, l) for n, l in out]


def fasta_sequence(ref_fa):
    """The sequences of all records, concatenated (positions as in the pipeline)."""
    with rcompat.r_open(ref_fa) as fh:
        return "".join(line.strip() for line in fh if not line.startswith(">"))


def genbank_records(ref_gb):
    """[(LOCUS name, length, the record's lines)] from the LOCUS lines (third word: the length)."""
    with rcompat.r_open(ref_gb) as fh:
        lines = fh.read().split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    starts = [k for k, line in enumerate(lines) if line.startswith("LOCUS")]
    out = []
    for i, k in enumerate(starts):
        block = lines[k:starts[i + 1] if i + 1 < len(starts) else len(lines)]
        toks = [t for t in block[0].split(" ") if t != ""]
        length = rcompat.r_as_numeric(toks[2]) if len(toks) >= 3 else None
        if length is None:
            r_stop("Error retrieving the length of GenBank record ", i + 1, " from its LOCUS line: ", block[0], "\n")
        out.append((toks[1] if len(toks) > 1 else "", int(length), block))
    return out


def check(ref_fa, ref_gb):
    """(errors, warnings) about the FASTA and GenBank files describing the same records."""
    errors, warnings = [], []
    fa = fasta_records(ref_fa)
    gb = genbank_records(ref_gb)
    if not fa:
        errors.append(f"no FASTA record in {ref_fa}")
    if not gb:
        errors.append(f"no LOCUS record in {ref_gb}")
    if errors:
        return errors, warnings
    if len(fa) != len(gb):
        errors.append(f"{ref_fa} has {len(fa)} records but {ref_gb} has {len(gb)}: they must describe the same "
                      "records, in the same order")
        return errors, warnings
    for k, ((fname, flen), (gname, glen, block)) in enumerate(zip(fa, gb), start=1):
        if flen != glen:
            errors.append(f"record {k}: {fname} is {flen} bases in {ref_fa} but {gname} is {glen} in {ref_gb}")
        names = {gname} | {line.split()[1] for line in block[:40]
                           if line.startswith(("VERSION", "ACCESSION")) and len(line.split()) > 1}
        if fname not in names and fname.split(".")[0] not in names:
            warnings.append(f"record {k} is called {fname} in {ref_fa} but {gname} in {ref_gb} (matched by order "
                            "and length)")
    return errors, warnings


def records(ref_gb):
    """The records, named by their LOCUS names, with their offsets."""
    out, offset = [], 0
    for name, length, block in genbank_records(ref_gb):
        out.append(Record(name, length, offset, name, block))
        offset += length
    return out


def n_records(ref_gb):
    with rcompat.r_open(ref_gb) as fh:
        return sum(1 for line in fh if line.startswith("LOCUS"))


def total_length(ref_gb):
    return sum(length for _, length, _ in genbank_records(ref_gb))


def record_of(recs, pos):
    """The record holding (global) position pos, or None."""
    for r in recs:
        if r.start <= pos <= r.end:
            return r
    return None


def genes(ref_gb):
    """The CDS features at global positions, names made unique, sorted by start, with a "record"
    column when there are several records (sequence_functions.reorder_reference_gbk)."""
    import sequence_functions as sf
    return sf.reorder_reference_gbk(ref_gb=ref_gb)


class Region:
    """A gene or intergenic region: name, ranges [(start, end)] of global positions, strand,
    record, kind ("gene", "intergenic" or "wrap")."""
    def __init__(self, name, ranges, strand, record, kind):
        self.name, self.ranges, self.strand, self.record, self.kind = name, ranges, strand, record, kind


def regions(ref, recs):
    """Genes (in ref's order), then each record's intergenic regions in order, ending with its
    wrap-round region "z:" (records without genes have none)."""
    out = [Region(n, [(float(s), float(e))], float(st), r, "gene")
           for n, s, e, st, r in zip(ref["name"], ref["start"], ref["end"], ref["strand"],
                                     ref["record"] if "record" in ref.columns else [recs[0].name] * len(ref))]
    for rec in recs:
        g = [x for x in out if x.record == rec.name and x.kind == "gene"]
        if not g:
            continue
        for a, b in zip(g, g[1:]):
            s, e = a.ranges[0][1] + 1, b.ranges[0][0] - 1
            if s <= e:
                out.append(Region(a.name + ":" + b.name, [(s, e)], 1.0, rec.name, "intergenic"))
        last, first = g[-1], g[0]
        ranges = []
        if last.ranges[0][1] + 1 <= rec.end:
            ranges.append((last.ranges[0][1] + 1, float(rec.end)))
        if first.ranges[0][0] - 1 >= rec.start:
            ranges.append((float(rec.start), first.ranges[0][0] - 1))
        out.append(Region(last.name + ":", ranges or [(last.ranges[0][1] + 1, float(rec.end))], 1.0, rec.name,
                          "wrap"))
    return out
