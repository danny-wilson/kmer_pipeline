"""reference.py: references with several records (D6)."""
import os

import pytest

from conftest import EXAMPLE_DIR

import reference
import sequence_functions as sf


def genbank_record(name, length, cds, topology="circular"):
    """A minimal GenBank record: cds is [(name, start, end, strand)]."""
    lines = [f"LOCUS       {name}            {length} bp    DNA     {topology}   BCT 01-JAN-2026",
             f"DEFINITION  {name} test record.", f"ACCESSION   {name}", f"VERSION     {name}.1",
             "FEATURES             Location/Qualifiers",
             f"     source          1..{length}"]
    for gene, start, end, strand in cds:
        loc = f"{start}..{end}" if strand > 0 else f"complement({start}..{end})"
        lines += [f"     CDS             {loc}", f'                     /gene="{gene}"',
                  f'                     /locus_tag="{gene}_tag"', f'                     /product="{gene} protein"']
    lines += ["ORIGIN", "//"]
    return "\n".join(lines) + "\n"


@pytest.fixture
def two_records(tmp_path):
    gb = tmp_path / "ref.gb"
    gb.write_text(genbank_record("chrom", 1000, [("dnaA", 101, 400, 1), ("tnpA", 501, 700, -1)])
                  + genbank_record("plas", 300, [("repA", 51, 150, 1), ("tnpA", 201, 260, 1)]))
    fa = tmp_path / "ref.fa"
    fa.write_text(">chrom.1 test\n" + "A" * 1000 + "\n>plas.1\n" + "C" * 300 + "\n")
    return str(fa), str(gb)


def test_records(two_records):
    fa, gb = two_records
    recs = reference.records(gb)
    assert [(r.name, r.length, r.offset, r.start, r.end) for r in recs] == \
        [("chrom", 1000, 0, 1, 1000), ("plas", 300, 1000, 1001, 1300)]
    assert reference.total_length(gb) == 1300 and reference.n_records(gb) == 2
    assert reference.check(fa, gb) == ([], [])
    assert reference.fasta_sequence(fa) == "A" * 1000 + "C" * 300


def test_check_mismatches(tmp_path, two_records):
    fa, gb = two_records
    short = tmp_path / "short.fa"
    short.write_text(">chrom.1\n" + "A" * 1000 + "\n>plas.1\n" + "C" * 299 + "\n")
    errors, _ = reference.check(str(short), gb)
    assert "299 bases" in errors[0]
    renamed = tmp_path / "renamed.fa"
    renamed.write_text(">contig_1\n" + "A" * 1000 + "\n>contig_2\n" + "C" * 300 + "\n")
    errors, warnings = reference.check(str(renamed), gb)
    assert errors == [] and len(warnings) == 2
    one = tmp_path / "one.fa"
    one.write_text(">chrom.1\n" + "A" * 1000 + "\n")
    assert "1 records" in reference.check(str(one), gb)[0][0]


def test_genes_and_regions(two_records):
    fa, gb = two_records
    ref = reference.genes(gb)
    assert list(zip(ref["name"], ref["start"], ref["end"], ref["record"])) == [
        ("dnaA", 101, 400, "chrom"), ("tnpA@chrom", 501, 700, "chrom"),
        ("repA", 1051, 1150, "plas"), ("tnpA@plas", 1201, 1260, "plas")]
    regs = reference.regions(ref, reference.records(gb))
    got = [(r.name, r.kind, r.ranges) for r in regs if r.kind != "gene"]
    assert got == [
        ("dnaA:tnpA@chrom", "intergenic", [(401, 500)]),
        ("tnpA@chrom:", "wrap", [(701, 1000), (1, 100)]),
        ("repA:tnpA@plas", "intergenic", [(1151, 1200)]),
        ("tnpA@plas:", "wrap", [(1261, 1300), (1001, 1050)])]


def test_single_record_genes_as_before():
    """The example reference: genes() equals reorder_reference_gbk (one record)."""
    fa = os.path.join(EXAMPLE_DIR, "Mtub_H37Rv_NC000962.3.fasta")
    gb = os.path.join(EXAMPLE_DIR, "Mtub_H37Rv_NC000962.3.gb")
    assert reference.n_records(gb) == 1 and reference.total_length(gb) == 4411532
    old = sf.reorder_reference_gbk(gb)
    new = reference.genes(gb)
    assert list(old["name"]) == list(new["name"])
    assert list(old["start"].astype(float)) == list(new["start"]) and list(old["end"].astype(float)) == list(new["end"])
