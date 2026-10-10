"""Per-kmer QC warnings (sequence_functions.kmer_qc_*) and the MAF of kmers without a BLAST result."""
import pytest

from conftest import SCRIPTS_DIR  # noqa: F401

import alignmentfunctions as af
import sequence_functions as sf


def test_homopolymer_is_flagged_and_ordinary_kmer_is_not():
    assert sf.kmer_qc_flags("A" * 31) == ["low_complexity"]
    assert sf.kmer_qc_flags("A" * 30 + "C") == ["low_complexity"]
    assert sf.kmer_qc_flags("ACGTACGTTGCATGCAAGCTTCGAGCTAGCT") == []
    assert sf.kmer_qc_flags("") == []


def test_threshold_is_at_ninety_percent():
    assert sf.kmer_qc_flags("A" * 9 + "C") == ["low_complexity"]
    assert sf.kmer_qc_flags("A" * 8 + "CG") == []


def test_html_marks_flagged_kmers_only():
    assert "<sup" in sf.kmer_qc_html("A" * 31) and sf.kmer_qc_html("A" * 31).startswith("A" * 31)
    assert sf.kmer_qc_html("ACGTACGTTGCATGCAAGCTTCGAGCTAGCT") == "ACGTACGTTGCATGCAAGCTTCGAGCTAGCT"
    assert sf.kmer_qc_html(None) == "NA"


def test_legend_only_when_something_is_flagged():
    assert sf.kmer_qc_legend(["ACGTACGTTGCATGCAAGCTTCGAGCTAGCT"]) == ""
    assert "low complexity" in sf.kmer_qc_legend(["A" * 31, "ACGT"])


def test_kmers_without_a_blast_result_get_a_maf(tmp_path):
    kmers = af.Table(["kmer", "negLog10", "beta", "mac"], [["AAAC", 5.0, 0.1, 3], ["CCCG", 9.0, -0.2, 6], ["GGGT", 1.0, 0.3, 2]])
    res = af.Table(["kmer", "origkmer"], [["GGGT", "GGGT"]])
    t = af.get_kmers_noresult([res], kmers, "p", 1, ["gene"], "nucleotide", 4, str(tmp_path) + "/", "ref", 30)
    assert t.col("kmer") == ["CCCG", "AAAC"]  # most significant first
    assert t.col("maf") == pytest.approx([6 / 30, 3 / 30])


def test_a_gene_whose_only_threshold_kmers_lack_a_blast_result_still_has_points():
    res = af.Table(["negLog10", "mac", "maf"], [[2.0, 1, 1 / 30]])
    none = af.Table(["kmer", "negLog10", "beta", "mac", "maf"], [["AAAC", 50.0, 0.1, 6, 0.2]])
    ypos, which = af._manhattan_values(res, none, "maf", 0.01)
    assert ypos == [2.0, 50.0] and which == [0, 1]
    ypos, which = af._manhattan_values(res, none, "maf", 0.1)
    assert which == [1]


def test_an_empty_threshold_plot_warns_instead_of_stopping(capsys, monkeypatch):
    res = af.Table(["negLog10", "mac", "maf"], [[2.0, 1, 1 / 30]])
    monkeypatch.setattr(af, "FIGURES", type("F", (), {"expect": lambda self, path: None})())
    af.run_manhattan_single(res, None, "/p", "geneX", 5.0, "maf", 0.5, "")
    out = capsys.readouterr().out
    assert "no k-mers to plot for" in out and "geneX" in out
