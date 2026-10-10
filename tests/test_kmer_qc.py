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


def test_get_kmers_noresult_returns_none_when_everything_matched(tmp_path):
    kmers = af.Table(["kmer", "negLog10", "beta", "mac"], [["GGGT", 1.0, 0.3, 2]])
    res = af.Table(["kmer", "origkmer"], [["GGGT", "GGGT"]])
    assert af.get_kmers_noresult([res], kmers, "p", 1, ["gene"], "nucleotide", 4, str(tmp_path) + "/", "ref", 30) is None


def test_get_kmers_noresult_writes_the_table_it_returns(tmp_path):
    kmers = af.Table(["kmer", "negLog10", "beta", "mac"], [["AAAC", 5.0, 0.1, 3], ["GGGT", 1.0, 0.3, 2]])
    res = af.Table(["kmer", "origkmer"], [["GGGT", "GGGT"]])
    out_dir = str(tmp_path) + "/"
    t = af.get_kmers_noresult([res], kmers, "p", 1, ["my_gene"], "nucleotide", 4, out_dir, "ref", 30)
    path = out_dir + "p_nucleotide4_ref_top_gene_1_my_gene_no_blast_result_or_poor_alignment.txt"
    written = af.read_table_header(path)
    assert written.columns == t.columns
    assert written.col("kmer") == t.col("kmer")
    assert written.col("maf") == pytest.approx(t.col("maf"))


def test_run_manhattan_allframes_includes_no_blast_kmers_under_threshold(monkeypatch):
    expected = []
    monkeypatch.setattr(af, "FIGURES", type("F", (), {"expect": lambda self, path: expected.append(path)})())
    res = af.Table(["negLog10", "maf"], [[5.0, 0.5]])
    # Without the no-BLAST k-mer, nothing in the thresholded set exceeds 100, so no ylim50 variant.
    af.run_manhattan_allframes([res], "/p", "gene", None, 0.1, "maf")
    assert not any("ylim50" in p for p in expected)
    # A no-BLAST k-mer passing the MAF threshold, with negLog10 > 100, must push the thresholded
    # set's max over 100 too -- it was dropped from the threshold list before the fix.
    expected.clear()
    no_result = af.Table(["negLog10", "maf"], [[150.0, 0.5]])
    af.run_manhattan_allframes([res], "/p", "gene", no_result, 0.1, "maf")
    assert any(p.endswith("_allframes_Manhattan_ylim50_maf0.1.png") for p in expected)


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
