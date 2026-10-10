"""plot_figures.R (the figures, drawn by R), run with Rscript."""
import math
import os
import re
import shutil
import subprocess
import sys

import pytest

from conftest import REPO_DIR, SCRIPTS_DIR

PLOT = next((p for p in (os.path.join(SCRIPTS_DIR, "plot_figures.R"), os.path.join(REPO_DIR, "plot_figures.R"))
             if os.path.exists(p)), None)
LAUNCHER = next((p for p in (os.path.join(SCRIPTS_DIR, "Rscript_launcher.R"), os.path.join(REPO_DIR, "Rscript_launcher.R"))
                 if os.path.exists(p)), None)
pytestmark = pytest.mark.skipif(PLOT is None or shutil.which("Rscript") is None,
                                reason="needs plot_figures.R and Rscript")

sys.path.insert(0, SCRIPTS_DIR)
import Manhattan_functions as mf  # noqa: E402


def run_r(code, tmp_path):
    (tmp_path / "t.R").write_text(f'source("{PLOT}")\n' + code)
    result = subprocess.run(["Rscript", "--vanilla", "t.R"], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return result.stdout


SINGLE_FRAME = """
captured = list()
plot_singleframe_manhattan_protein = function(beta_col, ...) captured[[length(captured) + 1]] <<- beta_col
res = data.frame(sstart = c(1, 5, 9), send = c(3, 7, 11), negLog10 = c(2, 5, 1), beta = c(0.5, -0.2, 0.3),
                 maf = c(0.2, 0.3, 0.4))
ref_gene_i = list(correct_frame = 1, length_correct = 100, all_translations = NULL, ref_start_i = 1, ref_end_i = 300)
run_manhattan_single_protein(which_kmers_no_result = {nores}, res = res, ref_gene_i = ref_gene_i, kmer_length = 3,
    prefix = "p", gene_name = "g", j = 1, bonferroni = 4, ref_length = 1000, kmer_type = "protein",
    minor_allele_threshold = 0.25, macormaf = "maf", output_dir = "", ref.name = "ref")
cat(captured[[1]], "\\n", captured[[2]], "\\n", sep = " ")
"""


def test_protein_single_frame_without_unaligned_kmers(tmp_path):
    """D2: with no unaligned k-mers, the last aligned k-mer keeps its beta colour."""
    out = run_r(SINGLE_FRAME.format(nores="NULL"), tmp_path).split("\n")
    assert out[0].split() == ["#838383", "#d3d3d3", "#838383"]  # all k-mers
    assert out[1].split() == ["#d3d3d3", "#838383"]  # MAF >= 0.25


def test_protein_single_frame_with_unaligned_kmers(tmp_path):
    nores = "data.frame(negLog10 = c(3, 1), beta = c(-0.1, 0.4), maf = c(0.3, 0.1))"
    out = run_r(SINGLE_FRAME.format(nores=nores), tmp_path).split("\n")
    assert out[0].split() == ["#838383", "#d3d3d3", "#838383", "#ffbc87", "#D55E00"]


def test_beta_colour_for_zero_and_nan_effects(tmp_path):
    """N10: an effect of exactly 0 or NaN is grey; before, no colour was returned."""
    out = run_r('cols = c("blue", "green", "orange", "red")\n'
                'b = rbind(c(0, -1, 0), c(0, 2, 0), c(0, 0, 0), c(0, NaN, 0), c(0, 1, 2))\n'
                'cat(apply(b, 1, function(x) get_betaCOL_2cols(x, cols, 1.96)), "\\n")\n', tmp_path)
    assert out.split() == ["blue", "red", "#d3d3d3", "#d3d3d3", "#d3d3d3"]


# --------------------------------------------------------------------------
# B1-R: an empty above-threshold set must not abort the single-frame plot
# --------------------------------------------------------------------------

SINGLE_FRAME_PROTEIN_COUNT = """
calls = 0
plot_singleframe_manhattan_protein = function(...) calls <<- calls + 1
res = data.frame(sstart = c(1, 5, 9), send = c(3, 7, 11), negLog10 = c(2, 5, 1), beta = c(0.5, -0.2, 0.3),
                 maf = c(0.01, 0.02, 0.03))
ref_gene_i = list(correct_frame = 1, length_correct = 100, all_translations = NULL, ref_start_i = 1, ref_end_i = 300)
run_manhattan_single_protein(which_kmers_no_result = NULL, res = res, ref_gene_i = ref_gene_i, kmer_length = 3,
    prefix = "p", gene_name = "g", j = 1, bonferroni = 4, ref_length = 1000, kmer_type = "protein",
    minor_allele_threshold = {thr}, macormaf = "maf", output_dir = "", ref.name = "ref")
cat(calls, "\\n")
"""

SINGLE_FRAME_NUCLEOTIDE_COUNT = """
calls = 0
plot_singleframe_manhattan_nucleotide = function(...) calls <<- calls + 1
res = data.frame(sstart = c(1, 5, 9), send = c(3, 7, 11), negLog10 = c(2, 5, 1), beta = c(0.5, -0.2, 0.3),
                 maf = c(0.01, 0.02, 0.03))
ref_gene_i = list(ref_gene_i = paste(rep("A", 300), collapse = ""), ref_start_i = 1, ref_end_i = 300)
run_manhattan_single_nucleotide(which_kmers_no_result = NULL, res = res, ref_gene_i = ref_gene_i, prefix = "p",
    gene_name = "g", bonferroni = 4, ref_length = 1000, kmer_type = "nucleotide", kmer_length = 31,
    minor_allele_threshold = {thr}, macormaf = "maf", output_dir = "", ref.name = "ref")
cat(calls, "\\n")
"""


@pytest.mark.parametrize("template", [SINGLE_FRAME_PROTEIN_COUNT, SINGLE_FRAME_NUCLEOTIDE_COUNT], ids=["protein", "nucleotide"])
def test_empty_threshold_set_draws_once_not_twice(template, tmp_path):
    """Before the fix, the threshold-only plot was attempted even with nothing passing the
    threshold; it is now skipped. A threshold above every k-mer's maf must give exactly 1 call
    (all k-mers only), not 2 (the empty threshold plot too)."""
    out = run_r(template.format(thr=0.5), tmp_path)
    assert out.strip() == "1"


@pytest.mark.parametrize("template", [SINGLE_FRAME_PROTEIN_COUNT, SINGLE_FRAME_NUCLEOTIDE_COUNT], ids=["protein", "nucleotide"])
def test_nonempty_threshold_set_still_draws_twice(template, tmp_path):
    """Sanity check for the test above: when some k-mers do pass, both plots are still drawn."""
    out = run_r(template.format(thr=0.015), tmp_path)
    assert out.strip() == "2"


# --------------------------------------------------------------------------
# B3a: the launcher contract (Rscript_launcher.R), and generic plot_figures.R properties
# (ported from the private porting project's harness, scrubbed of local paths/IDs)
# --------------------------------------------------------------------------

pytestmark_launcher = pytest.mark.skipif(LAUNCHER is None, reason="needs Rscript_launcher.R")


def run_launcher(args, cwd):
    return subprocess.run(["Rscript", "--vanilla", LAUNCHER] + args, cwd=cwd, capture_output=True, text=True)


@pytestmark_launcher
def test_launcher_exits_1_when_main_calls_stop(tmp_path):
    bad = tmp_path / "bad.R"
    bad.write_text('main = function(args) stop("boom")\n')
    p = run_launcher([str(bad)], tmp_path)
    assert p.returncode == 1
    assert "boom" in p.stderr
    assert re.search(r"bad\.R:\d+", p.stderr), p.stderr


@pytestmark_launcher
def test_errors_report_file_and_line(tmp_path):
    """A missing data directory and a table missing a column: exit status 1, the message,
    and a plot_figures.R file:line in the traceback."""
    p = run_launcher([PLOT, "--data-dir", str(tmp_path / "nowhere")], tmp_path)
    assert p.returncode == 1
    assert "Figure data directory doesn't exist" in p.stderr
    assert re.search(r"plot_figures\.R:\d+", p.stderr), p.stderr

    fd = mf.FigureData(str(tmp_path) + "/")
    fd.table("patterns", [("neglog10p", "numeric")], [(1.0,)])  # beta, maf, ma missing
    fd.table("kmers", [("kmer_index", "integer"), ("ma", "numeric")], [(1, 0.5)])
    for k, v in (("figures_dir", str(tmp_path) + "/"), ("output_prefix", "x"), ("kmer_type", "nucleotide"),
                 ("kmer_length", 31), ("ref_name", "R"), ("ref_length", 1000.0), ("macormaf", "maf"),
                 ("minor_allele_threshold", 0.01), ("bonferroni", 3.0), ("pheno_type", "continuous"), ("nsamples", 3),
                 ("override_signif", False), ("manhattan_stem", str(tmp_path) + "/x")):
        fd.param(k, v)
    fd.close()
    p = run_launcher([PLOT, "--data-dir", fd.dir], tmp_path)
    assert p.returncode == 1
    assert "Traceback" in p.stderr and re.search(r"plot_figures\.R:\d+", p.stderr), p.stderr

    p = run_launcher([PLOT], tmp_path)
    assert p.returncode == 1 and "Usage: plot_figures.R" in p.stderr


def test_round_trip(tmp_path):
    """Python writes a table of awkward values; R's reader gives them back exactly."""
    fd = mf.FigureData(str(tmp_path) + "/")
    strings = ["NA", "T", "F", "TRUE", "#comment", "'quoted'", '"double"', "", " lead", "trail ", "a,b", "NaN", "Inf",
               None, "\\n not a newline", "ACGT*-"]
    numbers = [1 / 3, -0.0, 0.0, 1e-300, 5e-324, 1.7976931348623157e308, math.inf, -math.inf, math.nan, None,
               5.62983468081783, 100000.0, 3.5611013836490559, 2 ** 53 + 1.0, -123.456, 1e5]
    ints = list(range(-3, len(strings) - 3))
    logic = [True, False, None] + [True] * (len(strings) - 3)
    rows = list(zip(strings, numbers, ints, logic))
    fd.table("t", [("s", "character"), ("x", "numeric"), ("i", "integer"), ("b", "logical")], rows)
    code = f'''source("{PLOT}")
d = read_fd("{fd.dir}", "t")
enc = function(v) ifelse(is.na(v) & !is.nan(v), "<NA>", v)
writeLines(c(enc(d$s), enc(sprintf("%a", d$x)), enc(as.character(d$i)), enc(as.character(d$b)), class(d$s), class(d$x),
             class(d$i), class(d$b)), "out.txt", useBytes = TRUE)
'''
    (tmp_path / "t.R").write_text(code)
    p = subprocess.run(["Rscript", "--vanilla", "t.R"], cwd=tmp_path, capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    out = (tmp_path / "out.txt").read_text().split("\n")
    n = len(rows)
    assert out[:n] == ["<NA>" if v is None else v for v in strings]
    rx = out[n:2 * n]
    for v, r in zip(numbers, rx):
        if v is None:  # sprintf("%a", NA) gives "NA"
            assert r in ("NA", "<NA>")
        elif math.isnan(v):  # Python's NaN is the pipeline's missing value: written as NA
            assert r in ("NA", "<NA>")
        else:
            assert float.fromhex(r.replace("Inf", "inf")) == v and math.copysign(1, float.fromhex(r.replace("Inf", "inf"))) \
                == math.copysign(1, v), (v, r)
    assert out[2 * n:3 * n] == [str(v) for v in ints]
    assert out[3 * n:4 * n] == ["<NA>" if v is None else ("TRUE" if v else "FALSE") for v in logic]
    assert out[4 * n:4 * n + 4] == ["character", "numeric", "integer", "logical"]


def test_writer_refuses_unwritable_strings(tmp_path):
    fd = mf.FigureData(str(tmp_path) + "/")
    for bad in ("a\tb", "a\nb", mf.NA_SENTINEL):
        with pytest.raises(Exception):
            fd.table("bad", [("s", "character")], [(bad,)])


def test_r_does_no_data_handling():
    """plot_figures.R only reads the figure-data files Python wrote for it; it must not open
    anything else (a GenBank file, another script, a shell command) on its own."""
    src = open(PLOT).read()
    code = "\n".join(l.split("#")[0] for l in src.split("\n"))
    for forbidden in (r"\bsystem2?\s*\(", r"read_dna_seg", r"library\s*\(", r"require\s*\(", r"\bpipe\s*\(",
                      r"read\.table\s*\(", r"read\.delim\s*\(", r"readRDS\s*\(", r"\bsource\s*\("):
        assert not re.search(forbidden, code), forbidden


def test_each_sample_is_seeded():
    """Every sample() call in plot_figures.R is preceded by set.seed(0), so figure layout
    (e.g. the unmapped-kmer x positions) is reproducible run to run."""
    lines = open(PLOT).read().split("\n")
    calls = [k for k, l in enumerate(lines) if re.search(r"\bsample\(", l.split("#")[0])]
    assert calls
    for k in calls:
        before = lines[k].split("sample(")[0]
        assert "set.seed(0)" in before or "set.seed(0)" in lines[k - 1], lines[k]
