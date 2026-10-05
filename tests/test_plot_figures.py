"""plot_figures.R (the figures, drawn by R), run with Rscript."""
import os
import shutil
import subprocess

import pytest

from conftest import REPO_DIR, SCRIPTS_DIR

PLOT = next((p for p in (os.path.join(SCRIPTS_DIR, "plot_figures.R"), os.path.join(REPO_DIR, "plot_figures.R"))
             if os.path.exists(p)), None)
pytestmark = pytest.mark.skipif(PLOT is None or shutil.which("Rscript") is None,
                                reason="needs plot_figures.R and Rscript")


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
