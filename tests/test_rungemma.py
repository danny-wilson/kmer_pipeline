"""rungemma.py (step 4) with the real GEMMA, through a wrapper that keeps a copy of the
phenotype file GEMMA is given and the library path it runs with."""
import gzip
import math
import os

import pytest

from conftest import EXAMPLE_DIR, run_script

import rungemma

pytestmark = pytest.mark.skipif(not os.path.exists("/usr/local/bin/gemma"),
                                reason="needs GEMMA (the kmer_pipeline image)")

WRAPPER = """#!/bin/bash
args=("$@")
for ((k=0; k<${{#args[@]}}; k++)); do
  [ "${{args[$k]}}" = -p ] && cp "${{args[$((k+1))]}}" {capture}/pheno.txt
done
echo "$LD_LIBRARY_PATH" > {capture}/ld_library_path.txt
exec /usr/local/bin/gemma "$@"
"""


def make_inputs(tmp_path, phenos, patterns, libraries=None):
    analysis = tmp_path / "analysis"
    analysis.mkdir()
    capture = tmp_path / "capture"
    capture.mkdir()
    ids = ["id\tpaths\tpheno"] + [f"S{i}\t/dev/null\t{v}" for i, v in enumerate(phenos)]
    (tmp_path / "ids.txt").write_text("\n".join(ids) + "\n")
    (analysis / "pre_protein5.patternmerge.patternKey.txt.gz").write_bytes(
        gzip.compress("".join(p + "\n" for p in patterns).encode()))
    (analysis / "pre_protein5.patternmerge.patternKeySize.txt").write_text(f"{len(patterns)}\n")
    n = len(phenos)
    kin = [[0.6 if i == j else 0.1 + 0.01 * ((i + j) % 3) for j in range(n)] for i in range(n)]
    (analysis / "pre_protein5.kinshipmerge.kinship.txt.gz").write_bytes(
        gzip.compress("".join(" ".join("%.17g" % v for v in r) + "\n" for r in kin).encode()))
    wrapper = tmp_path / "gemma.sh"
    wrapper.write_text(WRAPPER.format(capture=capture))
    wrapper.chmod(0o755)
    rows = open(os.path.join(EXAMPLE_DIR, "pipeline_software_location.txt")).read().rstrip("\n").split("\n")
    rows = [f"gemma\t{wrapper}" if r.startswith("gemma\t") else r for r in rows]
    if libraries:  # replaces the example's gemma_libraries entry (/usr/lib)
        rows = [r for r in rows if not r.startswith("gemma_libraries\t")] + [f"gemma_libraries\t{libraries}"]
    (tmp_path / "software.txt").write_text("\n".join(rows) + "\n")
    return analysis, capture


def run(tmp_path, analysis, env=None):
    old = dict(os.environ)
    os.environ.update(env or {})
    try:
        result = run_script("rungemma.py", "--task-id", "1", "--p", "1",
                            "--kmerfile-prefix", str(analysis / "pre_protein5"), "--id-file", str(tmp_path / "ids.txt"),
                            "--output-prefix", "pre", "--analysis-dir", str(analysis), "--kmertype", "protein",
                            "--kmer-length", "5", "--software-file", str(tmp_path / "software.txt"), cwd=tmp_path)
    finally:
        os.environ.clear()
        os.environ.update(old)
    assert result.returncode == 0, result.stdout[-2000:] + result.stderr[-2000:]
    return analysis / "proteinkmer5_gemma" / "output"


def test_phenotype_text():
    assert rungemma.gemma_phenotype_text(None) == "NA"
    assert rungemma.gemma_phenotype_text(math.nan) == "NA"
    assert rungemma.gemma_phenotype_text(math.inf) == "NA"
    assert rungemma.gemma_phenotype_text(3.0) == "3"
    assert rungemma.gemma_phenotype_text(-4.06) == "-4.0599999999999996"
    for v in (0.1234567891234, 123456.789, 1e-20, -1.7976931348623157e308):
        assert float(rungemma.gemma_phenotype_text(v)) == v


PHENOS = ["0.123456789", "123456.7891", "1e-05", "NA", "-2.5", "NaN", "3"]
PATTERNS = ["1010101", "0110011", "1111111", "0001110", "1000001"]


def test_full_precision_phenotypes_and_labelled_pvalues(tmp_path):
    """D3: GEMMA receives every finite phenotype exactly, NaN as NA. N2: the p-value file
    has a header and the pattern index of each tested pattern."""
    analysis, capture = make_inputs(tmp_path, PHENOS, PATTERNS)
    out = run(tmp_path, analysis)
    sent = (capture / "pheno.txt").read_text().split("\n")[:-1]
    assert [v == "NA" for v in sent] == [False, False, False, True, False, True, False]
    assert all(float(v) == float(p) for v, p in zip(sent, PHENOS) if v != "NA")
    with gzip.open(out / "pre_protein5.1-5.pval.txt.gz", "rt") as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh]
    assert rows[0] == ["rs", "p_lrt"]
    # Pattern 3 is present in every genome, so GEMMA does not test it
    assert [r[0] for r in rows[1:]] == ["1", "2", "4", "5"]
    assert all(0 <= float(r[1]) <= 1 for r in rows[1:])


def test_gemma_libraries_keep_the_library_path(tmp_path):
    """N3: a gemma_libraries entry adds to LD_LIBRARY_PATH rather than replacing it."""
    libs = tmp_path / "libs"
    libs.mkdir()
    analysis, capture = make_inputs(tmp_path, PHENOS, PATTERNS, libraries=libs)
    run(tmp_path, analysis, env={"LD_LIBRARY_PATH": "/opt/existing"})
    assert (capture / "ld_library_path.txt").read_text().strip() == f"/opt/existing:{libs}"
