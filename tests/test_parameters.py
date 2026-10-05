"""Workflow parameters: kmer_min_count and plot_min_genomes replace min_count (D4)."""
import os
import re
import shutil
import subprocess
import tempfile

import pytest

from conftest import REPO_DIR, SCRIPTS_DIR, run_script

NF = next((p for p in (os.path.join(REPO_DIR, "kmer_pipeline.nf"), os.path.join(SCRIPTS_DIR, "kmer_pipeline.nf"))
           if os.path.exists(p)), None)


def nextflow_runs():
    if NF is None or shutil.which("nextflow") is None:
        return False
    try:
        # run in a temporary folder: the launcher leaves nxf-tmp.* files in its working directory
        with tempfile.TemporaryDirectory() as d:
            return subprocess.run(["nextflow", "-version"], cwd=d, capture_output=True, timeout=120).returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        return False


NEXTFLOW_RUNS = nextflow_runs()


def nf_block(start, end):
    """Lines of kmer_pipeline.nf from the line starting with start to the one starting with end."""
    lines = open(NF).read().split("\n")
    i = next(k for k, line in enumerate(lines) if line.startswith(start))
    j = next(k for k in range(i, len(lines)) if lines[k].startswith(end))
    return "\n".join(lines[i:j + 1]) + "\n"


@pytest.mark.skipif(not NEXTFLOW_RUNS, reason="needs kmer_pipeline.nf and a working Nextflow")
@pytest.mark.parametrize("args, expected", [
    ([], "kmer_min_count=1 plot_min_genomes=1"),
    (["--min_count", "3"], "kmer_min_count=3 plot_min_genomes=3"),
    (["--kmer_min_count", "2"], "kmer_min_count=2 plot_min_genomes=1"),
    (["--plot_min_genomes", "4"], "kmer_min_count=1 plot_min_genomes=4"),
    (["--min_count", "3", "--plot_min_genomes", "2"], "error"),
])
def test_min_count_parameters(tmp_path, args, expected):
    block = nf_block("// D4: min_count", "println 'plot_min_genomes:")
    (tmp_path / "probe.nf").write_text(
        block + 'println "RESULT kmer_min_count=${params.kmer_min_count} plot_min_genomes=${params.plot_min_genomes}"\n')
    env = dict(os.environ, NXF_OPTS="-Dnxf.ansi.log=false")
    result = subprocess.run(["nextflow", "run", "probe.nf", *args], cwd=tmp_path, env=env, capture_output=True,
                            text=True)
    out = result.stdout + result.stderr
    assert "defined multiple times" not in out
    if expected == "error":
        assert result.returncode != 0 and "set only the new parameters" in out
    else:
        assert result.returncode == 0, out
        assert re.search("RESULT " + expected + r"\b", out), out
        assert ("min_count is replaced" in out) == ("--min_count" in args)


@pytest.mark.parametrize("script, new, old", [
    ("stringlist2patternandkinship.py", "--kmer-min-count", "--mincount"),
    ("plotManhattan.py", "--plot-min-genomes", "--min-count"),
    ("gen-report.py", "--plot-min-genomes", "--mincount"),
    ("gen-gene-report.py", "--plot-min-genomes", "--mincount"),
    ("gen-protein-report.py", "--plot-min-genomes", "--mincount"),
    ("gen-unmapped-report.py", "--plot-min-genomes", "--mincount"),
])
def test_script_option_names(script, new, old):
    """The new option name, with the old one still accepted."""
    result = run_script(script, "--help")
    assert result.returncode == 0
    assert re.search(re.escape(new) + r" \S+, " + re.escape(old) + r" \S+", result.stdout), result.stdout
