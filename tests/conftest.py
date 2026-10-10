"""Shared set-up for the kmer_pipeline tests.

The scripts under test are found in KMER_PIPELINE_SCRIPTS if it is set; otherwise in
../python (a checkout of the repository) or, failing that, /usr/local/bin (the
kmer_pipeline image, where the tests are installed in /usr/share/kmer_pipeline/tests).
Run them with the image's Python, e.g. from a checkout:
    apptainer exec --cleanenv kmer_pipeline.sif python3 -m pytest -p no:cacheprovider tests
"""
import glob
import os
import shutil
import subprocess
import sys

import pytest


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: a real pipeline run inside the image (NEXTFLOW_RUNS and gemma needed)")

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(TESTS_DIR)
EXAMPLE_DIR = os.path.join(REPO_DIR, "example")

# tests/e2e holds golden-output comparisons: they need an e2e config (local.conf
# or the same keys as environment variables) and are never run by a plain
# `pytest tests`. test_golden_hashes.py needs neither, so it lives under tests/
# instead of tests/e2e/, and is unaffected by this.
collect_ignore = []
if not os.environ.get("KMER_E2E_ROOT") and not os.path.exists(os.path.join(TESTS_DIR, "e2e", "local.conf")):
    collect_ignore.append("e2e")


def _scripts_dir():
    env = os.environ.get("KMER_PIPELINE_SCRIPTS")
    if env:
        return os.path.abspath(env)
    checkout = os.path.join(REPO_DIR, "python")
    if os.path.isdir(checkout):
        return checkout
    return "/usr/local/bin"


SCRIPTS_DIR = _scripts_dir()
if not os.path.isfile(os.path.join(SCRIPTS_DIR, "rcompat.py")):
    raise RuntimeError(f"pipeline scripts not found in {SCRIPTS_DIR}; set KMER_PIPELINE_SCRIPTS")
sys.path.insert(0, SCRIPTS_DIR)
sys.dont_write_bytecode = True


def pytest_report_header(config):
    return f"kmer_pipeline scripts under test: {SCRIPTS_DIR}"


@pytest.fixture
def scripts_dir():
    return SCRIPTS_DIR


@pytest.fixture
def example_dir():
    return EXAMPLE_DIR


def run_script(name, *args, cwd=None):
    """Run a pipeline script with the current Python; returns the CompletedProcess."""
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    return subprocess.run([sys.executable, os.path.join(SCRIPTS_DIR, name), *args],
                          capture_output=True, text=True, cwd=cwd, env=env)


# The image's compiled C++ tools (built separately; not part of this checkout). Kept in sync
# with the Dockerfile's build list.
CPP_TOOLS = ("kmerlist2pattern", "pattern2kinship", "patterncounts", "patternmerge", "sort_strings",
            "stringlist2count", "stringlist2pattern")


def stage_checkout(dest):
    """A flat scriptpath directory, as the image installs scripts in /usr/local/bin: this
    checkout's R and Python scripts, report assets, kmer_pipeline.nf, and symlinks to the
    image's C++ tools (not part of this checkout, so not staged from it). Python port of the
    private porting project's stage.sh, built from the working tree rather than `git archive`,
    so uncommitted test changes are exercised too."""
    os.makedirs(dest, exist_ok=True)
    for pattern in ("*.R", "*.Rscript"):
        for f in glob.glob(os.path.join(REPO_DIR, pattern)):
            out = os.path.join(dest, os.path.basename(f))
            shutil.copy(f, out)
            os.chmod(out, 0o755)
    for f in glob.glob(os.path.join(SCRIPTS_DIR, "*.py")):
        out = os.path.join(dest, os.path.basename(f))
        shutil.copy(f, out)
        os.chmod(out, 0o755)
    for name in ("report.css", "report.js", "kmer_pipeline.nf"):
        shutil.copy(os.path.join(REPO_DIR, name), os.path.join(dest, name))
    for tool in CPP_TOOLS:
        link = os.path.join(dest, tool)
        if not os.path.exists(link):
            os.symlink(os.path.join("/usr/local/bin", tool), link)
    return dest
