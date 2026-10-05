"""Shared set-up for the kmer_pipeline tests.

The scripts under test are found in KMER_PIPELINE_SCRIPTS if it is set; otherwise in
../python (a checkout of the repository) or, failing that, /usr/local/bin (the
kmer_pipeline image, where the tests are installed in /usr/share/kmer_pipeline/tests).
Run them with the image's Python, e.g. from a checkout:
    apptainer exec kmer_pipeline.sif python3 -m pytest -p no:cacheprovider tests
"""
import os
import subprocess
import sys

import pytest

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(TESTS_DIR)
EXAMPLE_DIR = os.path.join(REPO_DIR, "example")


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
