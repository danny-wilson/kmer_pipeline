"""Every script starts and every module imports."""
import importlib
import os

import pytest

from conftest import SCRIPTS_DIR, run_script

MODULES = ["rcompat", "sequence_functions", "alignmentfunctions", "Manhattan_functions", "inventory", "reference"]
SCRIPTS = sorted(f for f in os.listdir(SCRIPTS_DIR)
                 if f.endswith(".py") and f[:-3] not in MODULES)


@pytest.mark.parametrize("module", MODULES)
def test_module_imports(module):
    importlib.import_module(module)


@pytest.mark.parametrize("script", SCRIPTS)
def test_script_help(script):
    result = run_script(script, "--help")
    assert result.returncode == 0, result.stderr
    assert "usage:" in result.stdout

# get_ref_name.py was ported to pure Groovy (kmer_pipeline.nf's ref_name_from_fasta(), PLAN.md
# Part B, B3) and deleted; see tests/test_ref_name.py.
