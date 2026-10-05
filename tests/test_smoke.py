"""Every script starts and every module imports."""
import importlib
import os

import pytest

from conftest import SCRIPTS_DIR, run_script

MODULES = ["rcompat", "sequence_functions", "alignmentfunctions", "Manhattan_functions", "inventory"]
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


def test_get_ref_name(example_dir):
    result = run_script("get_ref_name.py", "--fasta-file",
                        os.path.join(example_dir, "Mtub_H37Rv_NC000962.3.fasta"))
    assert result.returncode == 0, result.stderr
    assert result.stderr == ""
    assert result.stdout == "NC_000962.3"
