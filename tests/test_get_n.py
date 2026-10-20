"""get_n() (PLAN.md Part B, B3): pure Groovy, DEPLOYMENT_FILES.id_list.id.size(), not a
container wc -l - and an explicit error on zero data rows (guards P21), not a silent fall-through
to Channel.of(1..0) (Groovy's descending-range reading of that call). Extracts the live function
from kmer_pipeline.nf, as test_ref_name.py does."""
import os
import subprocess

import pytest

from test_parameters import NEXTFLOW_RUNS, nf_block

pytestmark = pytest.mark.skipif(not NEXTFLOW_RUNS, reason="needs kmer_pipeline.nf and a working Nextflow")


def _run_with_ids(tmp_path, ids):
    block = nf_block("def get_n", "}")
    id_list_groovy = "[" + ", ".join(f"'{i}'" for i in ids) + "]"
    probe = tmp_path / "probe.nf"
    probe.write_text(f"""
DEPLOYMENT_FILES = [id_list: [id: {id_list_groovy}]]
{block}
try {{
	println "RESULT:" + get_n()
}} catch (Exception e) {{
	System.err.println(e.getMessage())
	System.exit(1)
}}
""")
    env = dict(os.environ, NXF_OPTS="-Dnxf.ansi.log=false")
    return subprocess.run(["nextflow", "run", str(probe)], cwd=tmp_path, env=env,
                          capture_output=True, text=True)


@pytest.mark.parametrize("ids", [["a"], ["a", "b"], ["a", "b", "c", "d", "e"]])
def test_counts_data_rows(tmp_path, ids):
    result = _run_with_ids(tmp_path, ids)
    assert result.returncode == 0, result.stderr
    line = next(l for l in result.stdout.split("\n") if l.startswith("RESULT:"))
    assert line == f"RESULT:{len(ids)}"


def test_zero_data_rows_errors_explicitly(tmp_path):
    """Guards P21: Channel.of(1..0) would otherwise silently become the descending range [1, 0]."""
    result = _run_with_ids(tmp_path, [])
    assert result.returncode != 0
    assert "no data rows" in result.stderr
