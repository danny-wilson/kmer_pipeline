"""read_id_file() (PLAN.md Part B, B6): relative paths in id_file's paths column now resolve
against id_file's own directory, not the launch directory - a deliberate, user-visible change
(guards P27). CRLF, a trailing blank line, and a wrong column count (naming the line) are all
handled explicitly; zero data rows is covered by test_get_n.py (guards P21), not duplicated here.
Extracts the live function from kmer_pipeline.nf via nf_block, as the other Part B test files do."""
import json
import os
import subprocess

import pytest

from test_parameters import NEXTFLOW_RUNS, nf_block

pytestmark = pytest.mark.skipif(not NEXTFLOW_RUNS, reason="needs kmer_pipeline.nf and a working Nextflow")


def _run_read_id_file(tmp_path, id_file_path):
    block = nf_block("def read_id_file", "}")
    probe = tmp_path / "probe.nf"
    probe.write_text(f"""
import groovy.json.JsonOutput
import java.nio.file.Paths
params.id_file = "{id_file_path}"
def absPath(p, canonical = false) {{ Paths.get(p.toString()).toAbsolutePath().normalize().toString() }}
{block}
try {{
	def result = read_id_file()
	println "RESULT:" + JsonOutput.toJson(result)
}} catch (Exception e) {{
	System.err.println(e.getMessage())
	System.exit(1)
}}
""")
    env = dict(os.environ, NXF_OPTS="-Dnxf.ansi.log=false")
    return subprocess.run(["nextflow", "run", str(probe)], cwd=tmp_path, env=env,
                          capture_output=True, text=True)


def _result(proc):
    line = next(l for l in proc.stdout.split("\n") if l.startswith("RESULT:"))
    return json.loads(line[len("RESULT:"):])


def test_relative_path_resolves_against_id_file_directory_not_launch_dir(tmp_path):
    """Guards P27: a deliberate, user-visible change from launch-directory resolution."""
    sub = tmp_path / "data"
    sub.mkdir()
    (sub / "a.fa").write_text(">a\nACGT\n")
    (sub / "id_file.txt").write_text("id\tpaths\tpheno\nA\ta.fa\t1\n")
    # Launch nextflow from tmp_path, NOT sub - if resolution used the launch dir, "a.fa" would
    # resolve to tmp_path/a.fa (which doesn't exist) instead of sub/a.fa.
    result = _result(_run_read_id_file(tmp_path, sub / "id_file.txt"))
    assert result["paths"] == [str(sub / "a.fa")]


def test_absolute_path_is_kept_as_is(tmp_path):
    sub = tmp_path / "data"
    sub.mkdir()
    abs_target = tmp_path / "elsewhere.fa"
    abs_target.write_text(">x\nACGT\n")
    (sub / "id_file.txt").write_text(f"id\tpaths\tpheno\nA\t{abs_target}\t1\n")
    result = _result(_run_read_id_file(tmp_path, sub / "id_file.txt"))
    assert result["paths"] == [str(abs_target)]


def test_crlf_line_endings(tmp_path):
    id_file = tmp_path / "id_file.txt"
    id_file.write_bytes(b"id\tpaths\tpheno\r\nA\ta.fa\t1\r\nB\tb.fa\t2\r\n")
    result = _result(_run_read_id_file(tmp_path, id_file))
    assert result["id"] == ["A", "B"]
    assert result["pheno"] == ["1", "2"]  # no stray \r stuck to the last column


def test_trailing_blank_line_is_dropped(tmp_path):
    id_file = tmp_path / "id_file.txt"
    id_file.write_text("id\tpaths\tpheno\nA\ta.fa\t1\n\n")
    result = _result(_run_read_id_file(tmp_path, id_file))
    assert result["id"] == ["A"]


def test_wrong_column_count_names_the_line(tmp_path):
    id_file = tmp_path / "id_file.txt"
    id_file.write_text("id\tpaths\tpheno\nA\ta.fa\t1\nB\tb.fa\n")  # line 3: only 2 columns
    proc = _run_read_id_file(tmp_path, id_file)
    assert proc.returncode != 0
    assert "line 3" in proc.stderr
    assert "2" in proc.stderr  # names the actual column count found


def test_wrong_header_is_a_clear_error(tmp_path):
    id_file = tmp_path / "id_file.txt"
    id_file.write_text("id\tpath\tphenotype\nA\ta.fa\t1\n")
    proc = _run_read_id_file(tmp_path, id_file)
    assert proc.returncode != 0
    assert "line 1" in proc.stderr and "header" in proc.stderr


def test_empty_pheno_column_is_not_dropped(tmp_path):
    """split('\\t', -1), not split('\\t'): Java's default split drops a trailing empty column,
    which would otherwise make a legitimately-empty pheno look like a wrong column count."""
    id_file = tmp_path / "id_file.txt"
    id_file.write_text("id\tpaths\tpheno\nA\ta.fa\t\n")
    result = _result(_run_read_id_file(tmp_path, id_file))
    assert result["pheno"] == [""]
