"""finish() (PLAN.md Part B, B3): a pure-Groovy port of preflight.py's now-deleted finish(),
called from workflow.onComplete. Only guard (a) - normalizing workflow.sessionId (a UUID) and
the manifest's stored session (a String) to String before comparing - is in scope here; guard
(b), gating the onComplete call itself on a READY signal so a refused run can never overwrite a
previous successful run's manifest, needs Part B4's process-based preflight to exist first and
is deferred there (STATE.md/DECISIONS.md record this explicitly, not silently)."""
import json
import os
import subprocess

import pytest

from test_parameters import NEXTFLOW_RUNS, nf_block

pytestmark = pytest.mark.skipif(not NEXTFLOW_RUNS, reason="needs kmer_pipeline.nf and a working Nextflow")

PREFIX, TYPE, K = "tb20", "nucleotide", 31


def _manifest_path(analysis_dir):
    return os.path.join(analysis_dir, f"{PREFIX}_{TYPE}{K}.run_manifest.json")


def _run_finish(tmp_path, analysis_dir, session_id, status):
    block = nf_block("def finish", "}")
    probe = tmp_path / "probe.nf"
    probe.write_text(f"""
import groovy.json.JsonOutput
params.output_prefix = "{PREFIX}"
params.kmer_type = "{TYPE}"
params.kmer_length = {K}
ANALYSIS_DIR = "{analysis_dir}"
{block}
finish("{session_id}", "{status}")
""")
    env = dict(os.environ, NXF_OPTS="-Dnxf.ansi.log=false")
    return subprocess.run(["nextflow", "run", str(probe)], cwd=tmp_path, env=env,
                          capture_output=True, text=True)


def _write_manifest(analysis_dir, manifest):
    os.makedirs(analysis_dir, exist_ok=True)
    with open(_manifest_path(analysis_dir), "w") as f:
        json.dump(manifest, f)


@pytest.mark.parametrize("status", ["finished", "failed"])
def test_matching_session_updates_status(tmp_path, status):
    analysis = str(tmp_path / "analysis")
    _write_manifest(analysis, {"session": "s1", "status": "running"})
    result = _run_finish(tmp_path, analysis, "s1", status)
    assert result.returncode == 0, result.stderr
    manifest = json.loads(open(_manifest_path(analysis)).read())
    assert manifest["status"] == status
    assert "finished" in manifest  # the timestamp field, set alongside status


def test_session_mismatch_leaves_manifest_unchanged(tmp_path):
    analysis = str(tmp_path / "analysis")
    _write_manifest(analysis, {"session": "s1", "status": "running"})
    before = open(_manifest_path(analysis)).read()
    result = _run_finish(tmp_path, analysis, "s2", "finished")
    assert result.returncode == 0, result.stderr
    after = open(_manifest_path(analysis)).read()
    assert before == after


def test_sessionid_uuid_vs_manifest_string_still_matches(tmp_path):
    """Guards P5: workflow.sessionId is a UUID, not a String - finish() must normalize both
    sides before comparing, or a real run (always passing a UUID) could never match."""
    analysis = str(tmp_path / "analysis")
    uuid_str = "123e4567-e89b-12d3-a456-426614174000"
    _write_manifest(analysis, {"session": uuid_str, "status": "running"})
    block = nf_block("def finish", "}")
    probe = tmp_path / "probe.nf"
    probe.write_text(f"""
import groovy.json.JsonOutput
params.output_prefix = "{PREFIX}"
params.kmer_type = "{TYPE}"
params.kmer_length = {K}
ANALYSIS_DIR = "{analysis}"
{block}
finish(UUID.fromString("{uuid_str}"), "finished")
""")
    env = dict(os.environ, NXF_OPTS="-Dnxf.ansi.log=false")
    result = subprocess.run(["nextflow", "run", str(probe)], cwd=tmp_path, env=env,
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    manifest = json.loads(open(_manifest_path(analysis)).read())
    assert manifest["status"] == "finished"


def test_no_manifest_is_a_silent_no_op(tmp_path):
    analysis = str(tmp_path / "analysis")
    os.makedirs(analysis)
    result = _run_finish(tmp_path, analysis, "s1", "finished")
    assert result.returncode == 0, result.stderr
    assert not os.path.exists(_manifest_path(analysis))
