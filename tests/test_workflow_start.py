"""kmer_pipeline.nf up to its checks (container_type "none", the scripts of this checkout): the
parameters reach preflight.py, and Nextflow never ignores a parameter assignment (it warns
"defined multiple times" and keeps the first, which silently disabled precomputed_dir once)."""
import os
import shutil
import subprocess

import pytest

from conftest import EXAMPLE_DIR, REPO_DIR, SCRIPTS_DIR, SOFTWARE_LOCATION_TEMPLATE
from test_parameters import NEXTFLOW_RUNS
from test_preflight import STEP_FILES, P

NF = os.path.join(REPO_DIR, "kmer_pipeline.nf")
pytestmark = pytest.mark.skipif(not NEXTFLOW_RUNS or not os.path.exists(NF), reason="needs Nextflow and the .nf")


def make_base(tmp_path, extra=""):
    base = tmp_path / "base"
    (base / "tb20").mkdir(parents=True)
    rows = ["id\tpaths\tpheno"]
    for k, sid in enumerate(("702", "725", "791", "10084")):
        src = os.path.join(EXAMPLE_DIR, f"{sid}.velvet_assembly.contigs.fa.gz")
        shutil.copy(src, base / "tb20" / os.path.basename(src))
        rows.append(f"{sid}\t{base}/tb20/{os.path.basename(src)}\t{k + 1}")
    (base / "tb20" / "id_file.txt").write_text("\n".join(rows) + "\n")
    for f in ("Mtub_H37Rv_NC000962.3.fasta", "Mtub_H37Rv_NC000962.3.gb"):
        shutil.copy(os.path.join(EXAMPLE_DIR, f), base / "tb20" / f)
    sw = open(SOFTWARE_LOCATION_TEMPLATE).read().rstrip("\n").split("\n")
    sw.append(f"scriptpath\t{SCRIPTS_DIR}")
    (base / "software.txt").write_text("\n".join(sw) + "\n")
    (base / "nextflow.config").write_text(f"""params {{
	base_dir = "{base}"
	output_prefix = "tb20"
	analysis_dir = "{base}/tb20/kmergwas"
	kmer_type = "nucleotide"
	kmer_length = 31
	id_file = "{base}/tb20/id_file.txt"
	ref_fa = "{base}/tb20/Mtub_H37Rv_NC000962.3.fasta"
	ref_gb = "{base}/tb20/Mtub_H37Rv_NC000962.3.gb"
	maxp = 2
	software_file = "{base}/software.txt"
{extra}
}}
""")
    return base


def run(base):
    # Run the repo's own kmer_pipeline.nf by its real path (PLAN.md Part B, B8), rather than
    # copying just the .nf into base: Nextflow's projectDir - and therefore bin/'s automatic
    # PATH entry, and every bare-name script call that relies on it - is always the main
    # script's own directory, regardless of cwd (which stays base, for base_dir-relative
    # resolution). Copying only the .nf would lose bin/ entirely, silently falling back to
    # whatever script of the same name happens to be on PATH otherwise (the toolchain image's
    # own stale copy, in every real run).
    env = dict(os.environ, NXF_OPTS="-Dnxf.ansi.log=false")
    r = subprocess.run(["nextflow", "run", NF, "-ansi-log", "false"], cwd=base, env=env,
                       capture_output=True, text=True, timeout=600)
    return r.returncode, r.stdout + r.stderr


def test_precomputed_dir_and_pheno_file_reach_the_checks(tmp_path):
    base = make_base(tmp_path)
    (base / "pre").mkdir()
    (base / "tb20" / "pheno.txt").write_text("id\tpheno\n702\t1\n725\t0\n791\t1\n10084\t0\n")
    with open(base / "nextflow.config", "a") as fh:
        fh.write(f'params.precomputed_dir = "{base}/pre"\nparams.pheno_file = "{base}/tb20/pheno.txt"\n')
    rc, out = run(base)
    assert "defined multiple times" not in out, out
    assert "steps 1-3 and 5 skipped" in out
    assert rc != 0 and "lacks outputs this analysis reads" in out, out[-3000:]
    assert not (base / "tb20" / "kmergwas").exists()  # nothing written before the checks pass


def test_running_a_precomputed_step_is_an_error(tmp_path):
    base = make_base(tmp_path, extra="\tskip3 = false")
    (base / "pre").mkdir()
    with open(base / "nextflow.config", "a") as fh:
        fh.write(f'params.precomputed_dir = "{base}/pre"\n')
    rc, out = run(base)
    assert rc != 0 and "skip3 = false: with precomputed_dir" in out, out[-3000:]


def test_rerun_into_an_existing_analysis_stops(tmp_path):
    base = make_base(tmp_path)
    kmer_dir = base / "tb20" / "kmergwas" / "nucleotidekmer31"
    kmer_dir.mkdir(parents=True)
    (kmer_dir / "702.kmer31.txt.gz").write_text("x\n")
    rc, out = run(base)
    assert "defined multiple times" not in out
    assert rc != 0 and "set overwrite = true" in out, out[-3000:]
    assert (kmer_dir / "702.kmer31.txt.gz").exists()


def test_symlinked_base_dir_is_accepted(tmp_path):
    """base_dir given as a symlink, with the inputs given in the real (target) spelling: the
    two spellings of the same place must be treated as equal, not rejected as outside base_dir."""
    base = make_base(tmp_path)
    link = tmp_path / "link"
    link.symlink_to(base)
    cfg_path = base / "nextflow.config"
    cfg_path.write_text(cfg_path.read_text().replace(f'base_dir = "{base}"', f'base_dir = "{link}"'))
    # Reuse the existing-analysis trick so the run stops deterministically at a known preflight
    # check, proving path resolution succeeded before ever reaching it.
    kmer_dir = base / "tb20" / "kmergwas" / "nucleotidekmer31"
    kmer_dir.mkdir(parents=True)
    (kmer_dir / "702.kmer31.txt.gz").write_text("x\n")
    rc, out = run(base)
    assert "Error converting from user_path" not in out, out[-3000:]
    assert "is not beneath base_dir" not in out, out[-3000:]
    assert rc != 0 and "set overwrite = true" in out, out[-3000:]


def _trace_task_names(trace_path):
    lines = trace_path.read_text().strip().split("\n")
    header = lines[0].split("\t")
    name_col = header.index("name")
    return [line.split("\t")[name_col] for line in lines[1:]]


@pytest.mark.parametrize("first_step", [1, 2, 3, 4, 5, 6, 7, "precomputed_dir"])
def test_preflight_error_leaves_a_trace_of_only_preflight(tmp_path, first_step):
    """A forced preflight error (an earlier run's output already in analysis_dir) must leave a
    trace containing only the preflight task - no other process may have started, whichever step
    is effectively the first one to run. Parametrized over every "first running step" (guards
    P11: skip1..skip3-only coverage would miss the roots reached via prepareGemma/plotManhattan
    that precomputed_dir forces by skipping 1, 2, 3 and 5)."""
    base = make_base(tmp_path)
    extra_lines = []
    if first_step == "precomputed_dir":
        pre = base / "pre"
        pre.mkdir()
        extra_lines.append(f'params.precomputed_dir = "{pre}"')
        trigger_step = 4  # 1, 2, 3, 5 are forced-skipped; prepareGemma (step 4) is the real root
    else:
        for s in range(1, first_step):
            extra_lines.append(f'params.skip{s} = true')
        trigger_step = first_step
    with open(base / "nextflow.config", "a") as fh:
        fh.write("\n".join(extra_lines) + "\n")
    analysis = base / "tb20" / "kmergwas"
    for f in STEP_FILES[trigger_step]:
        path = analysis / f
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("x\n")
    trace = base / "trace.tsv"
    env = dict(os.environ, NXF_OPTS="-Dnxf.ansi.log=false")
    r = subprocess.run(["nextflow", "run", NF, "-with-trace", str(trace)],
                       cwd=base, env=env, capture_output=True, text=True, timeout=600)
    out = r.stdout + r.stderr
    assert r.returncode != 0 and "set overwrite = true" in out, out[-3000:]
    assert trace.exists(), out[-3000:]
    assert _trace_task_names(trace) == ["preflight"], _trace_task_names(trace)


def test_overwrite_deletes_before_the_rerun_step_starts(tmp_path):
    """overwrite = true must delete a stale step's outputs inside the preflight process, before
    READY is ever emitted - not after deployment_write() has already written fresh files that a
    later deletion pass might then confuse for stale ones. Confirmed two ways: no "set overwrite =
    true" error (deletion succeeded) and the step's own stale marker file is gone immediately,
    without needing that step to actually run again."""
    base = make_base(tmp_path)
    with open(base / "nextflow.config", "a") as fh:
        fh.write("params.overwrite = true\n")
        for s in range(2, 8):
            fh.write(f"params.skip{s} = true\n")
    analysis = base / "tb20" / "kmergwas"
    stale_marker = analysis / STEP_FILES[1][0]
    for f in STEP_FILES[1]:
        path = analysis / f
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("stale\n")
    assert stale_marker.exists()
    rc, out = run(base)
    assert "set overwrite = true" not in out, out[-3000:]
    assert "Deleted" in out, out[-3000:]
    # Whatever countkmers itself does next (it is step 1, so it runs for real) is not this test's
    # concern - only that preflight's own deletion already happened before it could start.
    assert not stale_marker.exists() or stale_marker.read_bytes() != b"stale\n", out[-3000:]


def test_preflight_stderr_lands_in_command_err(tmp_path):
    """Guards P34: output: stdout only parses the preflight process's stdout (the JSON result's
    last line, exactly as the old Groovy .execute() parsed sout) - its stderr is no longer
    captured by that parsing at all, unlike the old code's explicit waitForProcessOutput(sout,
    serr). Confirms the Nextflow guarantee this depends on still holds for this exact process:
    anything the script writes to stderr lands in its own .command.err, inspectable even though
    it never reaches the .map{} closure. preflight.py itself has no stderr-writing code path to
    force today, so a copy of the real kmer_pipeline.nf has one spliced into the preflight
    process's own script (everything else about the process, and the rest of the pipeline,
    unchanged) - this is still the real process definition, not a hand-written substitute."""
    base = make_base(tmp_path)
    nf_text = open(NF).read()
    marker = 'preflight.py --analysis-dir "${ANALYSIS_DIR}"'
    assert marker in nf_text
    nf_text = nf_text.replace(marker, 'echo "a preflight warning, on stderr" >&2\n' + marker)
    (base / "kmer_pipeline.nf").write_text(nf_text)
    # This copy needs modified content, so it can't run NF by path (B8); projectDir is then base
    # itself, so bin/ needs to be reachable there too, for the bare preflight.py call to resolve
    # to this checkout rather than the toolchain image's own stale copy.
    (base / "bin").symlink_to(SCRIPTS_DIR)
    with open(base / "nextflow.config", "a") as fh:
        for s in range(2, 8):
            fh.write(f"params.skip{s} = true\n")
    trace = base / "trace.tsv"
    env = dict(os.environ, NXF_OPTS="-Dnxf.ansi.log=false")
    r = subprocess.run(["nextflow", "run", "kmer_pipeline.nf", "-with-trace", str(trace)],
                       cwd=base, env=env, capture_output=True, text=True, timeout=600)
    out = r.stdout + r.stderr
    lines = trace.read_text().strip().split("\n")
    header = lines[0].split("\t")
    hash_col, name_col = header.index("hash"), header.index("name")
    preflight_rows = [l.split("\t") for l in lines[1:] if l.split("\t")[name_col] == "preflight"]
    assert len(preflight_rows) == 1, out[-3000:]
    prefix, rest = preflight_rows[0][hash_col].split("/")
    candidates = list((base / "work" / prefix).glob(f"{rest}*"))
    assert len(candidates) == 1, candidates
    command_err = (candidates[0] / ".command.err").read_text()
    assert "a preflight warning, on stderr" in command_err, command_err


@pytest.mark.parametrize("key, bad_value, expected", [
    ("base_dir", "/no/such/base_dir", "base_dir /no/such/base_dir does not exist"),
    ("id_file", "/no/such/id_file.txt", "id_file /no/such/id_file.txt does not exist"),
    ("ref_fa", "/no/such/ref.fa", "ref_fa /no/such/ref.fa does not exist"),
])
def test_three_representative_error_messages_are_unchanged(tmp_path, key, bad_value, expected):
    """These three deployment() existence checks kept their exact wording through the B4
    rewrite (they are not container-path-translation concerns, just Files.exists() checks) -
    confirmed directly against the real pipeline, not just by reading the source."""
    base = make_base(tmp_path)
    cfg_path = base / "nextflow.config"
    text = cfg_path.read_text()
    if key == "base_dir":
        text = text.replace(f'base_dir = "{base}"', f'base_dir = "{bad_value}"')
    elif key == "id_file":
        text = text.replace(f'id_file = "{base}/tb20/id_file.txt"', f'id_file = "{bad_value}"')
    elif key == "ref_fa":
        text = text.replace(f'ref_fa = "{base}/tb20/Mtub_H37Rv_NC000962.3.fasta"', f'ref_fa = "{bad_value}"')
    cfg_path.write_text(text)
    rc, out = run(base)
    assert rc != 0 and expected in out, out[-3000:]


@pytest.mark.parametrize("key, value", [
    ("container_type", "singularity"),
    ("container_file", "/some/old/kmer_pipeline_2026-10-06.sif"),
    ("container_args", "--bind /scratch"),
    ("container_cmd", "singularity exec"),
    ("container_mount", "/mnt"),
    ("container_script_dir", "/usr/local/bin"),
])
def test_legacy_container_key_is_a_guided_migration_error(tmp_path, key, value):
    """PLAN.md Part C, C2 (OQ-S3-5: a clean break with a guided error, not silent acceptance):
    any of the six pre-refactor container_* keys stops the run before anything else happens -
    no analysis_dir, no workdir/logdir - with exit status 1 and a message naming the key and
    pointing at -profile/-params-file/-c and the manual's migration section."""
    base = make_base(tmp_path, extra=f'\t{key} = "{value}"')
    rc, out = run(base)
    assert rc == 1, (rc, out[-3000:])
    assert key in out and "-profile" in out and "docs/manual.md" in out, out[-3000:]
    assert not (base / "tb20" / "kmergwas").exists(), "analysis_dir must not be created"


def test_resume_reruns_preflight_not_cached(tmp_path):
    """preflight must never show CACHED under -resume, since a cached preflight would silently
    skip the resume-validity checks. All steps skipped, so the first run succeeds with nothing but
    preflight ever running; -resume must still show it COMPLETED.

    Checked directly (not assumed): dropping cache false from the process does not, on its own,
    make this test fail, because workflow.sessionId is interpolated into the --session-id
    argument of preflight's own script, so the rendered script - and therefore the task hash
    Nextflow caches against - differs on every invocation regardless of the cache setting; there
    is nothing for a cache hit to match even without it. cache false is kept anyway, as explicit,
    correct documentation of intent that does not depend on sessionId always being part of the
    script (guards a future edit that might remove it from there)."""
    base = make_base(tmp_path)
    with open(base / "nextflow.config", "a") as fh:
        for s in range(1, 8):
            fh.write(f"params.skip{s} = true\n")
    env = dict(os.environ, NXF_OPTS="-Dnxf.ansi.log=false")
    trace1 = base / "trace1.tsv"
    r1 = subprocess.run(["nextflow", "run", NF, "-with-trace", str(trace1)],
                        cwd=base, env=env, capture_output=True, text=True, timeout=120)
    out1 = r1.stdout + r1.stderr
    assert r1.returncode == 0, out1[-3000:]
    status_col = trace1.read_text().strip().split("\n")[0].split("\t").index("status")
    row1 = trace1.read_text().strip().split("\n")[1].split("\t")
    assert row1[status_col] == "COMPLETED", out1[-3000:]

    trace2 = base / "trace2.tsv"
    r2 = subprocess.run(["nextflow", "run", NF, "-with-trace", str(trace2), "-resume"],
                        cwd=base, env=env, capture_output=True, text=True, timeout=120)
    out2 = r2.stdout + r2.stderr
    assert r2.returncode == 0, out2[-3000:]
    row2 = trace2.read_text().strip().split("\n")[1].split("\t")
    assert row2[status_col] == "COMPLETED", out2[-3000:]


def test_relative_paths_and_trailing_slash_are_accepted(tmp_path):
    """analysis_dir with a trailing slash, and analysis_dir/id_file/ref_fa given as paths
    relative to the launch directory (base, since run() invokes nextflow with cwd=base):
    absPath() must resolve all of these correctly, not just an absolute as-written path
    (guards P47)."""
    base = make_base(tmp_path)
    cfg_path = base / "nextflow.config"
    text = cfg_path.read_text()
    text = text.replace(f'analysis_dir = "{base}/tb20/kmergwas"', 'analysis_dir = "tb20/kmergwas/"')
    text = text.replace(f'id_file = "{base}/tb20/id_file.txt"', 'id_file = "tb20/id_file.txt"')
    text = text.replace(f'ref_fa = "{base}/tb20/Mtub_H37Rv_NC000962.3.fasta"', 'ref_fa = "tb20/Mtub_H37Rv_NC000962.3.fasta"')
    cfg_path.write_text(text)
    # Reuse the existing-analysis trick so the run stops deterministically at a known preflight
    # check, proving path resolution succeeded before ever reaching it.
    kmer_dir = base / "tb20" / "kmergwas" / "nucleotidekmer31"
    kmer_dir.mkdir(parents=True)
    (kmer_dir / "702.kmer31.txt.gz").write_text("x\n")
    rc, out = run(base)
    assert "Error converting from user_path" not in out, out[-3000:]
    assert "is not beneath base_dir" not in out, out[-3000:]
    assert "does not exist" not in out, out[-3000:]
    assert rc != 0 and "set overwrite = true" in out, out[-3000:]


def test_input_outside_base_dir_fails_clearly(tmp_path):
    base = make_base(tmp_path)
    outside = tmp_path / "outside"
    outside.mkdir()
    bad_ref = outside / "Mtub_H37Rv_NC000962.3.fasta"
    shutil.copy(os.path.join(EXAMPLE_DIR, "Mtub_H37Rv_NC000962.3.fasta"), bad_ref)
    cfg_path = base / "nextflow.config"
    cfg_path.write_text(cfg_path.read_text().replace(
        f'ref_fa = "{base}/tb20/Mtub_H37Rv_NC000962.3.fasta"', f'ref_fa = "{bad_ref}"'))
    rc, out = run(base)
    assert rc != 0 and "is not beneath base_dir" in out, out[-3000:]
