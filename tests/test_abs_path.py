"""absPath() (B9): the base_dir prefix of a path must be rewritten to BASE_DIR's own canonical
form, matching what nextflow.config actually binds into the container -- found as a real bug (the
container could not find id_file at all) by running a profile-based invocation with base_dir
reached through a symlink, exactly as on a shared HPC filesystem."""
import os
import subprocess

import pytest

from test_parameters import NEXTFLOW_RUNS, nf_block

PROBE = """
import java.nio.file.Paths
import java.nio.file.Files
""" + nf_block("def resolvePath(p) {", "}") + "\n" + nf_block("def absPath(p, canonical = false) {", "}") + """
BASE_DIR = absPath(params.base_dir, true)
println "RESULT " + absPath(params.user_path)
"""

pytestmark = pytest.mark.skipif(not NEXTFLOW_RUNS, reason="needs kmer_pipeline.nf and a working Nextflow")


def run_probe(tmp_path, base_dir, user_path):
    (tmp_path / "probe.nf").write_text(PROBE)
    env = dict(os.environ, NXF_OPTS="-Dnxf.ansi.log=false")
    r = subprocess.run(["nextflow", "run", "probe.nf", "--base_dir", base_dir, "--user_path", user_path],
                       cwd=tmp_path, env=env, capture_output=True, text=True, timeout=120)
    out = r.stdout + r.stderr
    assert r.returncode == 0, out
    line = next(l for l in out.split("\n") if l.startswith("RESULT "))
    return line[len("RESULT "):]


def test_plain_base_dir_is_unaffected(tmp_path):
    real = tmp_path / "base"
    real.mkdir()
    got = run_probe(tmp_path, str(real), str(real / "tb20" / "id_file.txt"))
    assert got == str(real / "tb20" / "id_file.txt")


def test_symlinked_base_dir_rewrites_the_prefix_to_the_real_spelling(tmp_path):
    """A path given via "$base_dir/..." interpolation (every id_file/ref_fa default in this
    repo's own example configs) inherits base_dir's own spelling verbatim -- if that spelling is
    a symlink, the file must still resolve to wherever the symlink's REAL target is, since that
    real target (not the symlink) is what nextflow.config's containerOptions actually binds."""
    real = tmp_path / "real"
    real.mkdir()
    link = tmp_path / "link"
    link.symlink_to(real)
    got = run_probe(tmp_path, str(link), str(link / "tb20" / "id_file.txt"))
    assert got == str(real / "tb20" / "id_file.txt")


def test_a_path_outside_base_dir_as_written_is_left_unresolved(tmp_path):
    """A symlinked assembly in a shared store elsewhere (bound via extra_binds, not base_dir)
    must not be silently followed -- unchanged from before this fix (guards P36/P52)."""
    base = tmp_path / "base"
    base.mkdir()
    elsewhere_real = tmp_path / "elsewhere_real"
    elsewhere_real.mkdir()
    elsewhere_link = tmp_path / "elsewhere_link"
    elsewhere_link.symlink_to(elsewhere_real)
    asm = elsewhere_link / "assembly.fa"
    asm.write_text("x")
    got = run_probe(tmp_path, str(base), str(asm))
    assert got == str(asm)
