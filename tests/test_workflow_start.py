"""kmer_pipeline.nf up to its checks (container_type "none", the scripts of this checkout): the
parameters reach preflight.py, and Nextflow never ignores a parameter assignment (it warns
"defined multiple times" and keeps the first, which silently disabled precomputed_dir once)."""
import os
import shutil
import subprocess

import pytest

from conftest import EXAMPLE_DIR, REPO_DIR, SCRIPTS_DIR, SOFTWARE_LOCATION_TEMPLATE
from test_parameters import NEXTFLOW_RUNS

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
    shutil.copy(NF, base / "kmer_pipeline.nf")
    env = dict(os.environ, NXF_OPTS="-Dnxf.ansi.log=false")
    r = subprocess.run(["nextflow", "run", "kmer_pipeline.nf", "-ansi-log", "false"], cwd=base, env=env,
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
