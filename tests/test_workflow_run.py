"""Real pipeline runs, inside the image: the fixes that can only be exercised by Nextflow
actually running the workflow (not a preflight-only run, and not a pure Python/R unit test).

Needs the in-image toolchain (dsk, GEMMA, BLAST, MUMmer, R, the compiled C++ tools) and a
working Nextflow; most tests here are marked `slow` since they run the full 30-genome example.
Run with: apptainer exec --cleanenv kmer_pipeline.sif python3 -m pytest tests -m "not slow"
to skip them, or without -m to run everything.
"""
import glob
import json
import os
import shutil
import subprocess

import pytest

from conftest import EXAMPLE_DIR, REPO_DIR, SCRIPTS_DIR, stage_checkout
from test_parameters import NEXTFLOW_RUNS

NF = os.path.join(REPO_DIR, "kmer_pipeline.nf")
GEMMA_PRESENT = os.path.exists("/usr/local/bin/gemma")
pytestmark = pytest.mark.skipif(not NEXTFLOW_RUNS or not os.path.exists(NF),
                                reason="needs Nextflow and the .nf")

GENOME_IDS = ("702", "725", "791", "10084")  # the 4-genome subset test_workflow_start.py also uses


def run_nf(base, timeout=600):
    env = dict(os.environ, NXF_OPTS="-Dnxf.ansi.log=false")
    r = subprocess.run(["nextflow", "run", "kmer_pipeline.nf", "-ansi-log", "false"], cwd=base, env=env,
                       capture_output=True, text=True, timeout=timeout)
    return r.returncode, r.stdout + r.stderr


def write_software_file(base, scriptpath):
    sw = open(os.path.join(EXAMPLE_DIR, "pipeline_software_location.txt")).read().split("\n")
    sw = [f"scriptpath\t{scriptpath}" if line.startswith("scriptpath\t") else line for line in sw]
    (base / "software.txt").write_text("\n".join(sw))


def write_config(base, params_block):
    (base / "nextflow.config").write_text("params {\n" + params_block + "}\n")
    shutil.copy(NF, base / "kmer_pipeline.nf")


def read_top_genes(analysis_dir, n=None):
    """The gene names plotManhattan.py ranked (Manhattan_functions.top20genes), read from its
    own output file, written before any close-up figure or report is drawn. Deliberately not
    read from the close-up figures/reports B3b/B13 check, so a bug confined to those (as B3b's
    mutant is) cannot also corrupt the expected list (guards P25)."""
    matches = glob.glob(os.path.join(str(analysis_dir), "**", "*_top20genes_toppvals_*.txt"), recursive=True)
    assert len(matches) == 1, f"expected exactly one top20genes file under {analysis_dir}, found {matches}"
    genes = [line.split("\t")[0] for line in open(matches[0]) if line.strip()]
    return genes[:n] if n is not None else genes


def figures_dir_for(analysis_dir, kmer_type, kmer_length):
    matches = glob.glob(os.path.join(str(analysis_dir), "**", f"*{kmer_type}kmer{kmer_length}_kmergenealign_figures"),
                        recursive=True)
    assert len(matches) == 1, f"expected exactly one figures dir, found {matches}"
    return matches[0]


# --------------------------------------------------------------------------
# B12: min_contig_length actually drops short contigs, inside a real Nextflow run
# --------------------------------------------------------------------------


def make_minimal_base(tmp_path, extra=""):
    """The 4-genome subset, scriptpath pointing at a staged (not bare python/) directory, since
    step 1 (countkmers) needs sort_strings and dsk colocated with the scripts."""
    base = tmp_path / "base"
    (base / "tb20").mkdir(parents=True)
    rows = ["id\tpaths\tpheno"]
    for k, sid in enumerate(GENOME_IDS):
        src = os.path.join(EXAMPLE_DIR, f"{sid}.velvet_assembly.contigs.fa.gz")
        shutil.copy(src, base / "tb20" / os.path.basename(src))
        rows.append(f"{sid}\t{base}/tb20/{os.path.basename(src)}\t{k + 1}")
    (base / "tb20" / "id_file.txt").write_text("\n".join(rows) + "\n")
    for f in ("Mtub_H37Rv_NC000962.3.fasta", "Mtub_H37Rv_NC000962.3.gb"):
        shutil.copy(os.path.join(EXAMPLE_DIR, f), base / "tb20" / f)
    staged = stage_checkout(str(tmp_path / "staging"))
    write_software_file(base, staged)
    write_config(base, f"""\tbase_dir = "{base}"
\toutput_prefix = "tb20"
\tanalysis_dir = "{base}/tb20/kmergwas"
\tkmer_type = "nucleotide"
\tkmer_length = 31
\tid_file = "{base}/tb20/id_file.txt"
\tref_fa = "{base}/tb20/Mtub_H37Rv_NC000962.3.fasta"
\tref_gb = "{base}/tb20/Mtub_H37Rv_NC000962.3.gb"
\tmaxp = 2
\tsoftware_file = "{base}/software.txt"
{extra}""")
    return base


def test_min_contig_length_drops_short_contigs(tmp_path):
    threshold = 400
    base_filtered = make_minimal_base(tmp_path / "filtered", extra=f"""\tmin_contig_length = {threshold}
\tskip2 = true
\tskip3 = true
\tskip4 = true
\tskip5 = true
\tskip6 = true
\tskip7 = true
""")
    rc, out = run_nf(base_filtered, timeout=300)
    assert rc == 0, out[-3000:]
    logdir = base_filtered / "tb20" / "kmergwas" / "log.tb20_nucleotide31"
    logs = "".join(open(p).read() for p in glob.glob(str(logdir / "countkmers.*.log")))
    assert f"Minimum contig length {threshold}" in logs and "dropped" in logs, logs[-3000:]

    base_unfiltered = make_minimal_base(tmp_path / "unfiltered", extra="""\tskip2 = true
\tskip3 = true
\tskip4 = true
\tskip5 = true
\tskip6 = true
\tskip7 = true
""")
    rc2, out2 = run_nf(base_unfiltered, timeout=300)
    assert rc2 == 0, out2[-3000:]

    kmer_dir_filtered = base_filtered / "tb20" / "kmergwas" / "nucleotidekmer31"
    kmer_dir_unfiltered = base_unfiltered / "tb20" / "kmergwas" / "nucleotidekmer31"
    def total_kmers(path):
        return float(open(path).read().strip().split("\t")[1])  # "Total\t<n>"

    for sid in GENOME_IDS:
        filtered_count = total_kmers(kmer_dir_filtered / f"{sid}.kmer31.total.txt")
        unfiltered_count = total_kmers(kmer_dir_unfiltered / f"{sid}.kmer31.total.txt")
        assert filtered_count < unfiltered_count, (sid, filtered_count, unfiltered_count)
    assert not list(tmp_path.glob("**/*.minlen.fa"))  # the filtered copy is removed on exit
