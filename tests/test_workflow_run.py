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
import re
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


def run_nf(base, timeout=600, extra_args=()):
    env = dict(os.environ, NXF_OPTS="-Dnxf.ansi.log=false")
    r = subprocess.run(["nextflow", "run", "kmer_pipeline.nf", "-ansi-log", "false", *extra_args], cwd=base, env=env,
                       capture_output=True, text=True, timeout=timeout)
    return r.returncode, r.stdout + r.stderr


def add_params(base, lines):
    """Insert extra `params {}` lines into an existing nextflow.config (written by write_config,
    which always ends the file with a bare "}\\n" closing the params block), before rerunning."""
    config = (base / "nextflow.config").read_text()
    assert config.endswith("}\n"), config[-40:]
    (base / "nextflow.config").write_text(config[:-2] + "".join(line + "\n" for line in lines) + "}\n")


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


# --------------------------------------------------------------------------
# full_run / full_run_protein: the 30-genome example, module-scoped so B3b, B5, B6, B7 and B13
# share one real run each (nucleotide and protein) rather than paying for it per test.
# --------------------------------------------------------------------------

ALL_GENOME_IDS = [line.split("\t")[0] for line in open(os.path.join(EXAMPLE_DIR, "id_file.txt")).read().split("\n")[1:]
                  if line.strip()]


def make_full_base(base_dir, kmer_type, kmer_length, ntopgenes, staged, extra=""):
    base_dir.mkdir(parents=True)
    (base_dir / "tb20").mkdir()
    lines = open(os.path.join(EXAMPLE_DIR, "id_file.txt")).read().split("\n")
    header, rows = lines[0], [l for l in lines[1:] if l.strip()]
    out_rows = [header]
    for row in rows:
        sid, _, pheno = row.split("\t")
        src = os.path.join(EXAMPLE_DIR, f"{sid}.velvet_assembly.contigs.fa.gz")
        shutil.copy(src, base_dir / "tb20" / os.path.basename(src))
        out_rows.append(f"{sid}\t{base_dir}/tb20/{os.path.basename(src)}\t{pheno}")
    (base_dir / "tb20" / "id_file.txt").write_text("\n".join(out_rows) + "\n")
    for f in ("Mtub_H37Rv_NC000962.3.fasta", "Mtub_H37Rv_NC000962.3.gb"):
        shutil.copy(os.path.join(EXAMPLE_DIR, f), base_dir / "tb20" / f)
    write_software_file(base_dir, staged)
    write_config(base_dir, f"""\tbase_dir = "{base_dir}"
\toutput_prefix = "tb20"
\tanalysis_dir = "{base_dir}/tb20/kmergwas"
\tkmer_type = "{kmer_type}"
\tkmer_length = {kmer_length}
\tid_file = "{base_dir}/tb20/id_file.txt"
\tref_fa = "{base_dir}/tb20/Mtub_H37Rv_NC000962.3.fasta"
\tref_gb = "{base_dir}/tb20/Mtub_H37Rv_NC000962.3.gb"
\tmaxp = 2
\tntopgenes = {ntopgenes}
\tsoftware_file = "{base_dir}/software.txt"
{extra}""")
    return base_dir


NTOPGENES_BASE = 5
FULL_RUN_TIMEOUT = 30 * 60  # PLAN.md: at least 30 min, expected ~10-12 min
GEMMA_GATE = pytest.mark.skipif(not GEMMA_PRESENT, reason="needs the in-image gemma binary")


@pytest.fixture(scope="module")
def staged_checkout(tmp_path_factory):
    return stage_checkout(str(tmp_path_factory.mktemp("staging")))


@pytest.fixture(scope="module")
def full_run(tmp_path_factory, staged_checkout):
    base = make_full_base(tmp_path_factory.mktemp("full_run") / "base", "nucleotide", 31, NTOPGENES_BASE,
                          staged_checkout)
    rc, out = run_nf(base, timeout=FULL_RUN_TIMEOUT)
    assert rc == 0, out[-5000:]
    return base


@pytest.fixture(scope="module")
def full_run_protein(tmp_path_factory, staged_checkout):
    base = make_full_base(tmp_path_factory.mktemp("full_run_protein") / "base", "protein", 11, NTOPGENES_BASE,
                          staged_checkout)
    rc, out = run_nf(base, timeout=FULL_RUN_TIMEOUT)
    assert rc == 0, out[-5000:]
    return base


def read_trace(path):
    lines = open(path).read().strip().split("\n")
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def process_name(trace_row):
    return trace_row["name"].split(" ")[0].split("(")[0].strip()


def analysis_dir_of(base):
    return base / "tb20" / "kmergwas"


# --------------------------------------------------------------------------
# B3b: the baseline full run produced every report and figure, with no silently-skipped gene
# and no broken links. Must run first (file order): B7 and B5/B6 change the fixture's state.
#
# EXPECTED_CLOSEUP_PNGS is captured from one known-good run of this exact fixture (nucleotide
# k=31, maxp=2, ntopgenes=5, default blastident/minor_allele_threshold) -- not recomputed from
# the same figure output this test checks (that would be the P25 circularity the plan warns
# against), and not from first principles either (replicating plotManhattan's own alignment-window
# math was judged too large a reimplementation to safely duplicate correctly). A gene's count is
# deterministic given fixed example data and a fixed reference (the handoff's own evidence: 0 PNG
# differences run to run with a fixed image). If the example data or reference genome ever change,
# this needs recapturing -- see C7's identical precedent for the same reason.
# --------------------------------------------------------------------------

EXPECTED_CLOSEUP_PNGS = {"rpoB": 10, "PE_PGRS2": 2, "lpqW": 2, "PE_PGRS52": 2, "PE_PGRS38": 2}


@GEMMA_GATE
def test_full_run_has_every_report_and_figure(full_run):
    analysis_dir = analysis_dir_of(full_run)
    genes = read_top_genes(analysis_dir, n=NTOPGENES_BASE)
    assert len(genes) == NTOPGENES_BASE, genes
    assert set(genes) == set(EXPECTED_CLOSEUP_PNGS), (genes, "fixture's top genes changed; recapture the dict")

    figures_dir = figures_dir_for(analysis_dir, "nucleotide", 31)
    for gene in genes:
        for suffix in ("_Manhattan_allkmers.png", "_Manhattan_maf0.01.png"):
            assert glob.glob(os.path.join(figures_dir, f"*_{gene}{suffix}")), (gene, suffix)
        pngs = glob.glob(os.path.join(figures_dir, f"*_{gene}_*.png"))
        assert len(pngs) == EXPECTED_CLOSEUP_PNGS[gene], (gene, sorted(pngs))

        report = glob.glob(str(analysis_dir / "**" / f"*.report_{gene}.html"), recursive=True)
        assert len(report) == 1, (gene, report)
        html = open(report[0]).read()
        for m in re.finditer(r'(?:src|href)="([^"]+)"', html):
            target = m.group(1)
            if target.startswith(("http:", "https:", "#")):
                continue
            assert target != "NA", (gene, "src/href=\"NA\"")
            assert os.path.exists(os.path.join(os.path.dirname(report[0]), target)), (gene, target)

    main_report = glob.glob(str(analysis_dir / "**" / "*.report.html"), recursive=True)
    assert len(main_report) == 1, main_report
    unmapped_report = glob.glob(str(analysis_dir / "**" / "*.report_unmapped.html"), recursive=True)
    assert len(unmapped_report) == 1, unmapped_report

    for log_path in glob.glob(str(analysis_dir / "log.*" / "*.log")):
        text = open(log_path, errors="replace").read()
        assert not re.search(r"Execution halted|^Error in |^Error:", text, re.M), log_path
    assert not glob.glob(str(analysis_dir / "**" / "Rcoredump.rda"), recursive=True)


# --------------------------------------------------------------------------
# B7: -resume reuses steps 1-5's statistics when only ntopgenes/blastident change, and
# refuses a changed min_contig_length. Runs before B5/B6: it is what bumps ntopgenes to 6.
# --------------------------------------------------------------------------

NTOPGENES_AFTER_B7 = NTOPGENES_BASE + 1
CACHE_REUSED = {"countkmers", "createfullkmerlist", "stringlist2patternandkinship", "rungemma", "kmercontigalign",
               "kmercontigalignmerge"}
RERUN_ON_FIGURE_CHANGE = {"plotManhattan", "plotFigures", "genReport", "genGeneReport", "genUnmappedReport"}


@GEMMA_GATE
def test_resume_figure_only_parameter_reuses_statistics(full_run):
    add_params(full_run, [f"\tntopgenes = {NTOPGENES_AFTER_B7}"])
    rc, out = run_nf(full_run, timeout=300, extra_args=("-resume", "-with-trace", "trace2.txt"))
    assert rc == 0, out[-5000:]

    trace = read_trace(full_run / "trace2.txt")
    by_name = {}
    for r in trace:
        by_name.setdefault(process_name(r), []).append(r["status"])
    for name in CACHE_REUSED:
        assert name in by_name and all(s == "CACHED" for s in by_name[name]), (name, by_name.get(name))
    for name in RERUN_ON_FIGURE_CHANGE:
        assert name in by_name and all(s == "COMPLETED" for s in by_name[name]) and "CACHED" not in by_name[name], \
            (name, by_name.get(name))

    analysis_dir = analysis_dir_of(full_run)
    genes = read_top_genes(analysis_dir, n=NTOPGENES_AFTER_B7)
    assert len(genes) == NTOPGENES_AFTER_B7
    for gene in genes:
        report = glob.glob(str(analysis_dir / "**" / f"*.report_{gene}.html"), recursive=True)
        assert len(report) == 1, (gene, report)

    # Second resume: blastident changes too (closes the handoff's "not exercised end to end" gap).
    add_params(full_run, ["\tblastident = 80"])
    rc, out = run_nf(full_run, timeout=300, extra_args=("-resume",))
    assert rc == 0, out[-5000:]

    # Content, not just task status (guards P26): zero-CACHED rows on the resume proves the
    # report tasks reran, but a rerun could still read a stale intermediate file left over from
    # the previous blastident. A stale file is distinguishable from a correct one by its own
    # content: every row's BLAST percent identity (the "pident" column alignmentfunctions.py
    # keeps in the written table) must meet the *new* threshold -- a row below it is exactly
    # what a stale pre-change file would still contain.
    figures_dir = figures_dir_for(analysis_dir, "nucleotide", 31)
    blast_result_files = glob.glob(os.path.join(figures_dir, "*_blast_results.txt"))
    assert blast_result_files
    checked_rows = 0
    for path in blast_result_files:
        header, *rows = open(path).read().rstrip("\n").split("\n")
        cols = header.split("\t")
        pident_idx = cols.index("pident")
        for row in rows:
            pident = float(row.split("\t")[pident_idx])
            assert pident >= 80, (path, row, "pident below the new blastident=80 threshold -- stale content")
            checked_rows += 1
    assert checked_rows, "no BLAST-matched rows found to check at all"


# --------------------------------------------------------------------------
# B13: annotateGeneFile draws close-ups for named genes outside the current top set, without
# changing which genes the step-7B HTML reports cover (B8's documented limitation).
#
# Deliberately placed here (right after B7, before B5/B6), not last as PLAN.md's prose lists
# it: this needs -resume (see the comment below), and Nextflow refuses -resume across a
# run_steps (skip-flag) change. B7 leaves skip1-5=true/skip6=false/skip7=false; that is the
# state this test needs and must not disturb. B5/B6 (next) changes skip6 to true -- a
# deliberate run_steps change, made without -resume (overwrite=true instead, which Nextflow
# does allow), so it must come after, not before.
# --------------------------------------------------------------------------


@GEMMA_GATE
def test_annotate_gene_file_adds_close_ups_without_changing_reports(full_run):
    analysis_dir = analysis_dir_of(full_run)
    current_top = read_top_genes(analysis_dir, n=NTOPGENES_AFTER_B7)
    all_ranked = read_top_genes(analysis_dir)
    extra_genes = [g for g in all_ranked if g not in current_top][:2]
    assert len(extra_genes) >= 1, all_ranked

    figures_dir = figures_dir_for(analysis_dir, "nucleotide", 31)

    def closeups_for(gene):
        return set(glob.glob(os.path.join(figures_dir, f"*_{gene}_*.png")))

    for g in extra_genes:
        assert not closeups_for(g), (g, "already has close-ups before annotateGeneFile -- not outside the top set")

    report_files_before = sorted(glob.glob(str(analysis_dir / "**" / "*.report_*.html"), recursive=True))

    # -resume, not skip-flags + overwrite: overwrite deletes existing outputs first, and
    # annotateGeneFile replaces (does not add to) plot_closeup_alignments's gene list
    # (kmer_pipeline.nf passes only the annotateGeneFile genes to it when set -- checked in
    # plotManhattan.py). Deleting first would wipe the normal top-ntopgenes genes' own
    # close-up data, which this run would then never regenerate (it only processes the named
    # genes), breaking gen-gene-report.py for every top-ntopgenes gene. -resume reruns the
    # same steps without deleting anything first, which is also what the manual documents
    # ("changing annotateGeneFile under -resume reruns only the figure and report steps").
    # No skip-flag or overwrite override needed here: B7 (just above) already left
    # skip1-5=true/skip6=false/skip7=false/overwrite=false, exactly what this needs.
    annotate_file = full_run / "tb20" / "annotate.txt"
    annotate_file.write_text("\n".join(extra_genes) + "\n")
    add_params(full_run, [f'\tannotateGeneFile = "{annotate_file}"'])
    rc, out = run_nf(full_run, timeout=300, extra_args=("-resume",))
    assert rc == 0, out[-5000:]

    for g in extra_genes:
        assert closeups_for(g), (g, "no close-up figures were drawn for a named gene outside ntopgenes")

    report_files_after = sorted(glob.glob(str(analysis_dir / "**" / "*.report_*.html"), recursive=True))
    assert report_files_before == report_files_after, (report_files_before, report_files_after)
    reported_genes = {os.path.basename(p).split(".report_")[1][:-len(".html")]
                      for p in report_files_after if ".report_" in p} - {"unmapped"}
    assert reported_genes == set(current_top), (reported_genes, current_top)
    assert not (reported_genes & set(extra_genes)), (reported_genes, extra_genes)


# --------------------------------------------------------------------------
# B5 / B6: a stage-7-alone rerun (skip1..skip6, not skip7) still reports every gene, and
# submits nothing but the step-7 processes. Runs after B13 (see above): this makes its own
# deliberate run_steps change (skip6 true) via overwrite=true, not -resume, so it is free to
# run in either order relative to B13; placed last since nothing after it depends on its state.
# --------------------------------------------------------------------------


@GEMMA_GATE
def test_stage7_alone_reports_every_gene_and_submits_nothing_else(full_run):
    """B5 (process names) and B6 nucleotide (report count), run together since both read the
    same trace from the same rerun."""
    analysis_dir = analysis_dir_of(full_run)
    for html in glob.glob(str(analysis_dir / "**" / "*.report_*.html"), recursive=True):
        os.remove(html)
    genes = read_top_genes(analysis_dir, n=NTOPGENES_AFTER_B7)

    add_params(full_run, ["\tskip1 = true", "\tskip2 = true", "\tskip3 = true", "\tskip4 = true", "\tskip5 = true",
                          "\tskip6 = true", "\tskip7 = false", "\toverwrite = true"])
    rc, out = run_nf(full_run, timeout=300, extra_args=("-with-trace", "trace3.txt"))
    assert rc == 0, out[-5000:]

    trace = read_trace(full_run / "trace3.txt")
    names = {process_name(r) for r in trace}
    assert names <= {"genReport", "genGeneReport", "genUnmappedReport"}, names

    for gene in genes:
        report = glob.glob(str(analysis_dir / "**" / f"*.report_{gene}.html"), recursive=True)
        assert len(report) == 1, (gene, report)


@GEMMA_GATE
def test_stage7_alone_reports_every_gene_protein(full_run_protein):
    """The protein copy of the B6 wiring (genReport -> genProteinReport) has no prior coverage
    of its own (round 1, guards P4): it is a separate copy of the channel wiring, not shared
    code with the nucleotide branch."""
    analysis_dir = analysis_dir_of(full_run_protein)
    genes = read_top_genes(analysis_dir, n=NTOPGENES_BASE)
    assert len(genes) == NTOPGENES_BASE, genes
    for gene in genes:
        report = glob.glob(str(analysis_dir / "**" / f"*.report_{gene}.html"), recursive=True)
        assert len(report) == 1, (gene, report)

    for html in glob.glob(str(analysis_dir / "**" / "*.report_*.html"), recursive=True):
        os.remove(html)
    add_params(full_run_protein, ["\tskip1 = true", "\tskip2 = true", "\tskip3 = true", "\tskip4 = true",
                                  "\tskip5 = true", "\tskip6 = true", "\tskip7 = false", "\toverwrite = true"])
    rc, out = run_nf(full_run_protein, timeout=300, extra_args=("-with-trace", "trace3p.txt"))
    assert rc == 0, out[-5000:]

    trace = read_trace(full_run_protein / "trace3p.txt")
    names = {process_name(r) for r in trace}
    assert names <= {"genReport", "genProteinReport", "genUnmappedReport"}, names
    for gene in genes:
        report = glob.glob(str(analysis_dir / "**" / f"*.report_{gene}.html"), recursive=True)
        assert len(report) == 1, (gene, report)
