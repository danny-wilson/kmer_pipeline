"""preflight.py: overwrite checks, out-of-date outputs, reuse checks, -resume and parameters."""
import json
import os
import socket

import pytest

from conftest import run_script

import inventory

PREFIX, TYPE, K = "tb20", "nucleotide", 31
P = f"{PREFIX}_{TYPE}{K}"
GENOMES = ["702", "725", "791"]

# One file per step, as a run leaves them
STEP_FILES = {
    1: [f"{TYPE}kmer{K}/702.kmer{K}.txt.gz", f"{P}_kmers_filepaths.txt"],
    2: [f"{P}.kmermerge.txt.gz"],
    3: [f"{P}.patternmerge.patternKey.txt.gz", f"{TYPE}kmer{K}_patternbatches/{P}.1-9.patternKey.txt.gz"],
    4: [f"{TYPE}kmer{K}_gemma/output/{P}.1-9.assoc.txt.gz"],
    5: [f"{TYPE}kmer{K}_kmergenealign/{P}_NC_000962.3_kmergenecombination_filepaths.txt",
        f"{P}.NC_000962.3_t90.kmeralignmerge.txt.gz"],
    6: [f"{P}.summary.json", f"{TYPE}kmer{K}_kmergenealign_figures/figure_data/params.tsv"],
    7: [f"{P}.report.html", f"{P}.report_rpoB.html", "report.css"],
}


def make_run(analysis, steps=range(1, 8), genomes=GENOMES):
    for s in steps:
        for f in STEP_FILES[s]:
            path = analysis / f
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("x\n")
    if 1 in steps:
        (analysis / f"{P}_kmers_filepaths.txt").write_text(
            "".join(f"/home/jovyan/a/{TYPE}kmer{K}/{g}.kmer{K}.txt.gz\n" for g in genomes))


def id_file(tmp_path, genomes=GENOMES):
    path = tmp_path / "id_file.txt"
    path.write_text("id\tpaths\tpheno\n" + "".join(f"{g}\t/dev/null\t{k + 1}\n" for k, g in enumerate(genomes)))
    return str(path)


def preflight(tmp_path, analysis, run_steps=range(1, 8), overwrite=False, resume=False, session="s1",
              user_params=(), params=None, inputs=None, ids=None, finish=None, covs=""):
    args = ["--analysis-dir", str(analysis), "--output-prefix", PREFIX, "--kmer-type", TYPE, "--kmer-length", str(K),
            "--session-id", session, "--run-steps", ",".join(map(str, run_steps)),
            "--overwrite", str(overwrite).lower(), "--resume", str(resume).lower(), "--pid", "999999999",
            "--user-params", ",".join(user_params), "--params-json", json.dumps(params or {"kmer_min_count": 1}),
            "--input-files", json.dumps(inputs or {}), "--id-file", ids or id_file(tmp_path), "--covariate-file", covs]
    if finish:
        args += ["--finish", finish]
    result = run_script("preflight.py", *args)
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


def remaining(analysis):
    return sorted(os.path.relpath(os.path.join(r, n), analysis) for r, _, ns in os.walk(analysis) for n in ns)


def test_fresh_analysis_dir(tmp_path):
    analysis = tmp_path / "kmergwas"
    out = preflight(tmp_path, analysis)
    assert out["errors"] == [] and out["deleted"] == []
    manifest = json.loads((analysis / f"{P}.run_manifest.json").read_text())
    assert manifest["status"] == "running" and sorted(manifest["steps"]) == [str(s) for s in range(1, 8)]


def test_rerun_stops_by_default(tmp_path):
    analysis = tmp_path / "kmergwas"
    make_run(analysis)
    before = remaining(analysis)
    out = preflight(tmp_path, analysis)
    assert len(out["errors"]) == 1 and "overwrite = true" in out["errors"][0]
    assert remaining(analysis) == before


def test_overwrite_deletes_the_steps_that_run(tmp_path):
    analysis = tmp_path / "kmergwas"
    make_run(analysis)
    out = preflight(tmp_path, analysis, overwrite=True)
    assert out["errors"] == []
    assert sorted(out["deleted"]) == sorted(f for s in range(1, 8) for f in STEP_FILES[s])
    left = remaining(analysis)
    assert left[-1] == f"{P}.run_manifest.json" and len(left) == 2
    assert left[0].startswith(f"log.{P}/overwrite_deleted_")
    assert sorted((analysis / left[0]).read_text().split()) == sorted(out["deleted"])


def test_phenotype_rerun_keeps_the_alignments(tmp_path):
    """P108: steps 4, 6 and 7 rerun; steps 1-3 and 5 are kept (step 5 does not read step 4)."""
    analysis = tmp_path / "kmergwas"
    make_run(analysis)
    out = preflight(tmp_path, analysis, run_steps=[4, 6, 7], overwrite=True)
    assert out["errors"] == []
    assert sorted(out["deleted"]) == sorted(STEP_FILES[4] + STEP_FILES[6] + STEP_FILES[7])
    for s in (1, 2, 3, 5):
        for f in STEP_FILES[s]:
            assert (analysis / f).exists()


def test_out_of_date_outputs_of_skipped_steps_are_deleted(tmp_path):
    """Q5: rerunning step 6 alone also deletes step 7's reports, which it makes out of date."""
    analysis = tmp_path / "kmergwas"
    make_run(analysis)
    out = preflight(tmp_path, analysis, run_steps=[6], overwrite=True)
    assert sorted(out["deleted"]) == sorted(STEP_FILES[6] + STEP_FILES[7])


def test_reading_out_of_date_outputs_is_an_error(tmp_path):
    """Rerunning step 3 but skipping step 4, while step 6 runs, would use stale GEMMA results."""
    analysis = tmp_path / "kmergwas"
    make_run(analysis)
    out = preflight(tmp_path, analysis, run_steps=[3, 6, 7], overwrite=True)
    assert any("step 6 would read the outputs of step 4" in e for e in out["errors"])
    assert out["deleted"] == [] and (analysis / STEP_FILES[3][0]).exists()


def test_symbolic_links_are_removed_not_followed(tmp_path):
    analysis = tmp_path / "kmergwas"
    make_run(analysis, steps=[1, 2, 3])
    target = tmp_path / "elsewhere"
    target.mkdir()
    (target / "keep.txt").write_text("keep\n")
    os.symlink(target, analysis / f"{TYPE}kmer{K}_gemma")
    out = preflight(tmp_path, analysis, run_steps=[4], overwrite=True)
    assert out["errors"] == []
    assert not os.path.lexists(analysis / f"{TYPE}kmer{K}_gemma") and (target / "keep.txt").exists()


def test_other_analyses_files_are_left_alone(tmp_path):
    analysis = tmp_path / "kmergwas"
    make_run(analysis)
    other = analysis / "other_protein11.kmermerge.txt.gz"
    other.write_text("x\n")
    preflight(tmp_path, analysis, overwrite=True)
    assert other.exists()


def test_genomes_must_match_step1(tmp_path):
    analysis = tmp_path / "kmergwas"
    make_run(analysis, steps=[1, 2, 3])
    out = preflight(tmp_path, analysis, run_steps=[4, 5, 6, 7], ids=id_file(tmp_path, ["725", "702", "791"]))
    assert any("their order" in e for e in out["errors"])
    out = preflight(tmp_path, analysis, run_steps=[4, 5, 6, 7], ids=id_file(tmp_path, GENOMES))
    assert not any("their order" in e for e in out["errors"])


def test_reused_step3_made_with_another_kmer_min_count(tmp_path):
    analysis = tmp_path / "kmergwas"
    preflight(tmp_path, analysis, params={"kmer_min_count": 1})  # manifest of a full run
    make_run(analysis)
    preflight(tmp_path, analysis, finish="finished")
    out = preflight(tmp_path, analysis, run_steps=[4, 6, 7], overwrite=True, params={"kmer_min_count": 2})
    assert any("kmer_min_count" in e for e in out["errors"])


def test_unknown_parameters(tmp_path):
    out = preflight(tmp_path, tmp_path / "a", user_params=["kmerMinCount", "kmer-min-count", "queue", "maxp"])
    assert sum("did you mean 'kmer_min_count'" in e for e in out["errors"]) == 2
    assert any("'queue' is not used" in w for w in out["warnings"])


def test_the_kebab_case_alias_nextflow_adds_for_a_known_parameter_is_not_unknown(tmp_path):
    out = preflight(tmp_path, tmp_path / "a", user_params=["annotateGeneFile", "annotate-gene-file"])
    assert out["errors"] == [] and out["warnings"] == []
    # an alias of a name that is not a parameter is still reported
    out = preflight(tmp_path, tmp_path / "b", user_params=["kmerMinCount", "kmer-min-count"])
    assert sum("did you mean 'kmer_min_count'" in e for e in out["errors"]) == 2


def test_resume(tmp_path):
    analysis = tmp_path / "kmergwas"
    preflight(tmp_path, analysis, session="s1")
    make_run(analysis)
    assert preflight(tmp_path, analysis, resume=True, session="s1")["errors"] == []
    assert "another run" in preflight(tmp_path, analysis, resume=True, session="s2")["errors"][0]
    changed = preflight(tmp_path, analysis, resume=True, session="s1", params={"kmer_min_count": 3})
    assert "kmer_min_count" in changed["errors"][0]
    assert "cannot be used with -resume" in preflight(tmp_path, analysis, resume=True, overwrite=True)["errors"][0]
    assert all((analysis / f).exists() for s in range(1, 8) for f in STEP_FILES[s])


def test_a_live_run_blocks_overwrite(tmp_path):
    analysis = tmp_path / "kmergwas"
    make_run(analysis)
    (analysis / f"{P}.run_manifest.json").write_text(json.dumps(
        {"status": "running", "host": socket.gethostname(), "pid": os.getpid(), "session": "s0", "steps": {}}))
    out = preflight(tmp_path, analysis, overwrite=True)
    assert "may still be running" in out["errors"][0] and out["deleted"] == []


def test_finish(tmp_path):
    analysis = tmp_path / "kmergwas"
    preflight(tmp_path, analysis, session="s1")
    preflight(tmp_path, analysis, session="s1", finish="finished")
    assert json.loads((analysis / f"{P}.run_manifest.json").read_text())["status"] == "finished"


def test_every_workflow_parameter_is_known():
    import preflight as pf
    from test_parameters import NF
    if NF is None:
        pytest.skip("needs kmer_pipeline.nf")
    import re
    used = set(re.findall(r"params\.([A-Za-z_0-9]+)", open(NF).read())) - {"containsKey", "keySet", "container"}
    assert used <= pf.KNOWN_PARAMS, used - pf.KNOWN_PARAMS


def test_inventory_downstream():
    assert inventory.downstream([4]) == {6, 7}
    assert inventory.downstream([2]) == {3, 4, 5, 6, 7}
    assert inventory.downstream([5]) == {6, 7}


def raw_id_file(tmp_path, rows, header="id\tpaths\tpheno"):
    path = tmp_path / "ids_raw.txt"
    path.write_text(header + "\n" + "".join("\t".join(r) + "\n" for r in rows))
    return str(path)


def genomes(phenos, ids=None):
    ids = ids or [f"G{k}" for k in range(len(phenos))]
    return [[i, "/dev/null", p] for i, p in zip(ids, phenos)]


CONT = ["1.2", "-0.5", "3", "2.2", "0.1", "NA", "4", "1"]


def check(tmp_path, rows, run_steps=range(1, 8), header="id\tpaths\tpheno", covs=""):
    out = preflight(tmp_path, tmp_path / "a", run_steps=run_steps, ids=raw_id_file(tmp_path, rows, header), covs=covs)
    return out["errors"], out["warnings"]


def test_valid_inputs(tmp_path):
    assert check(tmp_path, genomes(CONT)) == ([], [])
    tf = ["TRUE", "FALSE"] * 12  # read as 1/0, as the scripts do
    assert check(tmp_path, genomes(tf))[0] == []


def test_id_file_header(tmp_path):
    errors, _ = check(tmp_path, genomes(CONT), header="ID\tpaths\tpheno")
    assert "lower case" in errors[0]


def test_ids(tmp_path):
    errors, _ = check(tmp_path, genomes(CONT, ["a", "b", "a", "c/d", "e", "f", "g", "h"]))
    assert any("duplicate IDs" in e for e in errors) and any("cannot contain /" in e for e in errors)
    errors, _ = check(tmp_path, genomes(CONT, ["1.1", "1.10", "3", "4", "5", "6", "7", "8"]))
    assert any("makes them equal" in e for e in errors)
    errors, warnings = check(tmp_path, genomes(CONT, ["007", "8", "9", "10", "11", "12", "13", "14"]))
    assert errors == [] and any("007 -> 7" in w for w in warnings)


def test_phenotypes(tmp_path):
    errors, _ = check(tmp_path, genomes(["1", "R", "1,5", "2", "3", "4", "NA", ""]))
    assert any("'R'" in e and "'1,5'" in e for e in errors)
    errors, warnings = check(tmp_path, genomes(["1", "NaN", "2", "3", "4", "5"]))  # worked before: missing
    assert errors == [] and any("NaN" in w for w in warnings)
    errors, _ = check(tmp_path, genomes(["1", "Inf", "2", "3", "4", "5"]))
    assert any("infinite" in e for e in errors)
    _, warnings = check(tmp_path, genomes(["1", "-9", "2", "3", "4", "5"]))
    assert any("-9" in w for w in warnings)
    errors, _ = check(tmp_path, genomes(["2"] * 6 + ["NA"]))
    assert any("same phenotype" in e for e in errors)


def test_phenotypes_not_checked_for_steps_1_to_3(tmp_path):
    """Steps 1-3 (and 5) do not use the phenotype: all NA is fine (precomputing)."""
    assert check(tmp_path, genomes(["NA"] * 8), run_steps=[1, 2, 3, 5])[0] == []
    assert check(tmp_path, genomes(["NA"] * 8), run_steps=[1, 2, 3, 4])[0] != []


def test_small_binary_group(tmp_path):
    _, warnings = check(tmp_path, genomes(["1"] * 3 + ["0"] * 20))
    assert any("only 3 genomes in the smaller group" in w for w in warnings)


def test_covariates(tmp_path):
    def cov(rows):
        path = tmp_path / "cov.txt"
        path.write_text("".join("\t".join(r) + "\n" for r in rows))
        return str(path)
    rows = genomes(CONT)
    assert check(tmp_path, rows, covs=cov([["1", str(k * 0.3 % 1)] for k in range(8)]))[0] == []
    errors, _ = check(tmp_path, rows, covs=cov([["1", "0.5"]] * 7))
    assert "7 rows but id_file has 8" in errors[0]
    errors, _ = check(tmp_path, rows, covs=cov([["2", str(k)] for k in range(8)]))
    assert "first column" in errors[0]
    errors, _ = check(tmp_path, rows, covs=cov([["1", "0.5"]] * 8))
    assert any("linearly dependent" in e for e in errors)


def test_precomputed_dir(tmp_path):
    pre = tmp_path / "pre"
    make_run(pre, steps=[1, 2, 3, 5])
    for f in (".patternmerge.patternKeySize.txt", ".patternmerge.patternIndex.txt.gz", ".kinshipmerge.kinship.txt.gz",
              ".NC_000962.3_t90.kmeralignmerge.count.txt.gz"):
        (pre / (P + f)).write_text("x\n")
    (pre / f"{TYPE}kmer{K}_kmergenealign/{P}_NC_000962.3_gene_id_name_lookup.txt").write_text("x\n")
    before = remaining(pre)

    def run(run_steps=(4, 6, 7), pre_dir=pre, analysis=tmp_path / "new", ids=None):
        args = ["--analysis-dir", str(analysis), "--output-prefix", PREFIX, "--kmer-type", TYPE,
                "--kmer-length", str(K), "--session-id", "s9", "--run-steps", ",".join(map(str, run_steps)),
                "--params-json", json.dumps({"kmer_min_count": "1", "nucmerident": "90"}),
                "--id-file", ids or id_file(tmp_path), "--precomputed-dir", str(pre_dir),
                "--precomputed-prefix", PREFIX, "--ref-name", "NC_000962.3"]
        return json.loads(run_script("preflight.py", *args).stdout)

    out = run()
    assert out["errors"] == [], out
    assert remaining(pre) == before  # precomputed_dir is only read
    assert any("cannot also run" in e for e in run(run_steps=[3, 4, 6, 7])["errors"])
    assert any("must differ" in e for e in run(analysis=pre)["errors"])
    assert any("their order" in e for e in run(ids=id_file(tmp_path, ["791", "725", "702"]))["errors"])
    os.remove(pre / (P + ".kinshipmerge.kinship.txt.gz"))
    assert any("lacks outputs" in e and "kinship" in e for e in run()["errors"])


def test_pheno_file_checks(tmp_path):
    rows = genomes(CONT)
    pf = tmp_path / "pheno.txt"
    pf.write_text("id\tpheno\nG0\t1\nG1\t2\nG2\tR\nG3\t4\nZZ\t5\n")
    out = preflight(tmp_path, tmp_path / "a", ids=raw_id_file(tmp_path, rows))
    args_pheno = ["--pheno-file", str(pf)]
    res = run_script("preflight.py", "--analysis-dir", str(tmp_path / "b"), "--output-prefix", PREFIX,
                     "--kmer-type", TYPE, "--kmer-length", str(K), "--session-id", "s", "--run-steps", "4,6,7",
                     "--id-file", raw_id_file(tmp_path, rows), *args_pheno)
    out = json.loads(res.stdout)
    assert any("pheno_file, ID G2: 'R'" in e for e in out["errors"])
    assert any("1 IDs of pheno_file are not in id_file" in w for w in out["warnings"])
    assert any("4 genomes of id_file have no phenotype" in w for w in out["warnings"])


def test_resume_figure_only_parameters_may_change(tmp_path):
    analysis = tmp_path / "kmergwas"
    base = {"kmer_min_count": 1, "ntopgenes": "20", "blastident": "70"}
    preflight(tmp_path, analysis, session="s1")
    make_run(analysis)
    out = preflight(tmp_path, analysis, resume=True, session="s1", params={**base, "ntopgenes": "5", "blastident": "80"})
    assert out["errors"] == []


def test_min_contig_length_must_be_a_whole_number(tmp_path):
    out = preflight(tmp_path, tmp_path / "a", params={"kmer_min_count": 1, "min_contig_length": "-5"})
    assert any("min_contig_length" in e for e in out["errors"])
    assert preflight(tmp_path, tmp_path / "b", params={"kmer_min_count": 1, "min_contig_length": "310"})["errors"] == []


def test_resume_refuses_a_changed_min_contig_length(tmp_path):
    analysis = tmp_path / "kmergwas"
    preflight(tmp_path, analysis, session="s1", params={"kmer_min_count": 1, "min_contig_length": "0"})
    make_run(analysis)
    out = preflight(tmp_path, analysis, resume=True, session="s1", params={"kmer_min_count": 1, "min_contig_length": "310"})
    assert "min_contig_length" in out["errors"][0]
