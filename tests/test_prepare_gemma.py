"""prepare_gemma.py (step 4's preparation): the analysed genomes, GEMMA's phenotype file, presence
counts and the decompressed kinship matrix; and how steps 4 and 6 use them."""
import gzip
import json
import math
import os

import pytest

from conftest import run_script

import Manhattan_functions
import prepare_gemma as pg
import rcompat


def test_analysed_set():
    pheno = [1.0, None, math.nan, 2.0, 3.0, math.inf]
    assert pg.analysed_set(pheno, []) == [True, False, False, True, True, False]
    cov = [[1, 0.5], [1, 0.1], [1, 0.2], [1, None], [1, 0.3], [1, 0.4]]
    assert pg.analysed_set(pheno, cov) == [True, False, False, False, True, False]


def test_check_analysed():
    assert pg.check_analysed([1.0, 2.0, 3.0], [True] * 3, []) == []
    assert "only 2 genomes" in pg.check_analysed([1.0, 2.0, None], [True, True, False], [])[0]
    assert "same phenotype" in pg.check_analysed([1.0, 1.0, 1.0, 1.0], [True] * 4, [])[0]
    pheno = [1.0, 2.0, 3.0, 4.0, 5.0]
    constant = [[1, 0.5]] * 5  # a covariate with one value: collinear with the intercept
    assert "linearly dependent" in pg.check_analysed(pheno, [True] * 5, constant)[0]
    near = [[1, z, round(2 * z, 3)] for z in (0.1234, 0.5, 0.25, 0.9, 0.6)]  # GEMMA fits these
    assert pg.check_analysed(pheno, [True] * 5, near) == []


def test_analysed_file_round_trip(tmp_path):
    path = str(tmp_path / "a.txt")
    pg.write_analysed(path, ["s1", "s2", "s3"], [1.5, None, math.nan], [True, False, False])
    assert pg.read_analysed(path) == (["s1", "s2", "s3"], [1.5, None, None])


def make_step3(analysis, patterns, n):
    analysis.mkdir(parents=True, exist_ok=True)
    pre = analysis / "pre_protein5"
    (analysis / "pre_protein5.patternmerge.patternKey.txt.gz").write_bytes(
        gzip.compress("".join(p + "\n" for p in patterns).encode()))
    (analysis / "pre_protein5.patternmerge.patternKeySize.txt").write_text(f"{len(patterns)}\n")
    kin = "".join(" ".join("%.17g" % (0.6 if i == j else 0.1) for j in range(n)) + "\n" for i in range(n))
    (analysis / "pre_protein5.kinshipmerge.kinship.txt.gz").write_bytes(gzip.compress(kin.encode()))
    return str(pre), kin


def prepare(tmp_path, phenos, covariates=None, cleanup=False):
    ids = tmp_path / "ids.txt"
    ids.write_text("id\tpaths\tpheno\n" + "".join(f"S{i}\t/dev/null\t{v}\n" for i, v in enumerate(phenos)))
    args = ["--kmerfile-prefix", str(tmp_path / "analysis" / "pre_protein5"), "--id-file", str(ids),
            "--analysis-dir", str(tmp_path / "analysis"), "--output-prefix", "pre", "--kmer-type", "protein",
            "--kmer-length", "5"]
    if covariates:
        (tmp_path / "cov.txt").write_text("".join("\t".join(r) + "\n" for r in covariates))
        args += ["--covariate-file", str(tmp_path / "cov.txt")]
    if cleanup:
        args.append("--cleanup")
    return run_script("prepare_gemma.py", *args)


def test_prepare_gemma(tmp_path):
    patterns = ["11001011", "01111000", "10000111"]
    _, kin = make_step3(tmp_path / "analysis", patterns, 8)
    phenos = ["1.5", "NA", "3", "NaN", "0.123456789012", "2", "-1", "4"]
    covariates = [["1", v] for v in ("0.1", "0.2", "NA", "0.4", "0.5", "0.7", "0.3", "0.9")]
    result = prepare(tmp_path, phenos, covariates)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "genome S2 has a missing covariate" in result.stdout
    gdir = tmp_path / "analysis" / "proteinkmer5_gemma"
    # Analysed: S0, S4, S5, S6, S7 (S1 NA, S2 covariate NA, S3 NaN)
    ids, pheno = pg.read_analysed(str(gdir / "pre_protein5.analysed_phenotypes.txt"))
    assert pheno == [1.5, None, None, None, 0.123456789012, 2.0, -1.0, 4.0]
    sent = (gdir / "pre_protein5_gemma_formatted_phenotype.txt").read_text().split("\n")[:-1]
    assert sent[1:4] == ["NA", "NA", "NA"] and float(sent[4]) == 0.123456789012
    with gzip.open(tmp_path / "analysis" / "pre_protein5.patternmerge.presenceCount.txt.gz", "rt") as fh:
        counts = fh.read().split()
    assert counts == ["4", "1", "4"]  # presence in columns 0, 4, 5, 6, 7 of each pattern
    assert (gdir / "pre_protein5.kinship.txt").read_text() == kin
    assert prepare(tmp_path, phenos, covariates, cleanup=True).returncode == 0
    assert not (gdir / "pre_protein5.kinship.txt").exists()


def test_prepare_gemma_stops_when_the_model_cannot_be_fitted(tmp_path):
    make_step3(tmp_path / "analysis", ["110010"], 6)
    result = prepare(tmp_path, ["1", "1", "1", "NA", "1", "1"])
    assert result.returncode != 0 and "same phenotype" in result.stderr


def test_summary_json_is_valid_with_missing_values(tmp_path):
    path = str(tmp_path / "s.json")
    Manhattan_functions.write_summary_json(path, n_kmers=10, n_patterns=4, n_untested_patterns=1,
                                           max_neglog10p=math.nan, minor_allele_threshold=0.01, macormaf="maf",
                                           n_tests=2, bonferroni=1.6, n_genomes=6, n_genomes_analysed=4,
                                           n_patterns_nan=1, pheno_type="binary", pheno=[0, 1, None, 1, 1, None])
    s = json.load(open(path))
    assert s["max_neglog10p"] is None and s["n_genomes_analysed"] == 4
    assert (s["n_cases"], s["n_controls"], s["case_value"], s["control_value"]) == (3, 1, 1, 0)


def gemma_log(path, n_individuals, n_snps, lognull="-75.06"):
    lines = [f"## line {k}" for k in range(20)]
    lines[5] = f"## number of analyzed individuals = {n_individuals}"
    lines[7] = f"## number of analyzed SNPs = {n_snps}"
    lines[12], lines[16] = "## lambda = 1", f"## log-likelihood under the null = {lognull}"
    with gzip.open(path, "wt") as fh:
        fh.write("\n".join(lines) + "\n")


def test_check_gemma_logs(tmp_path):
    log = str(tmp_path / "x.1-5.log.txt.gz")
    gemma_log(log, 28, 5)
    Manhattan_functions.check_gemma_logs([log], n_analysed=28, n_rows=5)
    with pytest.raises(rcompat.RError, match="GEMMA analysed 28 genomes"):
        Manhattan_functions.check_gemma_logs([log], n_analysed=29, n_rows=5)
    with pytest.raises(rcompat.RError, match="truncated"):
        Manhattan_functions.check_gemma_logs([log], n_analysed=28, n_rows=4)
    gemma_log(log, 28, 5, lognull="-nan")
    with pytest.raises(rcompat.RError, match="null model"):
        Manhattan_functions.check_gemma_logs([log], n_analysed=28, n_rows=5)


def test_analysed_phenotypes_for_an_earlier_release(tmp_path):
    """No analysed-phenotype file: rebuilt from id_file and checked against the GEMMA logs."""
    ids = tmp_path / "ids.txt"
    ids.write_text("id\tpaths\tpheno\nS0\t/dev/null\t1\nS1\t/dev/null\tNaN\nS2\t/dev/null\t2\n")
    log = str(tmp_path / "pre_protein5.1-5.log.txt.gz")
    gemma_log(log, 2, 5)
    pheno = pg.analysed_phenotypes(str(tmp_path), "pre", "protein", 5, str(ids), None, [log])
    assert pheno == [1.0, None, 2.0]
    gemma_log(log, 3, 5)
    with pytest.raises(rcompat.RError, match="not the inputs step 4 was run with"):
        pg.analysed_phenotypes(str(tmp_path), "pre", "protein", 5, str(ids), None, [log])


def test_pheno_file_is_matched_by_id_as_text(tmp_path):
    """N5, P123: 007 in a mixed id_file must match 007 in an all-numeric pheno_file, not 7."""
    ids = tmp_path / "ids.txt"
    ids.write_text("id\tpaths\tpheno\n007\t/x\t1\nX8\t/x\t2\n7\t/x\t3\nG4\t/x\t4\n")
    pheno_file = tmp_path / "pheno.txt"
    pheno_file.write_text("id\tpheno\n7\t70\n007\t0.07\n99\t9\n")
    report = {}
    got_ids, pheno = pg.read_phenotypes(str(ids), str(pheno_file), report)
    assert pheno == [0.07, None, 70.0, None]
    assert report == {"matched": ["007", "7"], "missing": ["X8", "G4"], "unknown": ["99"]}


def test_covariates_with_an_id_column(tmp_path):
    ids = tmp_path / "ids.txt"
    ids.write_text("id\tpaths\tpheno\nA\t/x\t1\nB\t/x\t2\nC\t/x\t3\n")
    cov = tmp_path / "cov.txt"
    cov.write_text("id\tintercept\tage\nC\t1\t30\nA\t1\t10\nZ\t1\t5\n")
    report = {}
    rows = pg.read_covariates(str(cov), str(ids), report)
    assert rows == [[1.0, 10.0], [None, None], [1.0, 30.0]]
    assert report == {"missing": ["B"], "unknown": ["Z"]}
    assert pg.analysed_set([1.0, 2.0, 3.0], rows) == [True, False, True]
    out = str(tmp_path / "gemma_cov.txt")
    pg.write_gemma_covariates(out, rows)
    assert open(out).read() == "1\t10\n1\tNA\n1\t30\n"
    plain = tmp_path / "plain.txt"
    plain.write_text("1\t10\n1\t20\n1\t30\n")
    assert pg.read_covariates(str(plain), str(ids)) == [[1.0, 10.0], [1.0, 20.0], [1.0, 30.0]]
