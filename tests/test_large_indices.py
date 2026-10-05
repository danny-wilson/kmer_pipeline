"""Counts, batch boundaries and pattern indices of 100,000 and more are written as plain
integers (R wrote 100000 as "1e+05", which head -n and tail -n reject), and results written
that way by earlier releases can still be read."""
import gzip
import os

import pytest

from conftest import EXAMPLE_DIR, run_script

import kmercontigalignmerge
import Manhattan_functions
import rcompat
import stringlist2patternandkinship


def test_parse_index():
    assert rcompat.parse_index("100000") == 100000
    assert rcompat.parse_index("1e+05") == 100000
    assert rcompat.parse_index("2e+05") == 200000
    with pytest.raises(ValueError):
        rcompat.parse_index("1.5")


def test_kmer_batch_numbers():
    batches = stringlist2patternandkinship.get_kmer_batch_numbers(n=200000, p=2)
    assert batches == [(1, 100000), (100001, 200000)]
    beg, end = batches[0]
    assert rcompat.r_paste0("pre.", beg, "-", end) == "pre.1-100000"
    assert rcompat.r_paste0(end - beg + 1) == "100000"


def test_count_batch_parameters():
    assert kmercontigalignmerge.get_count_batch_parameters(p=2, b=100000, n=200000) == \
        [(1, 100000), (100001, 200000)]
    assert kmercontigalignmerge.format_count(100000.0) == "100000"
    assert kmercontigalignmerge.format_count(1000000.0) == "1000000"


GEMMA_WRAPPER = """#!/bin/bash
# keep a copy of GEMMA's genotype file (-g), then run GEMMA
args=("$@")
for ((k=0; k<${{#args[@]}}; k++)); do
  [ "${{args[$k]}}" = -g ] && cp "${{args[$((k+1))]}}" {capture}/
done
exec /usr/local/bin/gemma "$@"
"""


def software_file(root, gemma):
    rows = open(os.path.join(EXAMPLE_DIR, "pipeline_software_location.txt")).read().split("\n")
    rows = [f"gemma\t{gemma}" if r.startswith("gemma\t") else r for r in rows]
    (root / "software.txt").write_text("\n".join(rows))
    return str(root / "software.txt")


@pytest.mark.skipif(not os.path.exists("/usr/local/bin/gemma"), reason="needs GEMMA (the kmer_pipeline image)")
def test_rungemma_200000_patterns(tmp_path):
    """Two GEMMA batches of exactly 100,000 patterns each, then read back in step 6."""
    npat, phenos = 200000, ["1.5", "-2", "3", "0.25", "7", "-1"]
    analysis = tmp_path / "analysis"
    analysis.mkdir()
    (tmp_path / "capture").mkdir()
    ids = ["id\tpaths\tpheno"] + [f"S{i}\t/dev/null\t{v}" for i, v in enumerate(phenos)]
    (tmp_path / "ids.txt").write_text("\n".join(ids) + "\n")
    pats = ["010011", "101100", "110001", "001110", "011010"]
    (analysis / "pre_protein5.patternmerge.patternKey.txt.gz").write_bytes(
        gzip.compress("".join(pats[k % len(pats)] + "\n" for k in range(npat)).encode()))
    (analysis / "pre_protein5.patternmerge.patternKeySize.txt").write_text(f"{npat}\n")
    kin = [[0.6 if i == j else 0.1 + 0.01 * ((i + j) % 3) for j in range(6)] for i in range(6)]
    (analysis / "pre_protein5.kinshipmerge.kinship.txt.gz").write_bytes(
        gzip.compress("".join(" ".join("%.17g" % v for v in r) + "\n" for r in kin).encode()))
    wrapper = tmp_path / "gemma.sh"
    wrapper.write_text(GEMMA_WRAPPER.format(capture=tmp_path / "capture"))
    wrapper.chmod(0o755)
    sw = software_file(tmp_path, wrapper)

    for t in (1, 2):
        result = run_script("rungemma.py", "--task-id", str(t), "--p", "2",
                            "--kmerfile-prefix", str(analysis / "pre_protein5"), "--id-file", str(tmp_path / "ids.txt"),
                            "--output-prefix", "pre", "--analysis-dir", str(analysis), "--kmertype", "protein",
                            "--kmer-length", "5", "--software-file", sw, cwd=tmp_path)
        assert result.returncode == 0, result.stdout[-2000:] + result.stderr[-2000:]

    geno = sorted(os.listdir(tmp_path / "capture"))
    assert geno == ["pre_protein5.1-100000.bimbam.txt", "pre_protein5.100001-200000.bimbam.txt"]
    with open(tmp_path / "capture" / geno[0]) as fh:
        rs = [line.split("\t", 1)[0] for line in fh]
    assert len(rs) == 100000 and rs[0] == "1" and rs[-1] == "100000"
    with open(tmp_path / "capture" / geno[1]) as fh:
        rs = [line.split("\t", 1)[0] for line in fh]
    assert rs[0] == "100001" and rs[99999] == "200000"

    out = str(analysis / "proteinkmer5_gemma" / "output") + "/"
    assoc = Manhattan_functions.read_gemma_files(input_dir=out, prefix="pre", kmer_type="protein", kmer_length=5,
                                                 nPatterns=npat)
    assert len(assoc) == npat
    assert assoc[99999][0] == "100000" and assoc[199999][0] == "200000"


def test_read_gemma_files_written_by_earlier_releases(tmp_path):
    """Batch file names and rs values like "1e+05" from earlier releases still match."""
    out = tmp_path / "output"
    out.mkdir()
    header = "chr\trs\tps\tn_miss\tbeta\tse\tl_remle\tl_mle\tp_wald\tp_lrt\tp_score\tlogl_H1\n"
    for beg, end, name in ((1, 100000, "1-1e+05"), (100001, 200000, "100001-2e+05")):
        rows = "".join(f"-9\t{rcompat.r_str(float(k), 7)}\t-9\t0\t0.5\t0.1\t1\t1\t0.01\t0.02\t0.03\t-10\n"
                       for k in range(beg, end + 1))
        (out / f"pre_protein5.{name}.assoc.txt.gz").write_bytes(gzip.compress((header + rows).encode()))
    log = [f"## line {k}" for k in range(20)]
    log[12], log[16] = "## lambda = 1", "## log-likelihood under the null = -12"
    (out / "pre_protein5.1-1e+05.log.txt.gz").write_bytes(gzip.compress(("\n".join(log) + "\n").encode()))
    assoc = Manhattan_functions.read_gemma_files(input_dir=str(out) + "/", prefix="pre", kmer_type="protein",
                                                 kmer_length=5, nPatterns=200000)
    assert assoc[99999][0] == "1e+05" and assoc[199999][0] == "2e+05"
    assert all(row is not None for row in assoc)


def test_read_gemma_files_with_a_nan_row(tmp_path):
    """N10: a pattern GEMMA cannot fit gives -nan; it is kept (as NaN), not a crash."""
    out = tmp_path / "output"
    out.mkdir()
    header = "chr\trs\tps\tn_miss\tbeta\tse\tl_remle\tl_mle\tp_wald\tp_lrt\tp_score\tlogl_H1\n"
    rows = ("-9\t1\t-9\t0\t0.75\t3.6\t279\t585\t0.83\t0.86\t0.87\t-75.05\n"
            "-9\t2\t-9\t0\t-nan\t-nan\t1e+05\t1e+05\tnan\t-nan\tnan\t-nan\n")
    (out / "pre_protein5.1-2.assoc.txt.gz").write_bytes(gzip.compress((header + rows).encode()))
    log = [f"## line {k}" for k in range(20)]
    log[12], log[16] = "## lambda = 1", "## log-likelihood under the null = -75.0673"
    (out / "pre_protein5.1-2.log.txt.gz").write_bytes(gzip.compress(("\n".join(log) + "\n").encode()))
    assoc = Manhattan_functions.read_gemma_files(input_dir=str(out) + "/", prefix="pre", kmer_type="protein",
                                                 kmer_length=5, nPatterns=2)
    nl = Manhattan_functions.assoc_column(assoc, 6)
    assert nl[0] > 0 and nl[1] != nl[1]
