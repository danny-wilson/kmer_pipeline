"""read_ref_name()'s pure-Groovy ref_name_from_fasta() (PLAN.md Part B, B3): a port of the
now-deleted get_ref_name.py, verified once against it directly (DECISIONS.md) across every case
below. Extracts the live function straight out of kmer_pipeline.nf (nf_block, as
test_parameters.py already does for min_count), so this test tracks the real implementation
rather than a hand-maintained duplicate that could silently drift from it."""
import gzip
import os

import pytest

from conftest import EXAMPLE_DIR
from test_parameters import NEXTFLOW_RUNS, NF, nf_block


def _write_probe(tmp_path):
    block = nf_block("def ref_name_from_fasta", "}")
    probe = tmp_path / "probe.nf"
    probe.write_text(block + """
try {
	println "RESULT:" + ref_name_from_fasta(params.test_fasta_path)
} catch (Exception e) {
	System.err.println(e.getMessage())
	System.exit(1)
}
""")
    return probe


def _result_line(stdout):
    line = next(l for l in stdout.split("\n") if l.startswith("RESULT:"))
    return line[len("RESULT:"):]


def _run(tmp_path, fasta_path):
    import subprocess
    probe = _write_probe(tmp_path)
    env = dict(os.environ, NXF_OPTS="-Dnxf.ansi.log=false")
    result = subprocess.run(["nextflow", "run", str(probe), "--test_fasta_path", str(fasta_path)],
                            cwd=tmp_path, env=env, capture_output=True, text=True)
    return result


pytestmark = pytest.mark.skipif(not NEXTFLOW_RUNS, reason="needs kmer_pipeline.nf and a working Nextflow")


@pytest.mark.parametrize("content, expected", [
    (b">chr1 description here\nACGT\nACGT\n", "chr1"),
    (b">chr1 first record\nACGT\n>chr2 second\nTTTT\n", "chr1"),  # multi-record
    (b">chr1 crlf\r\nACGT\r\n", "chr1"),  # CRLF
    (b">chr1\tdescription\nACGT\n", "chr1\tdescription"),  # tab-delimited header
    (b">\nACGT\n", ""),  # bare > header
    (b">chr1 \xff\xfe invalid\nACGT\n", "chr1"),  # non-UTF-8 bytes after the name
])
def test_matches_python_on_success(tmp_path, content, expected):
    fasta = tmp_path / "ref.fa"
    fasta.write_bytes(content)
    result = _run(tmp_path, fasta)
    assert result.returncode == 0, result.stderr
    assert _result_line(result.stdout) == expected


@pytest.mark.parametrize("content", [
    b"\n>chr1 after blank\nACGT\n",  # blank first line
    b" >chr1 leading space\nACGT\n",  # leading space
    b"\xef\xbb\xbf>chr1 bom\nACGT\n",  # UTF-8 BOM
    b"",  # empty file
    b"no header here\nACGT\n",
])
def test_matches_python_on_error(tmp_path, content):
    fasta = tmp_path / "ref.fa"
    fasta.write_bytes(content)
    result = _run(tmp_path, fasta)
    assert result.returncode != 0


def test_the_real_example_reference(tmp_path):
    """The exact check the deleted get_ref_name.py smoke test made, against the real example
    reference used throughout example/ and the e2e goldens."""
    fasta = os.path.join(EXAMPLE_DIR, "Mtub_H37Rv_NC000962.3.fasta")
    result = _run(tmp_path, fasta)
    assert result.returncode == 0, result.stderr
    assert _result_line(result.stdout) == "NC_000962.3"


def test_gzip_is_decompressed(tmp_path):
    fasta = tmp_path / "ref.fa.gz"
    fasta.write_bytes(gzip.compress(b">chr1 gzipped\nACGT\n"))
    result = _run(tmp_path, fasta)
    assert result.returncode == 0, result.stderr
    assert _result_line(result.stdout) == "chr1"


@pytest.mark.parametrize("magic, label", [
    (b"BZh91AY&SY", "bzip2"),
    (b"\xfd7zXZ\x00", "xz"),
])
def test_bzip2_and_xz_are_explicitly_unsupported(tmp_path, magic, label):
    """Accepted divergence from get_ref_name.py (DECISIONS.md): no bzip2/xz decoder is available
    on Nextflow's classpath, so these now fail with a clear, named error instead of decompressing."""
    fasta = tmp_path / f"ref.fa.{label}"
    fasta.write_bytes(magic + b"padding-so-the-file-is-long-enough-to-not-look-truncated")
    result = _run(tmp_path, fasta)
    assert result.returncode != 0
    assert label in result.stderr
