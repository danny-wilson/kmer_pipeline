"""countkmers.filter_short_contigs: contigs below --min-contig-length are left out of the k-mer counts."""
import gzip
import os

import pytest

from conftest import SCRIPTS_DIR  # noqa: F401  (puts the scripts on sys.path)

import countkmers

FASTA = ">long1 note\n" + "ACGT" * 30 + "\nACGT\n>short\n" + "ACGT" * 5 + "\n>long2\n" + "TTGA" * 25 + "\n"


def headers(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        return [line.strip() for line in fh if line.startswith(">")]


def test_zero_keeps_everything_and_copies_nothing(tmp_path):
    src = tmp_path / "c.fa"
    src.write_text(FASTA)
    assert countkmers.filter_short_contigs(str(src), 0, str(tmp_path / "out.fa")) == str(src)
    assert not (tmp_path / "out.fa").exists()


def test_short_contigs_are_dropped_and_sequence_lines_joined(tmp_path):
    src = tmp_path / "c.fa"
    src.write_text(FASTA)
    out = countkmers.filter_short_contigs(str(src), 100, str(tmp_path / "out.fa"))
    assert headers(out) == [">long1 note", ">long2"]
    assert "ACGT" * 30 + "ACGT" in open(out).read().replace("\n", "")  # a wrapped contig keeps all its lines


def test_threshold_is_inclusive(tmp_path):
    src = tmp_path / "c.fa"
    src.write_text(FASTA)
    assert ">short" in headers(countkmers.filter_short_contigs(str(src), 20, str(tmp_path / "out.fa")))
    assert ">short" not in headers(countkmers.filter_short_contigs(str(src), 21, str(tmp_path / "out2.fa")))


def test_gzipped_input_and_crlf(tmp_path):
    src = tmp_path / "c.fa.gz"
    with gzip.open(src, "wt", newline="") as fh:
        fh.write(FASTA.replace("\n", "\r\n"))
    out = countkmers.filter_short_contigs(str(src), 100, str(tmp_path / "out.fa"))
    assert headers(out) == [">long1 note", ">long2"]
    assert "\r" not in open(out).read()


def test_returned_path_is_absolute_so_a_chdir_does_not_lose_it(tmp_path, monkeypatch):
    src = tmp_path / "c.fa"
    src.write_text(FASTA)
    monkeypatch.chdir(tmp_path)
    out = countkmers.filter_short_contigs(str(src), 100, "relative_out.fa")
    assert os.path.isabs(out)
    os.chdir("/")
    assert os.path.exists(out)
