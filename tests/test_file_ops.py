"""rcompat's file operations that replace shell commands (D8)."""
import gzip
import os

import pytest

import rcompat


def test_size_touch_remove(tmp_path):
    f = tmp_path / "a.txt"
    assert rcompat.file_size_lines(str(f)) == []
    rcompat.touch(str(f))
    assert rcompat.file_size_lines(str(f)) == ["0"]
    f.write_text("hello\n")
    assert rcompat.file_size_lines(str(f)) == ["6"]
    rcompat.remove(str(f))
    assert not f.exists()
    with pytest.raises(RuntimeError, match="rm"):
        rcompat.remove(str(f))


def test_gzip_file(tmp_path):
    f = tmp_path / "a.txt"
    f.write_text("x\ny\n")
    rcompat.gzip_file(str(f))
    assert not f.exists() and gzip.open(str(f) + ".gz", "rt").read() == "x\ny\n"
    f.write_text("again\n")  # a rerun task (-resume) replaces its own earlier output
    rcompat.gzip_file(str(f))
    assert not f.exists() and gzip.open(str(f) + ".gz", "rt").read() == "again\n"


def test_counts(tmp_path):
    f = tmp_path / "a.txt"
    f.write_text("1\n2\n3")  # wc -l counts newlines
    assert rcompat.count_lines(str(f)) == 2
    g = tmp_path / "b.txt.gz"
    g.write_bytes(gzip.compress(b"1\n2\n3\n"))
    assert rcompat.count_lines(str(g)) == 3 and rcompat.count_bytes(str(g)) == 6
    assert rcompat.count_bytes(str(tmp_path / "missing.gz")) == 0


def test_cut_and_write(tmp_path):
    g = tmp_path / "t.txt.gz"
    g.write_bytes(gzip.compress(b"a\tb\tc\n1\t2\t3\n"))
    assert list(rcompat.cut_fields(str(g), (1, 3))) == ["a\tc", "1\t3"]
    out = str(tmp_path / "o.txt.gz")
    rcompat.write_gz_lines(out, rcompat.cut_fields(str(g), (2,)))
    assert gzip.open(out, "rt").read() == "b\n2\n"


def test_move_copy_ls(tmp_path):
    (tmp_path / "x.1.done").write_text("1")
    rcompat.copy(str(tmp_path / "x.1.done"), str(tmp_path / "x.2.done"))
    rcompat.move(str(tmp_path / "x.2.done"), str(tmp_path / "x.3.done"))
    assert [os.path.basename(p) for p in rcompat.ls(str(tmp_path / "x.*.done"))] == ["x.1.done", "x.3.done"]
    assert rcompat.ls(str(tmp_path / "none*")) == []
