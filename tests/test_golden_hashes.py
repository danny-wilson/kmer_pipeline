"""tests/e2e/golden_hashes.py needs no goldens and no image, so it lives under
tests/ (not tests/e2e/, which is collect-ignored without an e2e config) and
runs in a plain `pytest tests`."""
import gzip
import os
import sys

import pytest

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
E2E_DIR = os.path.join(os.path.dirname(TESTS_DIR), "tests", "e2e")
if E2E_DIR not in sys.path:
    sys.path.insert(0, E2E_DIR)

import golden_hashes  # noqa: E402

MANIFEST = """\
*.ignoreme\tignore
*.html\tstrip\t^STAMP:
*.png\tvisual
*\texact
"""


@pytest.fixture
def manifest(tmp_path):
    path = tmp_path / "manifest.tsv"
    path.write_text(MANIFEST)
    return str(path)


def write(path, data, mode="w"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, mode) as fh:
        fh.write(data)


def write_gz(path, text, mtime=0):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=mtime) as fh:
        fh.write(text.encode())


def make_run(root, extra=()):
    run = str(root / "run")
    write(os.path.join(run, "a.txt"), "hello\n")
    write(os.path.join(run, "page.html"), "STAMP: Mon Jan 1 00:00:00 2026.\nbody\n")
    write(os.path.join(run, "fig.png"), "\x89PNG-bytes", mode="w")
    write(os.path.join(run, "noise.ignoreme"), "whatever\n")
    write_gz(os.path.join(run, "data.txt.gz"), "content\n")
    for rel, content in extra:
        write(os.path.join(run, rel), content)
    return run


def test_gzip_files_with_different_mtimes_hash_equal(tmp_path, manifest):
    run_a = make_run(tmp_path / "a")
    run_b = str(tmp_path / "b" / "run")
    import shutil
    shutil.copytree(run_a, run_b)
    # Rewrite data.txt.gz in run_b with a different mtime but the same text
    write_gz(os.path.join(run_b, "data.txt.gz"), "content\n", mtime=12345)
    rules = golden_hashes.read_manifest(manifest)
    ha = golden_hashes.hash_run(run_a, rules)
    hb = golden_hashes.hash_run(run_b, rules)
    assert ha["data.txt.gz"] == hb["data.txt.gz"]


def test_strip_class_lines_are_removed_before_hashing(tmp_path, manifest):
    run_a = make_run(tmp_path / "a")
    run_b = str(tmp_path / "b" / "run")
    import shutil
    shutil.copytree(run_a, run_b)
    write(os.path.join(run_b, "page.html"), "STAMP: Tue Jan 2 00:00:00 2026.\nbody\n")
    rules = golden_hashes.read_manifest(manifest)
    ha = golden_hashes.hash_run(run_a, rules)
    hb = golden_hashes.hash_run(run_b, rules)
    assert ha["page.html"] == hb["page.html"]


def test_ignore_files_are_absent_from_the_hashes(tmp_path, manifest):
    run = make_run(tmp_path)
    rules = golden_hashes.read_manifest(manifest)
    hashed = golden_hashes.hash_run(run, rules)
    assert "noise.ignoreme" not in hashed


def test_an_unclassified_file_is_an_error(tmp_path, manifest):
    run = make_run(tmp_path, extra=[("sub/weird", "x")])
    # A manifest with no catch-all leaves "weird" unclassified
    strict = tmp_path / "strict_manifest.tsv"
    strict.write_text("*.txt\texact\n")
    rules = golden_hashes.read_manifest(str(strict))
    with pytest.raises(SystemExit):
        golden_hashes.hash_run(run, rules)


def test_a_single_changed_byte_gives_a_mismatch(tmp_path, manifest):
    run = make_run(tmp_path)
    write(os.path.join(run, "a.txt"), "hellx\n")
    rules = golden_hashes.read_manifest(manifest)
    golden_hashes.generate([("run", run)], manifest)  # smoke-check it doesn't crash
    tsv = tmp_path / "golden_hashes.tsv"
    import io
    import contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        golden_hashes.generate([("run", run)], manifest)
    tsv.write_text(buf.getvalue())
    write(os.path.join(run, "a.txt"), "CHANGED\n")
    problems, warnings = golden_hashes.verify(str(tsv), [("run", run)], manifest, visual_advisory=False)
    assert any("a.txt" in p for p in problems)


def test_png_mismatch_fails_by_default_and_only_warns_with_visual_advisory(tmp_path, manifest):
    run = make_run(tmp_path)
    import io
    import contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        golden_hashes.generate([("run", run)], manifest)
    tsv = tmp_path / "golden_hashes.tsv"
    tsv.write_text(buf.getvalue())
    write(os.path.join(run, "fig.png"), "different-bytes")

    problems, warnings = golden_hashes.verify(str(tsv), [("run", run)], manifest, visual_advisory=False)
    assert any("fig.png" in p for p in problems)

    problems, warnings = golden_hashes.verify(str(tsv), [("run", run)], manifest, visual_advisory=True)
    assert not any("fig.png" in p for p in problems)
    assert any("fig.png" in w for w in warnings)


def test_verify_fails_on_a_tsv_listed_path_missing_from_the_run(tmp_path, manifest):
    run = make_run(tmp_path)
    import io
    import contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        golden_hashes.generate([("run", run)], manifest)
    tsv = tmp_path / "golden_hashes.tsv"
    tsv.write_text(buf.getvalue())
    os.remove(os.path.join(run, "a.txt"))

    problems, warnings = golden_hashes.verify(str(tsv), [("run", run)], manifest, visual_advisory=False)
    assert any("a.txt" in p and "MISSING" in p for p in problems)


def test_verify_fails_on_a_classified_file_present_but_absent_from_the_tsv(tmp_path, manifest):
    run = make_run(tmp_path)
    import io
    import contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        golden_hashes.generate([("run", run)], manifest)
    tsv = tmp_path / "golden_hashes.tsv"
    tsv.write_text(buf.getvalue())
    write(os.path.join(run, "new_exact_file.txt"), "new content\n")

    problems, warnings = golden_hashes.verify(str(tsv), [("run", run)], manifest, visual_advisory=False)
    assert any("new_exact_file.txt" in p and "NOT IN TSV" in p for p in problems)


def test_a_listed_but_missing_png_still_fails_under_visual_advisory(tmp_path, manifest):
    run = make_run(tmp_path)
    import io
    import contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        golden_hashes.generate([("run", run)], manifest)
    tsv = tmp_path / "golden_hashes.tsv"
    tsv.write_text(buf.getvalue())
    os.remove(os.path.join(run, "fig.png"))

    problems, warnings = golden_hashes.verify(str(tsv), [("run", run)], manifest, visual_advisory=True)
    assert any("fig.png" in p and "MISSING" in p for p in problems)
