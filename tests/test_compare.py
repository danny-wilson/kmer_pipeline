"""tests/e2e/compare.py needs no goldens and no image, so it lives under
tests/ (not tests/e2e/, which is collect-ignored without an e2e config) and
runs in a plain `pytest tests`. Covers the 'path' class (B9): content is
compared after each side's own base_dir prefix is replaced by "@BASE_DIR@"."""
import os
import sys

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
E2E_DIR = os.path.join(os.path.dirname(TESTS_DIR), "tests", "e2e")
if E2E_DIR not in sys.path:
    sys.path.insert(0, E2E_DIR)

import compare  # noqa: E402

MANIFEST = "*.txt\tpath\n"


def write(root, rel, text):
    path = os.path.join(root, rel)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def make_manifest(tmp_path):
    path = tmp_path / "manifest.tsv"
    path.write_text(MANIFEST)
    return str(path)


def test_differing_prefixes_alone_do_not_cause_a_problem(tmp_path):
    manifest = make_manifest(tmp_path)
    dir_a, dir_b = str(tmp_path / "a"), str(tmp_path / "b")
    write(dir_a, "f.txt", "root=/run/a/base dir\nsame line\n")
    write(dir_b, "f.txt", "root=/elsewhere/b/base dir\nsame line\n")
    problems, counts, total = compare.compare(dir_a, dir_b, manifest,
                                               prefix_a="/run/a/base", prefix_b="/elsewhere/b/base")
    assert problems == []
    assert counts == {"path": 1}


def test_a_perturbed_byte_outside_the_prefix_still_fails(tmp_path):
    manifest = make_manifest(tmp_path)
    dir_a, dir_b = str(tmp_path / "a"), str(tmp_path / "b")
    write(dir_a, "f.txt", "root=/run/a/base dir\nsame line\n")
    write(dir_b, "f.txt", "root=/elsewhere/b/base dir\nDIFFERENT line\n")
    problems, counts, total = compare.compare(dir_a, dir_b, manifest,
                                               prefix_a="/run/a/base", prefix_b="/elsewhere/b/base")
    assert any("f.txt" in p for p in problems)


def test_with_no_prefix_given_the_raw_content_is_compared(tmp_path):
    manifest = make_manifest(tmp_path)
    dir_a, dir_b = str(tmp_path / "a"), str(tmp_path / "b")
    write(dir_a, "f.txt", "root=/run/a/base dir\n")
    write(dir_b, "f.txt", "root=/elsewhere/b/base dir\n")
    problems, counts, total = compare.compare(dir_a, dir_b, manifest)
    assert any("f.txt" in p for p in problems)


def test_double_slashes_in_the_prefixed_text_are_preserved_not_collapsed(tmp_path):
    manifest = make_manifest(tmp_path)
    dir_a, dir_b = str(tmp_path / "a"), str(tmp_path / "b")
    write(dir_a, "f.txt", "/run/a/base/kmergwas//out.txt\n")
    write(dir_b, "f.txt", "/elsewhere/b/base/kmergwas//out.txt\n")
    problems, counts, total = compare.compare(dir_a, dir_b, manifest,
                                               prefix_a="/run/a/base", prefix_b="/elsewhere/b/base")
    assert problems == []
