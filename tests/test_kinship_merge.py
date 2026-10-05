"""The kinship merge waits for batch files that appear late (D1b)."""
import gzip

import numpy as np

import stringlist2patternandkinship as s2p


def write_batch(prefix, kinship, weight):
    with gzip.open(prefix + ".kinship.txt.gz", "wt") as fh:
        fh.write("".join(" ".join("%.17g" % v for v in row) + "\n" for row in kinship))
    with open(prefix + ".kinshipWeight.txt", "w") as fh:
        fh.write("%d\n" % weight)
    open(prefix + ".kinship.completed.txt", "w").close()


def test_waits_for_a_late_batch(tmp_path, monkeypatch):
    batches_dir = str(tmp_path) + "/"
    nkmersbatch = [(1, 10), (11, 30)]
    prefixes = [batches_dir + f"pre_nucleotide31.{b}-{e}" for b, e in nkmersbatch]
    files = [p + ".kinship.txt.gz" for p in prefixes]
    k1, k2 = np.array([[1.0, 0.25], [0.25, 0.5]]), np.array([[0.5, 0.0], [0.0, 1.0]])
    write_batch(prefixes[0], k1, 10)
    sleeps = []

    def late_sleep(seconds):  # the second batch arrives while the merge waits
        sleeps.append(seconds)
        write_batch(prefixes[1], k2, 20)

    monkeypatch.setattr(s2p.time, "sleep", late_sleep)
    s2p.merge_kinship_matrices(t=1, n=2, b=2, p=2, imax=1, files=files, prefix="pre", nkmersbatch=nkmersbatch,
                               output_dir=batches_dir, kmertype="nucleotide", kmerlen=31, batches_dir=batches_dir)
    assert sleeps == [60]
    out = batches_dir + "pre_nucleotide31.kinshipmerge.j.1.1"
    merged = s2p.read_kinship(out + ".kinship.txt.gz", 2)
    np.testing.assert_allclose(merged, (10 * k1 + 20 * k2) / 30)
    assert open(out + ".kinshipWeight.txt").read() == "30\n"


def test_gives_up_after_merge_wait_minutes(tmp_path, monkeypatch):
    """N7: a batch that never arrives stops the merge after merge_wait_minutes, saying so."""
    import pytest
    import rcompat

    batches_dir = str(tmp_path) + "/"
    nkmersbatch = [(1, 10), (11, 30)]
    prefixes = [batches_dir + f"pre_nucleotide31.{b}-{e}" for b, e in nkmersbatch]
    write_batch(prefixes[0], np.eye(2), 10)
    sleeps = []
    monkeypatch.setattr(s2p.time, "sleep", lambda s: sleeps.append(s))
    monkeypatch.setattr(rcompat, "MERGE_WAIT_MINUTES", 100.0)
    rcompat.set_merge_wait("3")
    with pytest.raises(rcompat.RError, match="merge_wait_minutes"):
        s2p.merge_kinship_matrices(t=1, n=2, b=2, p=2, imax=1, files=[p + ".kinship.txt.gz" for p in prefixes],
                                   prefix="pre", nkmersbatch=nkmersbatch, output_dir=batches_dir,
                                   kmertype="nucleotide", kmerlen=31, batches_dir=batches_dir)
    assert sum(sleeps) == 3 * 60


def test_merge_wait_minutes_must_be_positive():
    import pytest
    import rcompat
    for bad in ("0", "-5", "soon"):
        with pytest.raises(rcompat.RError, match="positive"):
            rcompat.set_merge_wait(bad)


def test_merge_scripts_accept_merge_wait_minutes():
    from conftest import run_script
    for script in ("createfullkmerlist.py", "nucleotidekmermerge.py", "proteinkmermerge.py",
                   "stringlist2patternandkinship.py", "kmercontigalignmerge.py", "kmercontigalign.py"):
        assert "--merge-wait-minutes" in run_script(script, "--help").stdout, script


def test_alignment_only_has_no_merge_wait():
    from conftest import run_script
    assert "--merge-wait-minutes" not in run_script("kmercontigalignonly.py", "--help").stdout
