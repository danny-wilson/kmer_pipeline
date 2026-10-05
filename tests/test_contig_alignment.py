"""Contig alignment helpers in kmercontigalignonly.py."""
import kmercontigalignonly as kca


def test_covered_bases_pairs_each_alignment_start_with_its_end():
    """D1c: two passing alignments cover 1-10 and 21-30 of the contig (a reverse one too)."""
    alignfile = [[100, 109, 1, 10], [500, 509, 30, 21], [900, 909, 41, 50]]
    covered = kca.covered_contig_bases(alignfile, passing=[0, 1])
    assert covered == set(range(1, 11)) | set(range(21, 31))
    assert kca.covered_contig_bases(alignfile, passing=[]) == set()
