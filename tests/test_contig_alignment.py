"""Contig alignment helpers in kmercontigalignonly.py."""
import kmercontigalignonly as kca


def test_covered_bases_pairs_each_alignment_start_with_its_end():
    """D1c: two passing alignments cover 1-10 and 21-30 of the contig (a reverse one too)."""
    alignfile = [[100, 109, 1, 10], [500, 509, 30, 21], [900, 909, 41, 50]]
    covered = kca.covered_contig_bases(alignfile, passing=[0, 1])
    assert covered == set(range(1, 11)) | set(range(21, 31))
    assert kca.covered_contig_bases(alignfile, passing=[]) == set()


def test_kmer_gene_pairs_when_every_kmer_has_two_genes():
    """D1d: every k-mer with the same number (2) of genes keeps its k-mer in each pair."""
    genes, names = kca.kmer_gene_pairs({"AAC": [3, 4], "GGT": [7, 8]})
    assert genes == [3.0, 4.0, 7.0, 8.0]
    assert names == ["AAC", "AAC", "GGT", "GGT"]


def test_kmer_gene_pairs_mixed():
    genes, names = kca.kmer_gene_pairs({"AAC": [3], "GGT": [7, 8]})
    assert list(zip(names, genes)) == [("AAC", 3.0), ("GGT", 7.0), ("GGT", 8.0)]


def test_protein_kmer_windows():
    assert kca.protein_kmer_windows("ACDEF", 3, 3) == [[1, 2, 3], [2, 3, 4], [3, 4, 5]]


def test_protein_kmer_windows_frame_shorter_than_k():
    """D1e: a frame shorter than k has no k-mers and no windows."""
    assert kca.count_protein_kmers("ACD", 5) is None
    assert kca.protein_kmer_windows("ACD", 0, 5) == []
