import math

from structphylogeny.sdm import compute_sdm


def test_compute_sdm_self_like_pair_is_small():
    score = compute_sdm(aligned_residues=100, query_length=100, reference_length=100, rmsd=0.0)
    assert math.isclose(score, 0.0, abs_tol=1e-8)


def test_compute_sdm_increases_when_alignment_is_worse():
    better = compute_sdm(aligned_residues=120, query_length=150, reference_length=150, rmsd=1.0)
    worse = compute_sdm(aligned_residues=80, query_length=150, reference_length=150, rmsd=2.5)
    assert worse > better


def test_compute_sdm_uses_min_score_when_rmsd_exceeds_cap():
    score = compute_sdm(aligned_residues=142, query_length=128, reference_length=161, rmsd=3.36)
    assert math.isclose(score, 69.31471805599453)
