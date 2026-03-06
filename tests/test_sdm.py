import math

from structphylogeny.sdm import compute_sdm


def test_compute_sdm_self_like_pair_is_small():
    score = compute_sdm(query_length=100, reference_length=100, rmsd=0.0, topologically_equivalent_residues=100)
    assert math.isclose(score, 0.0, abs_tol=1e-8)


def test_compute_sdm_increases_when_alignment_is_worse():
    better = compute_sdm(query_length=150, reference_length=150, rmsd=1.0, topologically_equivalent_residues=120)
    worse = compute_sdm(query_length=150, reference_length=150, rmsd=2.5, topologically_equivalent_residues=80)
    assert worse > better


def test_compute_sdm_uses_min_score_when_rmsd_exceeds_cap():
    score = compute_sdm(query_length=128, reference_length=161, rmsd=3.36, topologically_equivalent_residues=142)
    assert math.isclose(score, 69.31471805599453)


def test_compute_sdm_prefers_topologically_equivalent_residues():
    score = compute_sdm(
        query_length=100,
        reference_length=120,
        rmsd=1.0,
        topologically_equivalent_residues=60,
        aligned_residues=90,
    )
    expected = compute_sdm(
        query_length=100,
        reference_length=120,
        rmsd=1.0,
        topologically_equivalent_residues=60,
    )
    assert math.isclose(score, expected)
