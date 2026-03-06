from __future__ import annotations

import math


def compute_sdm(
    query_length: int,
    reference_length: int,
    rmsd: float,
    topologically_equivalent_residues: int | None = None,
    aligned_residues: int | None = None,
    rmsd_cap: float = 3.0,
    min_score: float = 1.0e-6,
) -> float:
    """Compute the structural dissimilarity metric used in Lakshmi et al. 2015.

    The paper defines PFTE as:
    topologically equivalent residues / length of the smaller protein
    """
    equivalent_residues = topologically_equivalent_residues
    if equivalent_residues is None:
        equivalent_residues = aligned_residues

    if equivalent_residues is None or equivalent_residues <= 0 or query_length <= 0 or reference_length <= 0:
        raise ValueError("Aligned residues and structure lengths must be positive")

    pfte = min(1.0, equivalent_residues / min(query_length, reference_length))
    srms = max(0.0, 1.0 - (rmsd / rmsd_cap))

    w1 = ((1.0 - pfte) + (1.0 - srms)) / 2.0
    w2 = (pfte + srms) / 2.0
    score = max(min_score, (w1 * pfte) + (w2 * srms))
    return -100.0 * math.log(score)
