from __future__ import annotations

import math


def compute_sdm(
    aligned_residues: int,
    query_length: int,
    reference_length: int,
    rmsd: float,
    rmsd_cap: float = 3.0,
    min_score: float = 1.0e-6,
) -> float:
    """Compute the structural dissimilarity metric used in Lakshmi et al. 2015.

    GTalign does not report DALI PFTE directly, so PFTE is approximated here as:
    aligned_residues / min(query_length, reference_length)
    """
    if aligned_residues <= 0 or query_length <= 0 or reference_length <= 0:
        raise ValueError("Aligned residues and structure lengths must be positive")

    pfte = min(1.0, aligned_residues / min(query_length, reference_length))
    srms = max(0.0, 1.0 - (rmsd / rmsd_cap))

    w1 = ((1.0 - pfte) + (1.0 - srms)) / 2.0
    w2 = (pfte + srms) / 2.0
    score = max(min_score, (w1 * pfte) + (w2 * srms))
    return -100.0 * math.log(score)

