from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

import common  # noqa: F401
from structphylogeny.matrices import write_phylip_distance_matrix, write_square_matrix
from structphylogeny.sdm import compute_sdm


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--pairwise", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--matrix-output", type=Path, required=True)
    parser.add_argument("--phylip-output", type=Path, required=True)
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    pairwise = pd.read_csv(args.pairwise, sep="\t")
    manifest = pd.read_csv(args.manifest, sep="\t")
    lengths = manifest.set_index("sample")["residue_count"].to_dict()
    samples = sorted(lengths)

    matrix = pd.DataFrame(np.nan, index=samples, columns=samples, dtype=float)
    matrix.iloc[np.arange(len(samples)), np.arange(len(samples))] = 0.0

    for _, row in pairwise.iterrows():
        query = row["query"]
        reference = row["reference"]
        if query == reference:
            continue
        topologically_equivalent_residues = row.get("topologically_equivalent_residues")
        score = compute_sdm(
            query_length=int(lengths[query]),
            reference_length=int(lengths[reference]),
            rmsd=float(row["rmsd"]),
            topologically_equivalent_residues=(
                int(topologically_equivalent_residues) if pd.notna(topologically_equivalent_residues) else None
            ),
            aligned_residues=int(row["aligned_residues"]),
            rmsd_cap=float(config["structure"]["srms_rmsd_cap"]),
            min_score=float(config["structure"]["min_score"]),
        )
        matrix.loc[query, reference] = score
        matrix.loc[reference, query] = score

    missing_pairs: list[tuple[str, str]] = []
    for i, query in enumerate(samples):
        for reference in samples[i + 1 :]:
            if pd.isna(matrix.loc[query, reference]):
                missing_pairs.append((query, reference))

    if missing_pairs:
        preview = ", ".join(f"{query} vs {reference}" for query, reference in missing_pairs[:10])
        raise RuntimeError(
            "Missing GTalign pairwise results for "
            f"{len(missing_pairs)} structure pairs. "
            "This usually means GTalign filtering removed weak hits before SDM construction. "
            f"First missing pairs: {preview}"
        )

    args.matrix_output.parent.mkdir(parents=True, exist_ok=True)
    write_square_matrix(matrix, args.matrix_output)
    write_phylip_distance_matrix(matrix, args.phylip_output)


if __name__ == "__main__":
    main()
