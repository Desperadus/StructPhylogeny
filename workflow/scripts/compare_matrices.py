from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from scipy.stats import pearsonr, spearmanr

import common  # noqa: F401
from structphylogeny.matrices import read_square_matrix


def flatten_upper_triangle(matrix: pd.DataFrame) -> pd.Series:
    values = []
    labels = []
    for i, row_name in enumerate(matrix.index):
        for j, col_name in enumerate(matrix.columns):
            if j <= i:
                continue
            labels.append(f"{row_name}|{col_name}")
            values.append(matrix.iloc[i, j])
    return pd.Series(values, index=labels)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--structure-matrix", type=Path, required=True)
    parser.add_argument("--sequence-matrix", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    structure = read_square_matrix(args.structure_matrix)
    sequence = read_square_matrix(args.sequence_matrix)
    shared = sorted(set(structure.index) & set(sequence.index))
    structure = structure.loc[shared, shared]
    sequence = sequence.loc[shared, shared]

    structure_flat = flatten_upper_triangle(structure)
    sequence_flat = flatten_upper_triangle(sequence)
    pearson = pearsonr(structure_flat, sequence_flat)
    spearman = spearmanr(structure_flat, sequence_flat)

    summary = pd.DataFrame(
        [
            {"metric": "n_taxa", "value": len(shared)},
            {"metric": "n_pairs", "value": len(structure_flat)},
            {"metric": "pearson_r", "value": pearson.statistic},
            {"metric": "pearson_pvalue", "value": pearson.pvalue},
            {"metric": "spearman_rho", "value": spearman.statistic},
            {"metric": "spearman_pvalue", "value": spearman.pvalue},
        ]
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
