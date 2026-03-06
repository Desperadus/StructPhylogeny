from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from Bio import AlignIO


def read_square_matrix(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=0)


def write_square_matrix(matrix: pd.DataFrame, path: Path) -> None:
    matrix.to_csv(path, sep="\t", index=True)


def write_phylip_distance_matrix(matrix: pd.DataFrame, path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f"{matrix.shape[0]}\n")
        for sample, values in matrix.iterrows():
            name = sample[:10].ljust(10)
            row = " ".join(f"{value:.6f}" for value in values.to_list())
            handle.write(f"{name} {row}\n")


def p_distance(seq_a: str, seq_b: str) -> float:
    matches = 0
    mismatches = 0
    for aa, bb in zip(seq_a, seq_b):
        if aa == "-" or bb == "-":
            continue
        matches += 1
        if aa != bb:
            mismatches += 1
    if matches == 0:
        return 1.0
    return mismatches / matches


def alignment_distance_matrix(alignment_path: Path) -> pd.DataFrame:
    alignment = AlignIO.read(str(alignment_path), "fasta")
    names = [record.id for record in alignment]
    sequences = [str(record.seq) for record in alignment]
    matrix = np.zeros((len(names), len(names)), dtype=float)
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            distance = p_distance(sequences[i], sequences[j])
            matrix[i, j] = distance
            matrix[j, i] = distance
    return pd.DataFrame(matrix, index=names, columns=names)
