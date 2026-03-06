from __future__ import annotations

from io import StringIO
from pathlib import Path

import pandas as pd
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor


def matrix_to_distance_matrix(matrix: pd.DataFrame) -> DistanceMatrix:
    names = matrix.index.tolist()
    lower_triangle = []
    for i in range(len(names)):
        row = [float(matrix.iloc[i, j]) for j in range(i + 1)]
        lower_triangle.append(row)
    return DistanceMatrix(names, lower_triangle)


def write_neighbor_joining_tree(matrix: pd.DataFrame, output_path: Path) -> None:
    dm = matrix_to_distance_matrix(matrix)
    tree = DistanceTreeConstructor().nj(dm)
    with output_path.open("w", encoding="utf-8") as handle:
        Phylo.write(tree, handle, "newick")


def load_newick_terminals(path: Path) -> list[str]:
    tree = Phylo.read(str(path), "newick")
    return [clade.name for clade in tree.get_terminals()]


def canonical_newick(path: Path) -> str:
    tree = Phylo.read(str(path), "newick")
    buffer = StringIO()
    Phylo.write(tree, buffer, "newick")
    return buffer.getvalue().strip()
