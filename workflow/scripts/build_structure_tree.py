from __future__ import annotations

import argparse
from pathlib import Path

import common  # noqa: F401
from structphylogeny.matrices import read_square_matrix
from structphylogeny.trees import write_neighbor_joining_tree


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    matrix = read_square_matrix(args.matrix)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_neighbor_joining_tree(matrix, args.output)


if __name__ == "__main__":
    main()
