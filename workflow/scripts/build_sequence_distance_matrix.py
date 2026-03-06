from __future__ import annotations

import argparse
from pathlib import Path

import common  # noqa: F401
from structphylogeny.matrices import alignment_distance_matrix, write_square_matrix


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    matrix = alignment_distance_matrix(args.alignment)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_square_matrix(matrix, args.output)


if __name__ == "__main__":
    main()
