from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import seaborn as sns

import common  # noqa: F401
from structphylogeny.matrices import read_square_matrix


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--title", required=True)
    args = parser.parse_args()

    matrix = read_square_matrix(args.matrix)
    fig, ax = plt.subplots(figsize=(12, 10))
    sns.heatmap(matrix, cmap="viridis", ax=ax)
    ax.set_title(args.title)
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=200)


if __name__ == "__main__":
    main()
