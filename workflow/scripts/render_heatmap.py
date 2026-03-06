from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import seaborn as sns

import common  # noqa: F401
from structphylogeny.matrices import read_square_matrix


def figure_size_for_matrix(size: int) -> tuple[float, float]:
    edge = max(12.0, min(30.0, size * 0.45))
    return (edge+5, edge)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--title", required=True)
    args = parser.parse_args()

    matrix = read_square_matrix(args.matrix)
    fig, ax = plt.subplots(figsize=figure_size_for_matrix(len(matrix)))
    sns.heatmap(
        matrix,
        cmap="viridis",
        ax=ax,
        xticklabels=matrix.columns.to_list(),
        yticklabels=matrix.index.to_list(),
    )
    ax.set_xticks([index + 0.5 for index in range(len(matrix.columns))])
    ax.set_yticks([index + 0.5 for index in range(len(matrix.index))])
    ax.set_xticklabels(matrix.columns.to_list(), rotation=90, fontsize=20)
    ax.set_yticklabels(matrix.index.to_list(), rotation=0, fontsize=20)
    ax.set_title(args.title)
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=200)


if __name__ == "__main__":
    main()
