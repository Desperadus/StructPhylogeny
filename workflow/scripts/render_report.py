from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

import common  # noqa: F401
from structphylogeny.trees import load_newick_terminals


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--sequence-tree", type=Path, required=True)
    parser.add_argument("--structure-tree", type=Path, required=True)
    parser.add_argument("--correlation", type=Path, required=True)
    parser.add_argument("--sdm", type=Path, required=True)
    parser.add_argument("--sequence-distance", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, sep="\t")
    correlation = pd.read_csv(args.correlation, sep="\t").set_index("metric")["value"].to_dict()
    structure_taxa = load_newick_terminals(args.structure_tree)
    sequence_taxa = load_newick_terminals(args.sequence_tree)

    lines = [
        "# StructPhylogeny Report",
        "",
        "## Dataset",
        "",
        f"- Structures processed: {len(manifest)}",
        f"- Mean residue count: {manifest['residue_count'].mean():.1f}",
        f"- Sequence tree taxa: {len(sequence_taxa)}",
        f"- Structure tree taxa: {len(structure_taxa)}",
        "",
        "## Matrix agreement",
        "",
        f"- Pearson correlation between structural SDM and sequence p-distance: {float(correlation['pearson_r']):.4f}",
        f"- Spearman correlation between structural SDM and sequence p-distance: {float(correlation['spearman_rho']):.4f}",
        f"- Pair count compared: {int(float(correlation['n_pairs']))}",
        "",
        "## Key outputs",
        "",
        f"- Manifest: `{args.manifest}`",
        f"- Sequence tree: `{args.sequence_tree}`",
        f"- Structure tree: `{args.structure_tree}`",
        f"- Structural SDM matrix: `{args.sdm}`",
        f"- Sequence distance matrix: `{args.sequence_distance}`",
    ]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
