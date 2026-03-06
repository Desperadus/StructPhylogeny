from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

import common  # noqa: F401
from structphylogeny.io import load_records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--structures-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    records = load_records(args.structures_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(
        {
            "sample": [record.sample for record in records],
            "path": [str(record.path) for record in records],
            "chain_id": [record.chain_id for record in records],
            "residue_count": [record.residue_count for record in records],
        }
    ).sort_values("sample")
    df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
