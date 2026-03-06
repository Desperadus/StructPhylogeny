from __future__ import annotations

import argparse
from pathlib import Path

import common  # noqa: F401
from structphylogeny.io import load_records, write_fasta


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--structures-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    records = load_records(args.structures_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(records, args.output)


if __name__ == "__main__":
    main()
