from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import pandas as pd
import yaml

import common  # noqa: F401


def _extract_invalid_option(stderr: str) -> str | None:
    marker = "Invalid command-line option "
    for line in stderr.splitlines():
        if marker not in line:
            continue
        return line.split(marker, 1)[1].strip()
    return None


def _run_gtalign_with_compat_fallback(cmd: list[str]) -> None:
    current_cmd = list(cmd)
    removed_options: set[str] = set()

    while True:
        completed = subprocess.run(current_cmd, capture_output=True, text=True)
        if completed.returncode == 0:
            return

        invalid_option = _extract_invalid_option(completed.stderr)
        if invalid_option and invalid_option not in removed_options:
            current_cmd = [
                arg for arg in current_cmd if arg != invalid_option and not arg.startswith(f"{invalid_option}=")
            ]
            removed_options.add(invalid_option)
            continue

        raise subprocess.CalledProcessError(
            completed.returncode,
            current_cmd,
            output=completed.stdout,
            stderr=completed.stderr,
        )


def _normalize_hits(payload) -> list[dict]:
    if isinstance(payload, dict):
        if "hits" in payload:
            return [payload]
        if "results" in payload:
            return payload["results"]
        return [payload]
    if isinstance(payload, list):
        return payload
    return []


def _first_present(mapping: dict, keys: list[str], default=None):
    for key in keys:
        if key in mapping:
            return mapping[key]
    return default


def parse_gtalign_json(path: Path) -> list[dict]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if "gtalign_search" in payload:
        search = payload["gtalign_search"]
        query = search.get("query", {})
        query_name = Path(str(_first_present(query, ["description"], ""))).stem or None
        query_length = _first_present(query, ["length"])
        rows: list[dict] = []
        for record in search.get("search_results", []):
            hit = record.get("hit_record", {})
            alignment = hit.get("alignment", {})
            reference_name = Path(str(_first_present(hit, ["reference_description"], ""))).stem or None
            rows.append(
                {
                    "query": query_name,
                    "reference": reference_name,
                    "aligned_residues": int(_first_present(alignment, ["n_aligned", "n_matched"], 0)),
                    "rmsd": float(_first_present(alignment, ["rmsd"], 0.0)),
                    "tm_score_query": float(_first_present(alignment, ["tmscore_query"], 0.0)),
                    "tm_score_reference": float(_first_present(alignment, ["tmscore_refn"], 0.0)),
                    "query_length_reported": query_length,
                    "reference_length_reported": _first_present(hit, ["reference_length"]),
                }
            )
        return rows

    rows: list[dict] = []
    for query_entry in _normalize_hits(payload):
        query_name = _first_present(
            query_entry,
            ["query_name", "query", "query_description", "query_id"],
        )
        query_length = _first_present(query_entry, ["query_length", "query_len"])
        for hit in query_entry.get("hits", []):
            rows.append(
                {
                    "query": Path(str(query_name)).stem if query_name else None,
                    "reference": Path(str(_first_present(hit, ["hit_name", "reference_name", "reference", "description"]))).stem,
                    "aligned_residues": int(_first_present(hit, ["aligned_residues", "aligned_positions", "aligned_length"], 0)),
                    "rmsd": float(_first_present(hit, ["RMSD", "rmsd"], 0.0)),
                    "tm_score_query": float(_first_present(hit, ["TM_score_query_normalized", "tm_score_query", "TM2"], 0.0)),
                    "tm_score_reference": float(_first_present(hit, ["TM_score_ref_normalized", "tm_score_reference", "TM1"], 0.0)),
                    "query_length_reported": query_length,
                    "reference_length_reported": _first_present(hit, ["reference_length", "hit_length", "target_length"]),
                }
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    manifest = pd.read_csv(args.manifest, sep="\t")
    structures_dir = Path(config["structures_dir"])
    gtalign_bin = config["structure"]["gtalign_bin"]
    gtalign_args = config["structure"].get("gtalign_args", [])

    output_dir = args.output.parent / "gtalign_raw"
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        gtalign_bin,
        f"--qrs={structures_dir}",
        f"--rfs={structures_dir}",
        "-o",
        str(output_dir),
        *gtalign_args,
    ]
    _run_gtalign_with_compat_fallback(cmd)

    rows: list[dict] = []
    for json_path in sorted(output_dir.glob("*.json")):
        rows.extend(parse_gtalign_json(json_path))

    if not rows:
        raise RuntimeError("GTalign finished but no JSON results were parsed.")

    df = pd.DataFrame(rows)
    samples = set(manifest["sample"].tolist())
    df = df[df["query"].isin(samples) & df["reference"].isin(samples)].copy()
    df["pair_key"] = df.apply(lambda row: tuple(sorted((row["query"], row["reference"]))), axis=1)
    # Keep the best-scoring observation for each unordered pair.
    df = (
        df.sort_values(
            by=["pair_key", "tm_score_reference", "tm_score_query", "aligned_residues"],
            ascending=[True, False, False, False],
        )
        .drop_duplicates(subset=["pair_key"])
        .drop(columns=["pair_key"])
        .sort_values(["query", "reference"])
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
