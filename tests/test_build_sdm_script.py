import importlib.util
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / "workflow" / "scripts"
SCRIPT_PATH = SCRIPT_DIR / "build_sdm.py"

if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

spec = importlib.util.spec_from_file_location("build_sdm", SCRIPT_PATH)
module = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(module)


def test_build_sdm_script_writes_zero_diagonal(tmp_path, monkeypatch):
    pairwise_path = tmp_path / "pairwise.tsv"
    manifest_path = tmp_path / "manifest.tsv"
    matrix_output = tmp_path / "sdm.tsv"
    phylip_output = tmp_path / "sdm.phy"

    pd.DataFrame(
        [
            {"sample": "A", "residue_count": 100},
            {"sample": "B", "residue_count": 110},
        ]
    ).to_csv(manifest_path, sep="\t", index=False)
    pd.DataFrame(
        [
            {"query": "A", "reference": "B", "topologically_equivalent_residues": 90, "aligned_residues": 95, "rmsd": 1.2},
            {"query": "B", "reference": "A", "topologically_equivalent_residues": 90, "aligned_residues": 95, "rmsd": 1.2},
        ]
    ).to_csv(pairwise_path, sep="\t", index=False)

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "build_sdm.py",
            "--pairwise",
            str(pairwise_path),
            "--manifest",
            str(manifest_path),
            "--matrix-output",
            str(matrix_output),
            "--phylip-output",
            str(phylip_output),
            "--rmsd-cap",
            "3.0",
            "--min-score",
            "1.0e-6",
        ],
    )

    module.main()

    matrix = pd.read_csv(matrix_output, sep="\t", index_col=0)
    assert matrix.loc["A", "A"] == 0.0
    assert matrix.loc["B", "B"] == 0.0
    assert matrix.loc["A", "B"] == matrix.loc["B", "A"]

    phylip_lines = phylip_output.read_text(encoding="utf-8").splitlines()
    assert phylip_lines[0] == "2"
