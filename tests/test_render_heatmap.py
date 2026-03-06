import importlib.util
import sys
from pathlib import Path

import matplotlib
import pandas as pd


matplotlib.use("Agg")

ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / "workflow" / "scripts"
SCRIPT_PATH = SCRIPT_DIR / "render_heatmap.py"

if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

spec = importlib.util.spec_from_file_location("render_heatmap", SCRIPT_PATH)
module = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(module)


def test_render_heatmap_sets_all_axis_labels(tmp_path, monkeypatch):
    matrix_path = tmp_path / "matrix.tsv"
    output_path = tmp_path / "heatmap.png"
    names = ["sample_a", "sample_b", "sample_c"]
    pd.DataFrame(
        [
            [0.0, 1.0, 2.0],
            [1.0, 0.0, 3.0],
            [2.0, 3.0, 0.0],
        ],
        index=names,
        columns=names,
    ).to_csv(matrix_path, sep="\t")

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "render_heatmap.py",
            "--matrix",
            str(matrix_path),
            "--output",
            str(output_path),
            "--title",
            "Test Heatmap",
        ],
    )

    module.main()

    figure = module.plt.gcf()
    axis = figure.axes[0]

    assert [tick.get_text() for tick in axis.get_xticklabels()] == names
    assert [tick.get_text() for tick in axis.get_yticklabels()] == names
    assert output_path.exists()
