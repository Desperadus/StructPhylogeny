import importlib.util
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / "workflow" / "scripts"
SCRIPT_PATH = SCRIPT_DIR / "run_gtalign_all_vs_all.py"

if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

spec = importlib.util.spec_from_file_location("run_gtalign_all_vs_all", SCRIPT_PATH)
module = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(module)


def test_extract_invalid_option_parses_gtalign_stderr():
    stderr = (
        "[gtalign] ERROR: Invalid command-line option --nalns\n"
        "    (File: src/libutil/CLOptions.cpp; Line: 44; Function: AssignCLOptO_NALNS)\n"
    )
    assert module._extract_invalid_option(stderr) == "--nalns"


def test_extract_invalid_option_returns_none_for_other_errors():
    stderr = "some unrelated error"
    assert module._extract_invalid_option(stderr) is None


def test_parse_gtalign_json_supports_gtalign_search_schema(tmp_path):
    payload = {
        "gtalign_search": {
            "query": {"description": "data/example/FABP9_model_0_minimized.cif Chn:A0 (M:1)", "length": 142},
            "search_results": [
                {
                    "hit_record": {
                        "reference_description": "data/example/FABP4_model_0_minimized.cif Chn:A0 (M:1)",
                        "reference_length": 132,
                        "alignment": {
                            "n_aligned": 132,
                            "rmsd": 0.64,
                            "tmscore_query": 0.9121,
                            "tmscore_refn": 0.9799,
                        },
                    }
                }
            ],
        }
    }
    path = tmp_path / "gtalign.json"
    path.write_text(json.dumps(payload), encoding="utf-8")

    assert module.parse_gtalign_json(path) == [
        {
            "query": "FABP9_model_0_minimized",
            "reference": "FABP4_model_0_minimized",
            "aligned_residues": 132,
            "rmsd": 0.64,
            "tm_score_query": 0.9121,
            "tm_score_reference": 0.9799,
            "query_length_reported": 142,
            "reference_length_reported": 132,
        }
    ]


def test_parse_gtalign_json_supports_weak_hits_when_present(tmp_path):
    payload = {
        "gtalign_search": {
            "query": {"description": "data/example/FABP6_model_0_minimized.cif Chn:A0 (M:1)", "length": 128},
            "search_results": [
                {
                    "hit_record": {
                        "reference_description": "data/example/LCN16_model_0_minimized.cif Chn:A0 (M:1)",
                        "reference_length": 161,
                        "alignment": {
                            "n_aligned": 163,
                            "rmsd": 4.35,
                            "tmscore_query": 0.49817,
                            "tmscore_refn": 0.42247,
                        },
                    }
                }
            ],
        }
    }
    path = tmp_path / "gtalign.json"
    path.write_text(json.dumps(payload), encoding="utf-8")

    assert module.parse_gtalign_json(path) == [
        {
            "query": "FABP6_model_0_minimized",
            "reference": "LCN16_model_0_minimized",
            "aligned_residues": 163,
            "rmsd": 4.35,
            "tm_score_query": 0.49817,
            "tm_score_reference": 0.42247,
            "query_length_reported": 128,
            "reference_length_reported": 161,
        }
    ]
