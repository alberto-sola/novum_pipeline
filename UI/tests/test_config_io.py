from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # UI/ on path

import config_io


def test_roundtrip_preserves_structure(tmp_path):
    cfg = {
        "queries": {"ecoli": "a.fa"}, "targets": {"ecoli": "g.fna"},
        "threads": 16, "max_suboptimal_hits": 1, "max_total_energy": -20, "seed": None,
        "intarna": {"accessibility_variant": "on",
                    "output": {"columns": "id1,id2,E", "overlap": "B"}},
        "plots": {"type": None}, "results_dir": "Data/Results/",
    }
    p = tmp_path / "config.yaml"
    assert config_io.save_config(p, cfg)["ok"] is True
    assert config_io.load_config(p) == cfg

def test_special_characters_survive_roundtrip(tmp_path):
    # The #10 regression: ':' '#' and a leading '-' must not corrupt the YAML.
    cfg = {"results_dir": "a: b # c", "note": "-leading", "variant": "on"}
    p = tmp_path / "config.yaml"
    config_io.save_config(p, cfg)
    assert config_io.load_config(p) == cfg

def test_missing_file_returns_empty(tmp_path):
    assert config_io.load_config(tmp_path / "nope.yaml") == {}

def test_bad_yaml_returns_error(tmp_path):
    p = tmp_path / "config.yaml"
    p.write_text("queries: [unterminated\n")
    assert "__error__" in config_io.load_config(p)

def test_save_is_atomic_no_partial(tmp_path):
    p = tmp_path / "config.yaml"
    config_io.save_config(p, {"a": 1})
    # no leftover temp files in the directory
    assert [x.name for x in tmp_path.iterdir()] == ["config.yaml"]
