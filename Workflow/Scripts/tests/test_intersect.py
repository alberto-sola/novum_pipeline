from __future__ import annotations
import os, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # Workflow/Scripts on path

import pandas as pd
import intersect  # importable because the trailing call is guarded


def _rnahybrid_df():
    return pd.DataFrame({
        "miRNA": ["m1"], "Gene": ["g1"],
        "gene_name": ["G1"], "locus_tag": ["L1"], "protein_name": ["P1"],
        "protein_id": ["PID1"], "Gene_length": [300],
        "Energy": [-20.0], "P_value": [0.001], "Position": [0.1],
        "miRNA_unmatches": [""], "miRNA_matches": [""],
        "Target_matches": [""], "Target_unmatches": [""],
    })


def test_intersect_tolerates_missing_intarna_metadata(tmp_path, capsys):
    # IntaRNA side trimmed to README's stated minimum → no Pu*/ED*/seed*/dp cols.
    intarna = pd.DataFrame({
        "miRNA": ["m1"], "Gene": ["g1"],
        "E": [-9.5], "Start1": [5], "End1": [12], "Start2": [1], "End2": [8],
        "Position": [0.2],
    })
    rh = tmp_path / "rh.csv"; ia = tmp_path / "ia.csv"; out = tmp_path / "out.csv"
    _rnahybrid_df().to_csv(rh, index=False)
    intarna.to_csv(ia, index=False)

    intersect.intersect_annotations(rh, ia, out)          # must NOT raise KeyError

    merged = pd.read_csv(out)
    assert list(merged.columns) == intersect.CONSENSUS_COLUMNS  # full schema preserved
    assert merged["ED1"].isna().all()                          # absent col filled NA
    assert merged.loc[0, "E"] == -9.5
    note = capsys.readouterr().err
    assert "ED1" in note and "Pu1" in note                      # loud-but-graceful note


def test_intersect_full_columns_unchanged(tmp_path):
    intarna = pd.DataFrame({
        "miRNA": ["m1"], "Gene": ["g1"], "E": [-9.5],
        "E_hybrid": [-9.5], "ED1": [0.0], "ED2": [0.0], "Pu1": [1.0], "Pu2": [1.0],
        "Start1": [5], "End1": [12], "Start2": [1], "End2": [8], "Position": [0.2],
        "subseqDP": ["&"], "hybridDP": ["&"],
        "seedStart1": [5], "seedEnd1": [10], "seedE": [-6.0],
        "seedStart2": [1], "seedEnd2": [6],
    })
    rh = tmp_path / "rh.csv"; ia = tmp_path / "ia.csv"; out = tmp_path / "out.csv"
    _rnahybrid_df().to_csv(rh, index=False)
    intarna.to_csv(ia, index=False)
    intersect.intersect_annotations(rh, ia, out)
    merged = pd.read_csv(out)
    assert merged.loc[0, "ED1"] == 0.0 and merged.loc[0, "seedE"] == -6.0
