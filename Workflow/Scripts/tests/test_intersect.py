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


def _intarna_df(**overrides):
    """Minimal IntaRNA frame — E_hybrid is guaranteed by tidy_intarna's REQUIRED_COLS,
    everything else is trimmable. Override or add columns by keyword; single-value
    defaults broadcast to the length of the longest override, so multi-row cases only
    have to spell out the columns they actually vary."""
    cols = {
        "miRNA": ["m1"], "Gene": ["g1"],
        "E": [-9.5], "E_hybrid": [-9.5],
        "Start1": [5], "End1": [12], "Start2": [1], "End2": [8],
        "Position": [0.2],
    }
    cols.update(overrides)
    rows = max(len(v) for v in cols.values())
    return pd.DataFrame({k: v * rows if len(v) == 1 else v for k, v in cols.items()})


def _write(tmp_path, rh_df, ia_df):
    rh = tmp_path / "rh.csv"; ia = tmp_path / "ia.csv"; out = tmp_path / "out.csv"
    rh_df.to_csv(rh, index=False)
    ia_df.to_csv(ia, index=False)
    intersect.intersect_annotations(rh, ia, out)
    return pd.read_csv(out)


def test_intersect_tolerates_missing_intarna_metadata(tmp_path, capsys):
    merged = _write(tmp_path, _rnahybrid_df(), _intarna_df())   # must NOT raise KeyError

    assert list(merged.columns) == intersect.CONSENSUS_COLUMNS  # full schema preserved
    assert merged["ED1"].isna().all()                          # absent col filled NA
    assert merged.loc[0, "E"] == -9.5
    note = capsys.readouterr().err
    assert "ED1" in note and "Pu1" in note                      # loud-but-graceful note


def test_intersect_full_columns_unchanged(tmp_path):
    ia = _intarna_df(
        ED1=[0.0], ED2=[0.0], Pu1=[1.0], Pu2=[1.0],
        subseqDP=["&"], hybridDP=["&"],
        seedStart1=[5], seedEnd1=[10], seedE=[-6.0],
        seedStart2=[1], seedEnd2=[6],
    )
    merged = _write(tmp_path, _rnahybrid_df(), ia)
    assert merged.loc[0, "ED1"] == 0.0 and merged.loc[0, "seedE"] == -6.0


def test_intarna_representative_is_lowest_E_hybrid(tmp_path):
    # The yegH case. Site 81-104 wins on total E only because its target region is more
    # accessible (ED1 5.32 vs 10.33); site 1160-1187 is the stronger duplex and the one
    # comparable to RNAhybrid's MFE. Ranking by E carries the weaker duplex.
    ia = _intarna_df(
        E=[-17.31, -15.61], E_hybrid=[-22.63, -25.94],
        ED1=[5.32, 10.33], ED2=[0.0, 0.0], Pu1=[1.8e-4, 5.3e-8], Pu2=[1.0, 1.0],
        Start1=[81, 1160], End1=[104, 1187], Start2=[1, 1], End2=[21, 21],
        Position=[0.27, 0.9],
    )
    merged = _write(tmp_path, _rnahybrid_df(), ia)
    assert merged.loc[0, "E_hybrid"] == -25.94      # the stronger duplex, not the lower E
    assert merged.loc[0, "ED1"] == 10.33            # its real accessibility penalty


def test_rnahybrid_representative_is_lowest_energy_with_pvalue_tiebreak(tmp_path):
    rh = pd.DataFrame({
        "miRNA": ["m1", "m1", "m1"], "Gene": ["g1", "g1", "g1"],
        "gene_name": ["G1", "G1", "G1"], "locus_tag": ["L1", "L1", "L1"],
        "protein_name": ["P1", "P1", "P1"], "protein_id": ["PID1", "PID1", "PID1"],
        "Gene_length": [300, 300, 300],
        # Row 2 and row 3 tie on Energy (-30.0, the gated quantity) — row 3 has the
        # better (lower) P_value, so the tiebreak must pick it over row 2.
        "Energy": [-20.0, -30.0, -30.0], "P_value": [0.001, 0.5, 0.2],
        "Position": [0.1, 0.5, 0.6],
        "miRNA_unmatches": ["", "", ""], "miRNA_matches": ["", "", ""],
        "Target_matches": ["", "", ""], "Target_unmatches": ["", "", ""],
    })
    merged = _write(tmp_path, rh, _intarna_df())
    assert merged.loc[0, "Energy"] == -30.0         # gated quantity wins over P_value
    assert merged.loc[0, "P_value"] == 0.2          # tiebreak: better P_value wins the tie


def test_site_offset_is_annotated_not_gated(tmp_path):
    # The join is pair-granular: far-apart sites still join. Measure the distance instead
    # of dropping the row, so "do the tools agree on WHERE?" becomes answerable.
    ia = _intarna_df(Start1=[270], End1=[290], Position=[0.9])
    merged = _write(tmp_path, _rnahybrid_df(), ia)          # RNAhybrid Position = 0.1
    assert len(merged) == 1                                  # not gated away
    assert merged.loc[0, "Site_offset_nt"] == 240.0          # |0.9 - 0.1| * 300
    assert "Site_offset_nt" in intersect.CONSENSUS_COLUMNS
