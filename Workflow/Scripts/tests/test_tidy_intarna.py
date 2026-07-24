from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # Workflow/Scripts

import pandas as pd
import pytest

import tidy_intarna  # importable only once the trailing call is guarded


def _fasta(tmp_path):
    p = tmp_path / "t.fa"
    p.write_text(">g1 [gene=g1]\n" + "A" * 300 + "\n>g2 [gene=g2]\n" + "A" * 300 + "\n")
    return p


def _intarna_csv(tmp_path, rows):
    """IntaRNA writes ';'-separated CSV with its own lowercase column names."""
    p = tmp_path / "in.csv"
    pd.DataFrame(rows).to_csv(p, sep=";", index=False)
    return p


def _row(gene, mirna, e_hybrid, ed1=0.0, start1=10):
    return {
        "id1": gene, "id2": mirna,
        "start1": start1, "end1": start1 + 20, "start2": 1, "end2": 21,
        "E": e_hybrid + ed1, "E_hybrid": e_hybrid, "ED1": ed1, "ED2": 0.0,
    }


def _run(tmp_path, rows, **kw):
    out = tmp_path / "out.csv"
    tidy_intarna.tidy_intarna(_intarna_csv(tmp_path, rows), _fasta(tmp_path), out, **kw)
    return pd.read_csv(out)


def test_missing_e_hybrid_is_rejected(tmp_path):
    # E_hybrid is the gate quantity now, so a config that trimmed it out of
    # intarna.output.columns must fail loudly rather than silently skip the gate.
    rows = [_row("g1", "m1", -20.0)]
    del rows[0]["E_hybrid"]
    with pytest.raises(ValueError, match="E_hybrid"):
        _run(tmp_path, rows)


def test_gate_drops_rows_above_the_cutoff(tmp_path):
    df = _run(tmp_path, [_row("g1", "m1", -20.0), _row("g1", "m1", -10.0, start1=50)],
              max_hybrid_energy=-12.9)
    assert list(df["E_hybrid"]) == [-20.0]


def test_gate_is_inclusive_at_the_cutoff(tmp_path):
    df = _run(tmp_path, [_row("g1", "m1", -12.9)], max_hybrid_energy=-12.9)
    assert len(df) == 1


def test_null_cutoff_keeps_everything(tmp_path):
    df = _run(tmp_path, [_row("g1", "m1", -20.0), _row("g1", "m1", -2.0, start1=50)],
              max_hybrid_energy=None)
    assert len(df) == 2


def test_output_is_sorted_by_e_hybrid_not_total_e(tmp_path):
    # The yegH case: the strongest duplex is not the best total-E site, because its
    # target region is far less accessible.
    rows = [
        _row("g1", "m1", -22.63, ed1=5.32,  start1=81),    # total E -17.31
        _row("g1", "m1", -25.94, ed1=10.33, start1=1160),  # total E -15.61
    ]
    df = _run(tmp_path, rows)
    assert list(df["E_hybrid"]) == [-25.94, -22.63]


def test_cap_keeps_the_n_best_by_e_hybrid_per_pair(tmp_path):
    rows = [
        _row("g1", "m1", -22.63, ed1=5.32,  start1=81),
        _row("g1", "m1", -25.94, ed1=10.33, start1=1160),
        _row("g2", "m1", -19.0,  start1=10),
    ]
    df = _run(tmp_path, rows, max_suboptimal_hits=1)
    assert len(df) == 2                                  # one row per (miRNA, Gene)
    assert df.loc[df["Gene"] == "g1", "E_hybrid"].iloc[0] == -25.94


def test_gate_and_rank_use_the_same_field(tmp_path):
    # Gating on E_hybrid and ranking by E_hybrid makes order inconsequential: rows
    # passing the gate always sort ahead of failing rows. This case pins that ranking
    # has not reverted to total E; if it had, the row failing the gate (E_hybrid=-10.0,
    # total E=-10.0) would rank first and head(1) would return empty instead of [-20.0].
    rows = [
        _row("g1", "m1", -10.0, ed1=0.0,  start1=10),     # total E -10.0, fails the gate
        _row("g1", "m1", -20.0, ed1=15.0, start1=100),    # total E  -5.0, passes the gate
    ]
    df = _run(tmp_path, rows, max_hybrid_energy=-12.9, max_suboptimal_hits=1)
    assert list(df["E_hybrid"]) == [-20.0]


def test_gate_may_empty_the_output_without_raising(tmp_path):
    df = _run(tmp_path, [_row("g1", "m1", -5.0)], max_hybrid_energy=-12.9)
    assert len(df) == 0
    assert "E_hybrid" in df.columns                       # schema survives an empty gate
