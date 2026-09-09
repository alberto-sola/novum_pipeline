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


def _row(gene, mirna, e_hybrid, ed1=0.0, start1=10, pu1=1.0, pu2=1.0):
    return {
        "id1": gene, "id2": mirna,
        "start1": start1, "end1": start1 + 20, "start2": 1, "end2": 21,
        "E": e_hybrid + ed1, "E_hybrid": e_hybrid, "ED1": ed1, "ED2": 0.0,
        "Pu1": pu1, "Pu2": pu2,
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


def test_target_floor_drops_inaccessible_sites(tmp_path):
    rows = [_row("g1", "m1", -20.0, pu1=0.5),
            _row("g1", "m1", -21.0, start1=50, pu1=1e-6)]
    df = _run(tmp_path, rows, min_target_unpaired_probability=0.001)
    assert list(df["E_hybrid"]) == [-20.0]


def test_target_floor_is_inclusive_at_the_cutoff(tmp_path):
    # Mirrors test_gate_is_inclusive_at_the_cutoff: a floor of 0.001 admits exactly 0.001.
    df = _run(tmp_path, [_row("g1", "m1", -20.0, pu1=0.001)],
              min_target_unpaired_probability=0.001)
    assert len(df) == 1


def test_query_floor_is_independent_of_the_target_floor(tmp_path):
    # Pu1 and Pu2 sit orders of magnitude apart in real data, so the two floors must never
    # be collapsed into one shared value.
    rows = [_row("g1", "m1", -20.0, pu1=1e-9, pu2=0.9),
            _row("g1", "m1", -21.0, start1=50, pu1=0.9, pu2=1e-9)]
    df = _run(tmp_path, rows, min_query_unpaired_probability=0.001)
    assert list(df["E_hybrid"]) == [-20.0]      # kept on Pu2 alone; its Pu1 is ignored


def test_floor_without_its_column_is_rejected(tmp_path):
    # A trimmed intarna.output.columns is legal, but a floor against an absent column
    # would silently admit every row — the opposite of what was configured.
    rows = [_row("g1", "m1", -20.0)]
    del rows[0]["Pu1"]
    with pytest.raises(ValueError, match="Pu1"):
        _run(tmp_path, rows, min_target_unpaired_probability=0.001)


def test_null_floors_keep_every_row(tmp_path):
    df = _run(tmp_path, [_row("g1", "m1", -20.0, pu1=1e-12, pu2=1e-12)])
    assert len(df) == 1


def test_floor_reselects_the_site_before_the_cap(tmp_path):
    # The test that pins the design: strongest site inaccessible, runner-up not. Filtering
    # before the cap keeps the pair, represented by the runner-up; filtering after — once
    # each pair has collapsed to its best-E_hybrid row — would drop it entirely.
    rows = [
        _row("g1", "m1", -25.0, start1=100, pu1=1e-9),   # strongest, inaccessible
        _row("g1", "m1", -18.0, start1=400, pu1=0.20),   # weaker, accessible
    ]
    df = _run(tmp_path, rows, min_target_unpaired_probability=0.001, max_suboptimal_hits=1)
    assert len(df) == 1
    assert df.loc[0, "E_hybrid"] == -18.0
    assert df.loc[0, "Start1"] == 400


# --- batched reading: the file is consumed in chunks, but the result must not know that ---

@pytest.fixture
def tiny_chunks(monkeypatch):
    # Two rows per batch, so every multi-row fixture below spans several — the boundary
    # falls in a different place for each test, which is the point.
    monkeypatch.setattr(tidy_intarna, "READ_CHUNK_ROWS", 2)


def test_cap_is_global_not_per_chunk(tiny_chunks, tmp_path):
    # The failure batching invites: capping inside each batch would keep one row per batch
    # and emit 4, and the survivor would be whichever row led its own batch rather than the
    # best of the file. The cap has to see the whole surviving set.
    rows = [_row("g1", "m1", e, start1=10 * i) for i, e in enumerate([-15.0, -30.0, -20.0,
                                                                     -25.0, -18.0, -22.0,
                                                                     -17.0, -19.0])]
    df = _run(tmp_path, rows, max_suboptimal_hits=1)
    assert len(df) == 1
    assert df.loc[0, "E_hybrid"] == -30.0


def test_the_cap_bounds_what_the_loop_accumulates(tiny_chunks, tmp_path, monkeypatch):
    # The memory failure batching alone does not prevent: the gate is the only filter that
    # runs inside the loop, so with the floors null every gated row is held until the
    # concat. Measured on the shipped panel that is 16.6M rows (36% of the file), not the
    # 285k that reach the output — the cap is what makes the output small, and it used to
    # run only after everything had accumulated. Hoisting it into the loop bounds the
    # accumulation at cap x pairs per batch. The assertion is on rows reaching the concat
    # because that is the peak this rule was OOM-killed at.
    rows = [_row("g1", "m1", e, start1=10 * i) for i, e in enumerate([-15.0, -30.0, -20.0,
                                                                     -25.0, -18.0, -22.0,
                                                                     -17.0, -19.0])]
    accumulated = []
    real_concat = pd.concat

    def spy(frames, *args, **kwargs):
        frames = list(frames)
        accumulated.append(sum(len(frame) for frame in frames))
        return real_concat(frames, *args, **kwargs)

    monkeypatch.setattr(tidy_intarna.pd, "concat", spy)
    df = _run(tmp_path, rows, max_suboptimal_hits=1)

    assert accumulated == [4]                   # one survivor per batch, not all eight
    assert df.loc[0, "E_hybrid"] == -30.0       # and the global winner is still the answer


def test_the_cap_counts_per_pair_not_per_gene(tiny_chunks, tmp_path):
    # One gene, two miRNAs, both inside the same batch. The cap is per PAIR, so both rows
    # are each their pair's only hit and both must survive. A pre-cap grouping by Gene
    # alone would keep just the stronger, and the global cap would then never see m1 at
    # all — a pair silently lost inside the loop, which no amount of global capping undoes.
    rows = [_row("g1", "m1", -20.0, start1=10), _row("g1", "m2", -25.0, start1=20)]
    df = _run(tmp_path, rows, max_suboptimal_hits=1)
    assert sorted(df["miRNA"]) == ["m1", "m2"]


def test_ties_across_batches_keep_the_first_in_file_order(tiny_chunks, tmp_path):
    # Pre-capping per batch reorders each batch by E_hybrid. Tied rows must still resolve
    # the way the un-batched file would resolve them — first occurrence wins — or the cap
    # silently picks a different representative site for the pair.
    # The weaker row leads each batch, so the pre-cap's sort reorders both of them before
    # the tie between start1 20 and 40 is ever compared.
    rows = [_row("g1", "m1", -15.0, start1=10), _row("g1", "m1", -20.0, start1=20),
            _row("g1", "m1", -16.0, start1=30), _row("g1", "m1", -20.0, start1=40)]
    df = _run(tmp_path, rows, max_suboptimal_hits=1)
    assert list(df["Start1"]) == [20]


def test_sort_is_global_not_per_chunk(tiny_chunks, tmp_path):
    rows = [_row("g1", "m1", e, start1=10 * i) for i, e in enumerate([-15.0, -30.0, -20.0,
                                                                     -25.0, -18.0])]
    df = _run(tmp_path, rows)
    assert list(df["E_hybrid"]) == [-30.0, -25.0, -20.0, -18.0, -15.0]


def test_batching_does_not_change_the_result(tiny_chunks, tmp_path):
    # Same fixture as test_floor_reselects_the_site_before_the_cap, read in batches: the
    # floor still reselects across a boundary it now straddles.
    rows = [
        _row("g1", "m1", -25.0, start1=100, pu1=1e-9),
        _row("g1", "m1", -18.0, start1=400, pu1=0.20),
        _row("g2", "m1", -30.0, start1=10, pu1=1e-9),
    ]
    df = _run(tmp_path, rows, min_target_unpaired_probability=0.001, max_suboptimal_hits=1)
    assert list(df["E_hybrid"]) == [-18.0]      # g2 has no accessible site at all


def test_unknown_gene_is_rejected_even_when_the_row_would_be_gated_away(tiny_chunks, tmp_path):
    # Validation runs on the distinct genes BEFORE the filters, so a target FASTA that does
    # not cover the output fails whether or not the offending rows would have survived.
    rows = [_row("g1", "m1", -20.0), _row("ghost", "m1", -1.0, start1=50)]
    with pytest.raises(ValueError, match="ghost"):
        _run(tmp_path, rows, max_hybrid_energy=-12.9)


def test_seed_columns_written_as_NAN_are_read_as_missing(tmp_path):
    # Under `seed: null` IntaRNA emits the literal uppercase NAN, which pandas does not
    # treat as NA by default — left alone the five seed columns are strings, costing a
    # third of the frame and reaching the enhance report as the word "NAN".
    row = _row("g1", "m1", -20.0)
    row.update({c: "NAN" for c in ("seedStart1", "seedEnd1", "seedE",
                                   "seedStart2", "seedEnd2")})
    df = _run(tmp_path, [row])
    assert df["seedE"].isna().all()


def test_a_file_with_no_hits_still_emits_the_schema(tmp_path):
    # pandas yields a single EMPTY batch for a header-only file, which is what lets the
    # concat assume it always has at least one frame. Should an upgrade ever yield nothing
    # instead, this fails here rather than raising "No objects to concatenate" mid-run.
    src = tmp_path / "in.csv"
    pd.DataFrame(columns=list(_row("g1", "m1", -20.0))).to_csv(src, sep=";", index=False)
    out = tmp_path / "out.csv"
    tidy_intarna.tidy_intarna(src, _fasta(tmp_path), out)
    df = pd.read_csv(out)
    assert len(df) == 0
    assert "E_hybrid" in df.columns
