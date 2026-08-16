from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # Workflow/Scripts

import numpy as np
import pandas as pd
import pytest

import _length


def _pair(n_genes, mirna="m1", slope=-2.0, noise=0.0, seed=0):
    """One (bacterium x miRNA) pair whose energy is linear in ln(gene length)."""
    rng = np.random.default_rng(seed)
    length = np.exp(rng.uniform(np.log(150), np.log(4000), n_genes))
    energy = -20.0 + slope * (np.log(length) - np.log(800))
    if noise:
        energy = energy + rng.normal(0, noise, n_genes)
    return pd.DataFrame({
        "miRNA": mirna,
        "Gene": [f"g{i}" for i in range(n_genes)],
        "Gene_length": length.astype(int),
        "Energy": energy,
    })


def _top_length_ratio(df, column, frac=0.01):
    """Median gene length of the best `frac` of hits, over the pool's median."""
    k = max(1, int(round(frac * len(df))))
    selected = df["Gene_length"].to_numpy()[np.argsort(df[column].to_numpy())[:k]]
    return float(np.median(selected)) / float(df["Gene_length"].median())


def test_corrected_column_is_added_and_never_null():
    out = _length.add_length_corrected(_pair(600, noise=1.0), "Energy")
    assert "Energy_corrected" in out.columns
    assert out["Energy_corrected"].notna().all()


def test_pre_existing_columns_are_untouched():
    df = _pair(600, noise=1.0)
    out = _length.add_length_corrected(df, "Energy")
    for column in df.columns:
        assert out[column].equals(df[column])


def test_length_bias_is_removed():
    """The point of the module: raw energy picks long genes, corrected does not."""
    df = _pair(2000, slope=-2.0, noise=1.5, seed=1)
    out = _length.add_length_corrected(df, "Energy")
    assert _top_length_ratio(out, "Energy") > 1.3
    assert _top_length_ratio(out, "Energy_corrected") == pytest.approx(1.0, abs=0.15)


def test_correction_stays_on_the_energy_scale():
    """Quantile-normalised, not rank-transformed — the value is still kcal/mol."""
    out = _length.add_length_corrected(_pair(2000, noise=1.5, seed=2), "Energy")
    assert out["Energy_corrected"].median() == pytest.approx(out["Energy"].median(), abs=1.0)
    assert out["Energy_corrected"].min() >= out["Energy"].min() - 1e-9


def test_monotone_within_one_gene_length():
    """Why the per-pair cap and intersect.py's best-hit choice cannot change: two hits
    of equal gene length keep their order, so a per-gene argmin is preserved."""
    df = _pair(600, noise=1.0, seed=3)
    df.loc[:, "Gene_length"] = 800                      # collapse to one length bin
    out = _length.add_length_corrected(df, "Energy")
    order_raw = np.argsort(out["Energy"].to_numpy(), kind="stable")
    order_new = np.argsort(out["Energy_corrected"].to_numpy(), kind="stable")
    assert (out["Energy"].to_numpy()[order_raw] == out["Energy"].to_numpy()[order_new]).all()


def test_each_group_is_corrected_independently():
    """Two miRNAs in one file must not borrow each other's length distribution."""
    a = _pair(800, mirna="mA", slope=-2.0, noise=1.0, seed=4)
    b = _pair(800, mirna="mB", slope=-6.0, noise=1.0, seed=5)
    out = _length.add_length_corrected(pd.concat([a, b], ignore_index=True), "Energy")
    for mirna in ("mA", "mB"):
        group = out[out.miRNA == mirna]
        assert _top_length_ratio(group, "Energy_corrected") == pytest.approx(1.0, abs=0.2)


def test_small_group_falls_back_to_the_raw_value():
    """Below MIN_HITS a within-bin rank is noise; keep the uncorrected value so the
    column is still populated for downstream sorts."""
    df = _pair(_length.MIN_HITS - 1, noise=1.0, seed=6)
    out = _length.add_length_corrected(df, "Energy")
    assert out["Energy_corrected"].equals(out["Energy"].astype(float))


def test_empty_input_still_yields_the_column():
    df = pd.DataFrame({"miRNA": [], "Gene": [], "Gene_length": [], "Energy": []})
    out = _length.add_length_corrected(df, "Energy")
    assert "Energy_corrected" in out.columns
    assert len(out) == 0


def test_works_for_the_intarna_column_name():
    df = _pair(2000, noise=1.5, seed=7).rename(columns={"Energy": "E_hybrid"})
    out = _length.add_length_corrected(df, "E_hybrid")
    assert "E_hybrid_corrected" in out.columns
    assert _top_length_ratio(out, "E_hybrid_corrected") == pytest.approx(1.0, abs=0.15)
