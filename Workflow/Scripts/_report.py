"""Shared scaffold for the per-record human-readable reports (enhance_* rules).

Both enhancers write a `name: value` metadata block per hit followed by an
ASCII duplex; only the arm-specific sections (RNAhybrid's 4-line alignment vs
IntaRNA's energy/seed blocks + dot-bracket duplex) differ. This module owns the
shared pieces: the NA-aware formatter, the duplex label constants, and the
streaming driver that pairs columns to indices once and writes via itertuples to
keep peak memory flat on large hit tables. Imports pandas, so it is only pulled
in by the postprocess-env scripts.
"""

from pathlib import Path
import pandas as pd


# Strand labels for the ASCII duplex block. Two spaces after "miRNA" so it
# left-aligns with "target"; the inner indent matches the outer label width so
# the alignment rows stack column-for-column under each other.
LABEL_TARGET = "target 5' "
LABEL_MIRNA  = "miRNA  3' "
LABEL_INDENT = "          "
MARK_3PRIME  = "   3'"
MARK_5PRIME  = "   5'"
SEPARATOR    = "-" * 39


#----- Renders one cell as text, mapping pandas NA to the literal "NA" so the report never shows "nan" -----#
def fmt(value):
    return "NA" if pd.isna(value) else str(value)


#----- Per-record report driver: writes the metadata block (every column not in `skip_cols`) for each row, then defers the rest of the record to `render_record(fh, row, index_of)` -----#
def stream_records(annotated, output_path, skip_cols, render_record):
    cols = list(annotated.columns)
    index_of = {col: idx for idx, col in enumerate(cols)}
    metadata_pairs = [(col, index_of[col]) for col in cols if col not in skip_cols]

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as fh:
        for row in annotated.itertuples(index=False, name=None):
            for col, idx in metadata_pairs:
                fh.write(f"{col}: {fmt(row[idx])}\n")
            render_record(fh, row, index_of)
