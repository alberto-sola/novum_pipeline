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


LABEL_TARGET = "target 5' "
LABEL_MIRNA  = "miRNA  3' "
LABEL_INDENT = "          "
MARK_3PRIME  = "   3'"
MARK_5PRIME  = "   5'"
SEPARATOR    = "-" * 39


#----- Renders one cell as text, mapping pandas NA to the literal "NA" so the report never shows "nan" -----#
def fmt(value):
    return "NA" if pd.isna(value) else str(value)


#----- Raises unless every column in `columns` is present, naming all the missing ones at once -----#
def require_columns(annotated, columns, label):
    missing = [column for column in columns if column not in annotated.columns]
    if missing:
        raise ValueError(f"{label} missing required column(s): {', '.join(missing)}")


#----- Resolves `columns` to (name, itertuples index) pairs, skipping any that are absent. Hoists the lookup out of the per-row loop: `row` is a positional tuple, so a renderer that resolved names per row would pay for it on every record -----#
def column_pairs(annotated, columns):
    index_of = {col: idx for idx, col in enumerate(annotated.columns)}
    return [(column, index_of[column]) for column in columns if column in index_of]


#----- Writes one `name: value` line per pair — the report's metadata-block format, shared by the driver below and the arm-specific energy/seed blocks -----#
def write_pairs(fh, pairs, row):
    for column, idx in pairs:
        fh.write(f"{column}: {fmt(row[idx])}\n")


#----- Per-record report driver: writes the metadata block (every column not in `skip_cols`) for each row, then defers the rest of the record to `render_record(fh, row, index_of)` -----#
def stream_records(annotated, output_path, skip_cols, render_record):
    index_of = {col: idx for idx, col in enumerate(annotated.columns)}
    metadata_pairs = [(col, idx) for col, idx in index_of.items() if col not in skip_cols]

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as fh:
        for row in annotated.itertuples(index=False, name=None):
            write_pairs(fh, metadata_pairs, row)
            render_record(fh, row, index_of)


# Duplex-string columns each renderer consumes. Kept beside the renderers so all
# enhance_* scripts share one source of truth for their skip-lists and guards.
RNAHYBRID_DUPLEX_COLUMNS = ["Target_unmatches", "Target_matches", "miRNA_matches", "miRNA_unmatches"]
INTARNA_DUPLEX_COLUMNS   = ["subseqDP", "hybridDP"]


#----- Writes RNAhybrid's 4-line ASCII duplex (target 5'->3' on top, miRNA 3'->5' below) -----#
def render_rnahybrid_duplex(fh, target_unmatches, target_matches, mirna_matches, mirna_unmatches):
    fh.write(f"{LABEL_TARGET}{fmt(target_unmatches)}{MARK_3PRIME}\n")
    fh.write(f"{LABEL_INDENT}{fmt(target_matches)}\n")
    fh.write(f"{LABEL_INDENT}{fmt(mirna_matches)}\n")
    fh.write(f"{LABEL_MIRNA}{fmt(mirna_unmatches)}{MARK_5PRIME}\n")


#----- Writes IntaRNA's 3-line dot-bracket duplex from subseqDP/hybridDP (miRNA shown reversed) -----#
def render_intarna_duplex(fh, subseq_dp, hybrid_dp):
    # IntaRNA writes both fields as "target&query"; partition yields ("text", "", "") when
    # the separator is absent, so a malformed field degrades to an empty query strand.
    target_seq, _, query_seq = str(subseq_dp).partition("&")
    target_dp, _, _          = str(hybrid_dp).partition("&")
    n = max(len(target_seq), len(target_dp))
    indicator = "".join(
        "|" if i < len(target_dp) and target_dp[i] == "(" else " "
        for i in range(n)
    )
    query_display = query_seq[::-1]
    fh.write(f"{LABEL_TARGET}{target_seq}{MARK_3PRIME}\n")
    fh.write(f"{LABEL_INDENT}{indicator}\n")
    fh.write(f"{LABEL_MIRNA}{query_display}{MARK_5PRIME}\n")
