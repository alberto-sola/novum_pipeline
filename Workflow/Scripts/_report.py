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


#----- Splits IntaRNA's "target&query" field on '&'; returns (text, "") when there's no '&' -----#
def _split_target_query(field):
    text = str(field)
    if "&" not in text:
        return text, ""
    left, right = text.split("&", 1)
    return left, right


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
    target_seq, query_seq = _split_target_query(subseq_dp)
    target_dp, _query_dp  = _split_target_query(hybrid_dp)
    n = max(len(target_seq), len(target_dp))
    indicator = "".join(
        "|" if i < len(target_dp) and target_dp[i] == "(" else " "
        for i in range(n)
    )
    query_display = query_seq[::-1]
    fh.write(f"{LABEL_TARGET}{target_seq}{MARK_3PRIME}\n")
    fh.write(f"{LABEL_INDENT}{indicator}\n")
    fh.write(f"{LABEL_MIRNA}{query_display}{MARK_5PRIME}\n")
