from pathlib import Path
import pandas as pd


ALIGNMENT_COLUMNS = ["subseqDP", "hybridDP"]
ENERGY_COLUMNS    = ["E", "E_hybrid", "ED1", "ED2", "Pu1", "Pu2"]
SEED_COLUMNS      = ["seedStart1", "seedEnd1", "seedE", "seedStart2", "seedEnd2"]
SKIP_AS_METADATA  = set(ALIGNMENT_COLUMNS + ENERGY_COLUMNS + SEED_COLUMNS)

LABEL_TARGET = "target 5' "
LABEL_MIRNA  = "miRNA  3' "
LABEL_INDENT = "          "
MARK_3PRIME  = "   3'"
MARK_5PRIME  = "   5'"
SEPARATOR_WIDTH = 39


def _fmt(value):
    return "NA" if pd.isna(value) else str(value)


def _parse_dp(hybrid_dp_field):
    """Return (target_dp, query_dp) from the 'target&query' hybridDP field."""
    if "&" not in str(hybrid_dp_field):
        return str(hybrid_dp_field), ""
    left, right = str(hybrid_dp_field).split("&", 1)
    return left, right


def _parse_subseq(subseq_dp_field):
    """Return (target_seq, query_seq) from the 'target&query' subseqDP field."""
    if "&" not in str(subseq_dp_field):
        return str(subseq_dp_field), ""
    left, right = str(subseq_dp_field).split("&", 1)
    return left, right


def _render_duplex(subseq_dp, hybrid_dp):
    """
    Render an ASCII duplex from IntaRNA's hybridDP dot-bracket notation.

    target reads 5'→3' left-to-right; miRNA reads 3'→5' (displayed reversed).
    Base-paired positions (target_dp[i] == '(') are shown with '|'.
    """
    target_seq, query_seq = _parse_subseq(subseq_dp)
    target_dp,  query_dp  = _parse_dp(hybrid_dp)

    n = max(len(target_seq), len(target_dp))
    indicator = "".join(
        "|" if i < len(target_dp) and target_dp[i] == "(" else " "
        for i in range(n)
    )

    query_display = query_seq[::-1]

    lines = [
        f"{LABEL_TARGET}{target_seq}{MARK_3PRIME}",
        f"{LABEL_INDENT}{indicator}",
        f"{LABEL_MIRNA}{query_display}{MARK_5PRIME}",
    ]
    return "\n".join(lines)


#----- Streams the annotated CSV and emits a per-record human-readable report -----#
def enhance_results(annotated_csv_path, output_path):
    annotated = pd.read_csv(annotated_csv_path)

    if "subseqDP" not in annotated.columns or "hybridDP" not in annotated.columns:
        raise ValueError("annotated CSV must contain 'subseqDP' and 'hybridDP' columns")

    cols = list(annotated.columns)
    metadata_cols = [c for c in cols if c not in SKIP_AS_METADATA]
    energy_cols   = [c for c in ENERGY_COLUMNS if c in cols]
    seed_cols     = [c for c in SEED_COLUMNS   if c in cols]
    subseq_idx    = cols.index("subseqDP")
    hybrid_idx    = cols.index("hybridDP")

    meta_indices   = [cols.index(c) for c in metadata_cols]
    energy_indices = [cols.index(c) for c in energy_cols]
    seed_indices   = [cols.index(c) for c in seed_cols]

    separator = "-" * SEPARATOR_WIDTH

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as fh:
        for row in annotated.itertuples(index=False, name=None):
            for col, idx in zip(metadata_cols, meta_indices):
                fh.write(f"{col}: {_fmt(row[idx])}\n")

            if energy_cols:
                fh.write("\n")
                for col, idx in zip(energy_cols, energy_indices):
                    fh.write(f"{col}: {_fmt(row[idx])}\n")

            if seed_cols:
                fh.write("\n")
                for col, idx in zip(seed_cols, seed_indices):
                    fh.write(f"{col}: {_fmt(row[idx])}\n")

            fh.write("\n")
            fh.write(_render_duplex(row[subseq_idx], row[hybrid_idx]))
            fh.write(f"\n\n{separator}\n\n")


def run_from_snakemake(snakemake):
    enhance_results(
        annotated_csv_path=snakemake.input.annotated,
        output_path=snakemake.output.enhanced,
    )

run_from_snakemake(snakemake)
