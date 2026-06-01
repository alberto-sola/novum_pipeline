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


#----- Renders NaN as "NA", everything else via str() — keeps missing cells from printing "nan" -----#
def _fmt(value):
    return "NA" if pd.isna(value) else str(value)


#----- Splits IntaRNA's "target&query" field on '&'; returns (text, "") when there's no '&' -----#
def _split_target_query(field):
    text = str(field)
    if "&" not in text:
        return text, ""
    left, right = text.split("&", 1)
    return left, right


#----- Builds the 3-line ASCII duplex (target 5'→3', "|" pairing row, miRNA 3'→5' reversed) from subseqDP/hybridDP -----#
def _render_duplex(subseq_dp, hybrid_dp):
    # target reads 5'→3' left-to-right; miRNA reads 3'→5' (displayed reversed).
    # Base-paired positions are target_dp[i] == '('.
    target_seq, query_seq = _split_target_query(subseq_dp)
    target_dp,  query_dp  = _split_target_query(hybrid_dp)

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
    subseq_idx = cols.index("subseqDP")
    hybrid_idx = cols.index("hybridDP")

    # Pair each emitted column with its positional index once, so the per-row
    # loop reads itertuples fields by position without re-resolving names.
    metadata_pairs = [(c, cols.index(c)) for c in cols if c not in SKIP_AS_METADATA]
    energy_pairs   = [(c, cols.index(c)) for c in ENERGY_COLUMNS if c in cols]
    seed_pairs     = [(c, cols.index(c)) for c in SEED_COLUMNS   if c in cols]

    separator = "-" * SEPARATOR_WIDTH

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as fh:
        for row in annotated.itertuples(index=False, name=None):
            for col, idx in metadata_pairs:
                fh.write(f"{col}: {_fmt(row[idx])}\n")

            if energy_pairs:
                fh.write("\n")
                for col, idx in energy_pairs:
                    fh.write(f"{col}: {_fmt(row[idx])}\n")

            if seed_pairs:
                fh.write("\n")
                for col, idx in seed_pairs:
                    fh.write(f"{col}: {_fmt(row[idx])}\n")

            fh.write("\n")
            fh.write(_render_duplex(row[subseq_idx], row[hybrid_idx]))
            fh.write(f"\n\n{separator}\n\n")


#----- Snakemake entry point: unpacks the injected `snakemake` object and calls the pure logic above -----#
def run_from_snakemake(snakemake):
    enhance_results(
        annotated_csv_path=snakemake.input.annotated,
        output_path=snakemake.output.enhanced,
    )

run_from_snakemake(snakemake)
