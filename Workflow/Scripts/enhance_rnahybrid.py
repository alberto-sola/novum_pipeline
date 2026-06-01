from pathlib import Path
import pandas as pd


ALIGNMENT_COLUMNS = [
    "Target_unmatches",
    "Target_matches",
    "miRNA_matches",
    "miRNA_unmatches"
]

# Strand labels for the alignment block. Two spaces after "miRNA" so it left-aligns
# vertically with "target". The inner indent matches the outer label width so the
# four lines stack column-for-column under each other.
LABEL_TARGET_OUTER = "target 5' "
LABEL_MIRNA_OUTER  = "miRNA  3' "
LABEL_INNER        = "          "
MARKER_3PRIME      = "   3'"
MARKER_5PRIME      = "   5'"

SEPARATOR_WIDTH = 39


#----- Renders one cell as a string, mapping pandas NA to the literal "NA" so the report never shows "nan" -----#
def format_value(value):
    if pd.isna(value):
        return "NA"
    return str(value)


#----- Renders the annotated CSV as a per-record human-readable report (metadata block + RNAhybrid alignment ASCII art) -----#
def enhance_results(annotated_csv_path, output_path):
    annotated = pd.read_csv(annotated_csv_path)

    missing_columns = [column for column in ALIGNMENT_COLUMNS if column not in annotated.columns]
    if missing_columns:
        raise ValueError(f"Missing required alignment column(s): {', '.join(missing_columns)}")

    columns = list(annotated.columns)
    # Pair each metadata column with its positional index once, ahead of the loop.
    metadata_pairs = [(c, columns.index(c)) for c in columns if c not in ALIGNMENT_COLUMNS]
    target_loops_idx, target_paired_idx, mirna_paired_idx, mirna_loops_idx = (
        columns.index(c) for c in ALIGNMENT_COLUMNS
    )
    separator = "-" * SEPARATOR_WIDTH

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    # itertuples + per-record streaming write keeps peak memory flat for the
    # potentially huge annotated frame (one row per RNAhybrid hit).
    with output_path.open("w", encoding="utf-8") as fh:
        for row in annotated.itertuples(index=False, name=None):
            for col, idx in metadata_pairs:
                fh.write(f"{col}: {format_value(row[idx])}\n")
            fh.write("\n")
            # RNAhybrid's ASCII duplex: target on top reads 5'→3' left-to-right,
            # miRNA on bottom reads 3'→5' left-to-right (antiparallel pairing).
            target_loops  = format_value(row[target_loops_idx])
            target_paired = format_value(row[target_paired_idx])
            mirna_paired  = format_value(row[mirna_paired_idx])
            mirna_loops   = format_value(row[mirna_loops_idx])
            fh.write(f"{LABEL_TARGET_OUTER}{target_loops}{MARKER_3PRIME}\n")
            fh.write(f"{LABEL_INNER}{target_paired}\n")
            fh.write(f"{LABEL_INNER}{mirna_paired}\n")
            fh.write(f"{LABEL_MIRNA_OUTER}{mirna_loops}{MARKER_5PRIME}\n")
            fh.write(f"\n{separator}\n\n")


#----- Snakemake entry point: unpacks the injected `snakemake` object and calls the pure logic above -----#
def run_from_snakemake(snakemake):
    enhance_results(
        annotated_csv_path=snakemake.input.annotated,
        output_path=snakemake.output.enhanced
    )

run_from_snakemake(snakemake)