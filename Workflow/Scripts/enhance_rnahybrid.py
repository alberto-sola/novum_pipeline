from pathlib import Path
import pandas as pd


ALIGNMENT_COLUMNS = [
    "Target_unmatches",
    "Target_matches",
    "miRNA_matches",
    "miRNA_unmatches"
]

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
    metadata_columns = [c for c in columns if c not in ALIGNMENT_COLUMNS]
    metadata_idx = [columns.index(c) for c in metadata_columns]
    alignment_idx = [columns.index(c) for c in ALIGNMENT_COLUMNS]
    separator = "-" * SEPARATOR_WIDTH

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    # itertuples + per-record streaming write keeps peak memory flat for the
    # potentially huge annotated frame (one row per RNAhybrid hit).
    with output_path.open("w", encoding="utf-8") as fh:
        for row in annotated.itertuples(index=False, name=None):
            for col, idx in zip(metadata_columns, metadata_idx):
                fh.write(f"{col}: {format_value(row[idx])}\n")
            fh.write("\n")
            # The alignment columns form RNAhybrid's ASCII art (target on top,
            # miRNA on bottom); printed without column labels on purpose.
            for idx in alignment_idx:
                fh.write(f"{format_value(row[idx])}\n")
            fh.write(f"\n{separator}\n\n")


def run_from_snakemake(snakemake):
    enhance_results(
        annotated_csv_path=snakemake.input.annotated,
        output_path=snakemake.output.enhanced
    )

run_from_snakemake(snakemake)