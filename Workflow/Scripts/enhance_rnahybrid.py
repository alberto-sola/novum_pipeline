from pathlib import Path
import pandas as pd


ALIGNMENT_COLUMNS = [
    "Target_unmatches",
    "Target_matches",
    "miRNA_matches",
    "miRNA_unmatches"
]

SEPARATOR_WIDTH = 39


#----- If a field in ALIGNMENT_COLUMNS is empty (NA), print it as a string -----#
def format_value(value):
    if pd.isna(value):
        return "NA"
    return str(value)


def render_record(row):
    metadata_columns = [column for column in row.index if column not in ALIGNMENT_COLUMNS]

    lines = [f"{column}: {format_value(row[column])}" for column in metadata_columns]
    lines.append("")
    # The four alignment columns together form the RNAhybrid ASCII art (target on top,
    # miRNA on bottom); they are printed without column labels on purpose.
    lines.extend(format_value(row[column]) for column in ALIGNMENT_COLUMNS)
    lines.append("")
    lines.append("-" * SEPARATOR_WIDTH)
    lines.append("")

    return "\n".join(lines)


def enhance_results(annotated_csv_path, output_path):
    annotated = pd.read_csv(annotated_csv_path)

    missing_columns = [column for column in ALIGNMENT_COLUMNS if column not in annotated.columns]
    if missing_columns:
        missing = ", ".join(missing_columns)
        raise ValueError(f"Missing required alignment column(s): {missing}")

    rendered = [render_record(row) for _, row in annotated.iterrows()]

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(rendered) + "\n", encoding="utf-8")


def run_from_snakemake(snakemake):
    enhance_results(
        annotated_csv_path=snakemake.input.annotated,
        output_path=snakemake.output.enhanced
    )

run_from_snakemake(snakemake)