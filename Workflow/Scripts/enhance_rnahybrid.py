from pathlib import Path
import pandas as pd


ALIGNMENT_COLUMNS = [
    "Target_unmatches",
    "Target_matches",
    "miRNA_matches",
    "miRNA_unmatches"
]

#----- If a field in ALIGNMENT_COLUMNS is empty (NA), print it as a string -----#
def format_value(value):
    if pd.isna(value):
        return "NA"
    return str(value)


def render_record(row, index):
    #----- Write the column until the ALIGNMENT_COLUMNS are encountered -----#
    metadata_columns = [column for column in row.index if column not in ALIGNMENT_COLUMNS]

    # lines = [f"Alignment {index}"]
    lines = [f"{column}: {format_value(row[column])}" for column in metadata_columns]
    lines.append("")
    # lines.append("Alignment")
    lines.extend(format_value(row[column]) for column in ALIGNMENT_COLUMNS)
    lines.append("")
    lines.append("-" * 39)
    lines.append("")

    return "\n".join(lines)


def enhance_results(annotated_csv_path, output_path):
    annotated = pd.read_csv(annotated_csv_path)

    #----- Write which column of the ALIGNMENT_COLUMNS are not there, i.e. the variable doesn't exist -----#
    missing_columns = [column for column in ALIGNMENT_COLUMNS if column not in annotated.columns]

    #----- If there are missing columns, say it -----#
    if missing_columns:
        missing = ", ".join(missing_columns)
        raise ValueError(f"Missing required alignment column(s): {missing}")

    #----- WHAT'S THIS? -----#
    rendered = [render_record(row, index) for index, (_, row) in enumerate(annotated.iterrows(), start=1)]

    #----- WHAT'S THIS? -----#
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(rendered) + "\n", encoding="utf-8")


def run_from_snakemake(snakemake):
    enhance_results(
        annotated_csv_path=snakemake.input.annotated,
        output_path=snakemake.output.enhanced
    )

run_from_snakemake(snakemake)