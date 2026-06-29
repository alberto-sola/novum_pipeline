import pandas as pd

from _report import (
    fmt, stream_records,
    LABEL_TARGET, LABEL_MIRNA, LABEL_INDENT, MARK_3PRIME, MARK_5PRIME, SEPARATOR,
)


ALIGNMENT_COLUMNS = [
    "Target_unmatches",
    "Target_matches",
    "miRNA_matches",
    "miRNA_unmatches",
]


#----- Renders the annotated CSV as a per-record human-readable report (metadata block + RNAhybrid alignment ASCII art) -----#
def enhance_results(annotated_csv_path, output_path):
    annotated = pd.read_csv(annotated_csv_path)

    missing_columns = [column for column in ALIGNMENT_COLUMNS if column not in annotated.columns]
    if missing_columns:
        raise ValueError(f"Missing required alignment column(s): {', '.join(missing_columns)}")

    def render_record(fh, row, index_of):
        fh.write("\n")
        # RNAhybrid's ASCII duplex: target on top reads 5'->3' left-to-right,
        # miRNA on bottom reads 3'->5' left-to-right (antiparallel pairing).
        target_loops  = fmt(row[index_of["Target_unmatches"]])
        target_paired = fmt(row[index_of["Target_matches"]])
        mirna_paired  = fmt(row[index_of["miRNA_matches"]])
        mirna_loops   = fmt(row[index_of["miRNA_unmatches"]])
        fh.write(f"{LABEL_TARGET}{target_loops}{MARK_3PRIME}\n")
        fh.write(f"{LABEL_INDENT}{target_paired}\n")
        fh.write(f"{LABEL_INDENT}{mirna_paired}\n")
        fh.write(f"{LABEL_MIRNA}{mirna_loops}{MARK_5PRIME}\n")
        fh.write(f"\n{SEPARATOR}\n\n")

    stream_records(annotated, output_path, set(ALIGNMENT_COLUMNS), render_record)


#----- Snakemake entry point: unpacks the injected `snakemake` object and calls the pure logic above -----#
def run_from_snakemake(snakemake):
    enhance_results(
        annotated_csv_path=snakemake.input.annotated,
        output_path=snakemake.output.enhanced
    )

run_from_snakemake(snakemake)
