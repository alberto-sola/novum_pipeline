import pandas as pd

from _report import (
    require_columns, stream_records, render_rnahybrid_duplex,
    RNAHYBRID_DUPLEX_COLUMNS, SEPARATOR,
)


#----- Renders the annotated CSV as a per-record human-readable report (metadata block + RNAhybrid alignment ASCII art) -----#
def enhance_results(annotated_csv_path, output_path):
    annotated = pd.read_csv(annotated_csv_path)

    require_columns(annotated, RNAHYBRID_DUPLEX_COLUMNS, "annotated CSV")

    def render_record(fh, row, index_of):
        fh.write("\n")
        render_rnahybrid_duplex(
            fh,
            row[index_of["Target_unmatches"]],
            row[index_of["Target_matches"]],
            row[index_of["miRNA_matches"]],
            row[index_of["miRNA_unmatches"]],
        )
        fh.write(f"\n{SEPARATOR}\n\n")

    stream_records(annotated, output_path, set(RNAHYBRID_DUPLEX_COLUMNS), render_record)


#----- Snakemake entry point: unpacks the injected `snakemake` object and calls the pure logic above -----#
def run_from_snakemake(snakemake):
    enhance_results(
        annotated_csv_path=snakemake.input.annotated,
        output_path=snakemake.output.enhanced
    )

if "snakemake" in globals():
    run_from_snakemake(snakemake)
