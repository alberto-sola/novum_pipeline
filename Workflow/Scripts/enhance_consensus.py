import pandas as pd

from _report import (
    require_columns, stream_records, render_rnahybrid_duplex, render_intarna_duplex,
    RNAHYBRID_DUPLEX_COLUMNS, INTARNA_DUPLEX_COLUMNS,
)


# Duplex-string columns are rendered as ASCII art, not echoed in the metadata block.
DUPLEX_COLUMNS = RNAHYBRID_DUPLEX_COLUMNS + INTARNA_DUPLEX_COLUMNS


#----- Streams the consensus CSV: merged metadata block, then RNAhybrid ASCII duplex, then IntaRNA dot-bracket duplex, per record -----#
def enhance_results(annotated_csv_path, output_path):
    annotated = pd.read_csv(annotated_csv_path)

    require_columns(annotated, DUPLEX_COLUMNS, "consensus CSV")

    #----- Both arms' duplexes, one after the other. stream_records adds the separator -----#
    def render_record(fh, row, index_of):
        fh.write("\n")
        render_rnahybrid_duplex(fh, row, index_of)
        fh.write("\n")
        render_intarna_duplex(fh, row, index_of)

    stream_records(annotated, output_path, set(DUPLEX_COLUMNS), render_record)


#----- Snakemake entry point: unpacks the injected `snakemake` object and calls the pure logic above -----#
def run_from_snakemake(snakemake):
    enhance_results(
        annotated_csv_path=snakemake.input.annotated,
        output_path=snakemake.output.enhanced,
    )

if "snakemake" in globals():
    run_from_snakemake(snakemake)
