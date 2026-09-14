import pandas as pd

from _report import (
    column_blocks, require_columns, stream_records, write_blocks,
    render_intarna_duplex, INTARNA_DUPLEX_COLUMNS,
)
# The producer owns the column groups, so a column added to the tidy schema reaches its own
# block here instead of falling through to the generic metadata lines.
from tidy_intarna import ENERGY_COLUMNS, SEED_COLUMNS


SECTION_COLUMNS  = [ENERGY_COLUMNS, SEED_COLUMNS]
SKIP_AS_METADATA = set(INTARNA_DUPLEX_COLUMNS + ENERGY_COLUMNS + SEED_COLUMNS)


#----- Streams the annotated CSV and emits a per-record human-readable report -----#
def enhance_results(annotated_csv_path, output_path):
    annotated = pd.read_csv(annotated_csv_path)

    require_columns(annotated, INTARNA_DUPLEX_COLUMNS, "annotated CSV")
    blocks = column_blocks(annotated, SECTION_COLUMNS)

    #----- Energy block, seed block, then the duplex. stream_records adds the separator -----#
    def render_record(fh, row, index_of):
        write_blocks(fh, blocks, row)
        fh.write("\n")
        render_intarna_duplex(fh, row, index_of)

    stream_records(annotated, output_path, SKIP_AS_METADATA, render_record)


#----- Snakemake entry point: unpacks the injected `snakemake` object and calls the pure logic above -----#
def run_from_snakemake(snakemake):
    enhance_results(
        annotated_csv_path=snakemake.input.annotated,
        output_path=snakemake.output.enhanced,
    )

if "snakemake" in globals():
    run_from_snakemake(snakemake)
