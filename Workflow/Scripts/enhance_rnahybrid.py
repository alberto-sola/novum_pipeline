import pandas as pd

from _report import (
    column_blocks, require_columns, stream_records, write_blocks,
    render_rnahybrid_duplex, RNAHYBRID_DUPLEX_COLUMNS,
)
from tidy_rnahybrid import ENERGY_COLUMNS


SECTION_COLUMNS  = [ENERGY_COLUMNS]
SKIP_AS_METADATA = set(RNAHYBRID_DUPLEX_COLUMNS + ENERGY_COLUMNS)


#----- Renders the annotated CSV as a per-record human-readable report (metadata block + scoring block + RNAhybrid alignment ASCII art) -----#
def enhance_results(annotated_csv_path, output_path):
    annotated = pd.read_csv(annotated_csv_path)

    require_columns(annotated, RNAHYBRID_DUPLEX_COLUMNS, "annotated CSV")
    blocks = column_blocks(annotated, SECTION_COLUMNS)

    #----- The arm's scoring block, then the alignment. stream_records adds the separator -----#
    def render_record(fh, row, index_of):
        write_blocks(fh, blocks, row)
        fh.write("\n")
        render_rnahybrid_duplex(fh, row, index_of)

    stream_records(annotated, output_path, SKIP_AS_METADATA, render_record)


#----- Snakemake entry point: unpacks the injected `snakemake` object and calls the pure logic above -----#
def run_from_snakemake(snakemake):
    enhance_results(
        annotated_csv_path=snakemake.input.annotated,
        output_path=snakemake.output.enhanced
    )

if "snakemake" in globals():
    run_from_snakemake(snakemake)
