import pandas as pd

from _report import (
    column_blocks, require_columns, stream_records, write_blocks,
    render_rnahybrid_duplex, render_intarna_duplex,
    RNAHYBRID_DUPLEX_COLUMNS, INTARNA_DUPLEX_COLUMNS,
)
from tidy_intarna import ENERGY_COLUMNS as INTARNA_ENERGY_COLUMNS, SEED_COLUMNS
from tidy_rnahybrid import ENERGY_COLUMNS as RNAHYBRID_ENERGY_COLUMNS


# Duplex-string columns are rendered as ASCII art, not echoed in the metadata block.
DUPLEX_COLUMNS = RNAHYBRID_DUPLEX_COLUMNS + INTARNA_DUPLEX_COLUMNS

# One block per arm's scoring quantities, then IntaRNA's seed block — the same sectioning the
# two per-arm reports use, so no column reads as a section there and as metadata here.
SECTION_COLUMNS  = [RNAHYBRID_ENERGY_COLUMNS, INTARNA_ENERGY_COLUMNS, SEED_COLUMNS]
SKIP_AS_METADATA = set(DUPLEX_COLUMNS).union(*SECTION_COLUMNS)


#----- Streams the consensus CSV: merged metadata block, each arm's scoring block, then the RNAhybrid ASCII duplex and the IntaRNA dot-bracket duplex, per record -----#
def enhance_results(annotated_csv_path, output_path):
    annotated = pd.read_csv(annotated_csv_path)

    require_columns(annotated, DUPLEX_COLUMNS, "consensus CSV")
    blocks = column_blocks(annotated, SECTION_COLUMNS)

    #----- Both arms' scoring blocks, then both duplexes. stream_records adds the separator -----#
    def render_record(fh, row, index_of):
        write_blocks(fh, blocks, row)
        fh.write("\n")
        render_rnahybrid_duplex(fh, row, index_of)
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
