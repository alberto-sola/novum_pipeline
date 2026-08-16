import pandas as pd

from _report import (
    column_pairs, require_columns, stream_records, write_pairs,
    render_intarna_duplex, INTARNA_DUPLEX_COLUMNS,
)


ENERGY_COLUMNS   = ["E", "E_hybrid", "ED1", "ED2", "Pu1", "Pu2"]
SEED_COLUMNS     = ["seedStart1", "seedEnd1", "seedE", "seedStart2", "seedEnd2"]
SKIP_AS_METADATA = set(INTARNA_DUPLEX_COLUMNS + ENERGY_COLUMNS + SEED_COLUMNS)


#----- Streams the annotated CSV and emits a per-record human-readable report -----#
def enhance_results(annotated_csv_path, output_path):
    annotated = pd.read_csv(annotated_csv_path)

    require_columns(annotated, INTARNA_DUPLEX_COLUMNS, "annotated CSV")

    energy_pairs = column_pairs(annotated, ENERGY_COLUMNS)
    seed_pairs   = column_pairs(annotated, SEED_COLUMNS)

    #----- Energy block, seed block, then the duplex. stream_records adds the separator -----#
    def render_record(fh, row, index_of):
        for pairs in (energy_pairs, seed_pairs):
            if pairs:
                fh.write("\n")
                write_pairs(fh, pairs, row)
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
