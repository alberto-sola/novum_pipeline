import pandas as pd

from _report import (
    fmt, stream_records, render_intarna_duplex, INTARNA_DUPLEX_COLUMNS, SEPARATOR,
)


ENERGY_COLUMNS   = ["E", "E_hybrid", "ED1", "ED2", "Pu1", "Pu2"]
SEED_COLUMNS     = ["seedStart1", "seedEnd1", "seedE", "seedStart2", "seedEnd2"]
SKIP_AS_METADATA = set(INTARNA_DUPLEX_COLUMNS + ENERGY_COLUMNS + SEED_COLUMNS)


#----- Streams the annotated CSV and emits a per-record human-readable report -----#
def enhance_results(annotated_csv_path, output_path):
    annotated = pd.read_csv(annotated_csv_path)

    if "subseqDP" not in annotated.columns or "hybridDP" not in annotated.columns:
        raise ValueError("annotated CSV must contain 'subseqDP' and 'hybridDP' columns")

    cols = list(annotated.columns)
    energy_pairs = [(c, cols.index(c)) for c in ENERGY_COLUMNS if c in cols]
    seed_pairs   = [(c, cols.index(c)) for c in SEED_COLUMNS   if c in cols]

    def render_record(fh, row, index_of):
        if energy_pairs:
            fh.write("\n")
            for col, idx in energy_pairs:
                fh.write(f"{col}: {fmt(row[idx])}\n")
        if seed_pairs:
            fh.write("\n")
            for col, idx in seed_pairs:
                fh.write(f"{col}: {fmt(row[idx])}\n")
        fh.write("\n")
        render_intarna_duplex(fh, row[index_of["subseqDP"]], row[index_of["hybridDP"]])
        fh.write(f"\n{SEPARATOR}\n\n")

    stream_records(annotated, output_path, SKIP_AS_METADATA, render_record)


#----- Snakemake entry point: unpacks the injected `snakemake` object and calls the pure logic above -----#
def run_from_snakemake(snakemake):
    enhance_results(
        annotated_csv_path=snakemake.input.annotated,
        output_path=snakemake.output.enhanced,
    )

run_from_snakemake(snakemake)
