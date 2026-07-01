import pandas as pd

from _report import (
    stream_records, render_rnahybrid_duplex, render_intarna_duplex,
    RNAHYBRID_DUPLEX_COLUMNS, INTARNA_DUPLEX_COLUMNS, SEPARATOR,
)


# Duplex-string columns are rendered as ASCII art, not echoed in the metadata block.
DUPLEX_COLUMNS = RNAHYBRID_DUPLEX_COLUMNS + INTARNA_DUPLEX_COLUMNS


#----- Streams the consensus CSV: merged metadata block, then RNAhybrid ASCII duplex, then IntaRNA dot-bracket duplex, per record -----#
def enhance_results(annotated_csv_path, output_path):
    annotated = pd.read_csv(annotated_csv_path)

    missing = [c for c in DUPLEX_COLUMNS if c not in annotated.columns]
    if missing:
        raise ValueError(f"consensus CSV missing duplex column(s): {', '.join(missing)}")

    def render_record(fh, row, index_of):
        fh.write("\n")
        render_rnahybrid_duplex(
            fh,
            row[index_of["Target_unmatches"]],
            row[index_of["Target_matches"]],
            row[index_of["miRNA_matches"]],
            row[index_of["miRNA_unmatches"]],
        )
        fh.write("\n")
        render_intarna_duplex(fh, row[index_of["subseqDP"]], row[index_of["hybridDP"]])
        fh.write(f"\n{SEPARATOR}\n\n")

    stream_records(annotated, output_path, set(DUPLEX_COLUMNS), render_record)


#----- Snakemake entry point: unpacks the injected `snakemake` object and calls the pure logic above -----#
def run_from_snakemake(snakemake):
    enhance_results(
        annotated_csv_path=snakemake.input.annotated,
        output_path=snakemake.output.enhanced,
    )

run_from_snakemake(snakemake)
