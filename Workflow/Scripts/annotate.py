import re
import sys
import pandas as pd

from _common import iter_fasta_records, parse_header_id, ensure_parent


# Extracts key/value pairs from FASTA headers like:
# [gene=thrL] [locus_tag=b0001] [protein=thr operon leader peptide]
FIELD_RE = re.compile(r"\[([^=\]]+)=([^\]]*)\]")


#----- Pulls gene/locus/protein metadata out of FASTA headers and returns one deduped row per Gene -----#
def parse_fasta_annotations(fasta_path):
    records = []

    for header, _sequence in iter_fasta_records(fasta_path):
        gene_id = parse_header_id(header)
        fields = dict(FIELD_RE.findall(header))

        records.append({
            "Gene": gene_id,
            "gene_name": fields.get("gene"),
            "locus_tag": fields.get("locus_tag"),
            "protein_name": fields.get("protein") or fields.get("product"),
            "protein_id": fields.get("protein_id"),
        })

    annotations = pd.DataFrame(records)
    duplicates = annotations.duplicated(subset=["Gene"]).sum()
    if duplicates:
        print(
            f"annotate: dropped {duplicates} duplicate Gene row(s) from FASTA annotations.",
            file=sys.stderr,
        )
    annotations = annotations.drop_duplicates(subset=["Gene"])

    return annotations


#----- Left-joins the FASTA annotations onto the tidy table, inserting them after `insert_after` column -----#
def annotate_results(tidy_csv_path, fasta_path, output_path, insert_after):
    tidy = pd.read_csv(tidy_csv_path)

    # Checked before the FASTA is streamed: the anchor column depends on neither, so a
    # misconfigured insert_after should not first cost a pass over a genome-scale file.
    if insert_after not in tidy.columns:
        raise ValueError(f"tidy CSV is missing required column {insert_after!r}")

    annotations = parse_fasta_annotations(fasta_path)
    annotated = tidy.merge(annotations, on="Gene", how="left")

    fasta_columns = [column for column in annotations.columns if column != "Gene"]
    ordered_columns = []
    for column in tidy.columns:
        ordered_columns.append(column)
        if column == insert_after:
            ordered_columns.extend(fasta_columns)

    annotated = annotated[ordered_columns]

    annotated.to_csv(ensure_parent(output_path), index=False)


#----- Snakemake entry point: unpacks the injected `snakemake` object and calls the pure logic above -----#
def run_from_snakemake(snakemake):
    annotate_results(
        tidy_csv_path=snakemake.input.tidy,
        fasta_path=snakemake.input.target,
        output_path=snakemake.output.annotated,
        insert_after=snakemake.params.insert_after,
    )

if "snakemake" in globals():
    run_from_snakemake(snakemake)