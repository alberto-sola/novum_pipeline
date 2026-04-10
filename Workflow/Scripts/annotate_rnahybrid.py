from pathlib import Path
import re
import pandas as pd


# This regex extracts key/value pairs from FASTA headers such as:
# [gene=thrL] [locus_tag=b0001] [protein=thr operon leader peptide]
FIELD_RE = re.compile(r"\[([^=\]]+)=([^\]]*)\]")


#----- Takes metadata from the FASTA target genome and returns them as list -----#
def parse_fasta_annotations(fasta_path):
    records = []

    with open(fasta_path) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue

            header = line[1:].strip()
            # Keep the exact leading FASTA identifier as the join key.
            gene_id = header.split(" [", 1)[0]

            # re.findall() returns all [key=value] matches in the header.
            # Converting that list of tuples into a dict makes the fields easy to access.
            fields = dict(FIELD_RE.findall(header))

            records.append(
                {
                    "Gene": gene_id,
                    "gene_name": fields.get("gene"),
                    "locus_tag": fields.get("locus_tag"),
                    "protein_name": fields.get("protein") or fields.get("product"),
                    "protein_id": fields.get("protein_id"),
                }
            )

    annotations = pd.DataFrame(records).drop_duplicates(subset=["Gene"])
    # annotations["display_name"] = (annotations["gene_name"].fillna(annotations["locus_tag"]).fillna(annotations["protein_name"]))

    return annotations


#----- Adds metadata to the table -----#
def annotate_results(tidy_csv_path, fasta_path, output_path):
    tidy = pd.read_csv(tidy_csv_path)
    annotations = parse_fasta_annotations(fasta_path)

    annotated = tidy.merge(annotations, on="Gene", how="left")

    # Insert annotations after the P_value variable
    fasta_columns = [column for column in annotations.columns if column != "Gene"]
    insert_after = "P_value"
    ordered_columns = []

    for column in tidy.columns:
        ordered_columns.append(column)
        if column == insert_after:
            ordered_columns.extend(fasta_columns)

    # Fallback for unexpected input schemas where P_value is absent.
    if insert_after not in tidy.columns:
        ordered_columns = list(tidy.columns) + fasta_columns

    annotated = annotated[ordered_columns]

    #----- If the fallback variable/column is empty, then say it -----#
    # unmatched_rows = annotated["display_name"].isna().sum()
    # if unmatched_rows:
    #     print(
    #         f"Warning: {unmatched_rows} result row(s) could not be matched to FASTA annotations."
    #     )

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    annotated.to_csv(output_path, index=False)


def run_from_snakemake(snakemake):
    annotate_results(
        tidy_csv_path=snakemake.input.tidy,
        fasta_path=snakemake.input.target,
        output_path=snakemake.output.annotated,
    )

run_from_snakemake(snakemake)