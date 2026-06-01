from pathlib import Path
import pandas as pd


OUTPUT_COLUMNS = [
    "miRNA", "Gene", "E", "E_hybrid", "ED1", "ED2", "Pu1", "Pu2",
    "Gene_length", "Start1", "End1", "Start2", "End2", "Position",
    "subseqDP", "hybridDP",
    "seedStart1", "seedEnd1", "seedE", "seedStart2", "seedEnd2",
]

# Columns that must be present in IntaRNA CSV output for this script to work.
REQUIRED_COLS = {"id1", "id2", "start1", "end1", "start2", "end2", "E"}


#----- Streams the target FASTA once and returns a {Gene → cumulative sequence length} map for downstream position normalization -----#
def _parse_gene_lengths(fasta_path):
    lengths = {}
    current_id = None
    current_len = 0
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    lengths[current_id] = current_len
                header = line[1:].strip()
                current_id = header.split(" [", 1)[0]
                current_len = 0
            else:
                current_len += len(line)
    if current_id is not None:
        lengths[current_id] = current_len
    return lengths


#----- Renames IntaRNA's columns to the pipeline's canonical names, joins gene lengths, computes a 0–1 Position fraction, and emits an energy-sorted CSV -----#
def tidy_intarna(input_path, target_fasta_path, output_path):
    df = pd.read_csv(input_path, sep=";", dtype={"id1": str, "id2": str})

    missing = REQUIRED_COLS - set(df.columns)
    if missing:
        raise ValueError(f"IntaRNA CSV is missing required columns: {sorted(missing)}")

    df = df.rename(columns={"id1": "Gene", "id2": "miRNA",
                             "start1": "Start1", "end1": "End1",
                             "start2": "Start2", "end2": "End2"})

    gene_lengths = _parse_gene_lengths(target_fasta_path)
    df["Gene_length"] = df["Gene"].map(gene_lengths)

    missing_genes = df["Gene_length"].isna().sum()
    if missing_genes:
        bad = df.loc[df["Gene_length"].isna(), "Gene"].unique().tolist()
        raise ValueError(f"Gene_length missing for {missing_genes} record(s): {bad[:5]}")

    df["Gene_length"] = df["Gene_length"].astype(int)

    if (df["Gene_length"] <= 0).any():
        bad = df.loc[df["Gene_length"] <= 0, "Gene"].unique().tolist()
        raise ValueError(f"Non-positive Gene_length in target FASTA for: {bad}")

    # Normalize position to a 0-1 fraction for cross-pipeline comparability with tidy_rnahybrid.
    df["Position"] = df["Start1"].astype(float) / df["Gene_length"]

    df["E"] = df["E"].astype(float)

    present_cols = [c for c in OUTPUT_COLUMNS if c in df.columns]
    df = df.sort_values("E", kind="stable")[present_cols]

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)


#----- Snakemake entry point: unpacks the injected `snakemake` object and calls the pure logic above -----#
def run_from_snakemake(snakemake):
    tidy_intarna(
        input_path=snakemake.input.csv,
        target_fasta_path=snakemake.input.target,
        output_path=snakemake.output.tidy,
    )

run_from_snakemake(snakemake)
