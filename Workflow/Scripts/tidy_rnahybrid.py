import pandas as pd

from _common import ensure_parent


COLUMNS = [
    "Gene", "Gene_length", "miRNA", "miRNA_length",
    "Energy", "P_value", "Position",
    "Target_unmatches", "Target_matches", "miRNA_matches", "miRNA_unmatches",
]

OUTPUT_COLUMNS = [
    "miRNA", "Gene", "Energy", "P_value", "Gene_length", "Position",
    "miRNA_unmatches", "miRNA_matches", "Target_matches", "Target_unmatches",
]


#----- Parses RNAhybrid's colon-separated table, normalizes Position to a 0-1 fraction of gene length, sorts by Energy -----#
def tidy_rnahybrid(input_path, output_path):
    df = pd.read_csv(
        input_path,
        sep=":",
        header=None,
        names=COLUMNS,
        dtype={"Gene": str, "miRNA": str},
        keep_default_na=False,
    )

    df["Energy"]      = df["Energy"].astype(float)
    df["P_value"]     = df["P_value"].astype(float)
    df["Gene_length"] = df["Gene_length"].astype(int)

    if (df["Gene_length"] <= 0).any():
        bad = df.loc[df["Gene_length"] <= 0, "Gene"].tolist()
        raise ValueError(f"Non-positive Gene_length in RNAhybrid output for: {bad}")

    df["Position"]    = df["Position"].astype(float) / df["Gene_length"]

    # Ranks only — no gate, no per-pair cap, unlike tidy_intarna. RNAhybrid's -e filters
    # exactly the quantity we gate on and -b caps per pair, both at the tool, so repeating
    # either here would be redundant. IntaRNA has no flag that bounds E_hybrid (--outMaxE
    # bounds the total), which is why that arm has to do both downstream.
    df = df.sort_values("Energy", kind="stable")[OUTPUT_COLUMNS]

    df.to_csv(ensure_parent(output_path), index=False)


#----- Snakemake entry point: unpacks the injected `snakemake` object and calls the pure logic above -----#
def run_from_snakemake(snakemake):
    tidy_rnahybrid(
        input_path=snakemake.input.compact,
        output_path=snakemake.output.tidy,
    )

if "snakemake" in globals():
    run_from_snakemake(snakemake)