import pandas as pd

from _common import ensure_parent
from _length import add_length_corrected, corrected_name


COLUMNS = [
    "Gene", "Gene_length", "miRNA", "miRNA_length",
    "Energy", "P_value", "Position",
    "Target_unmatches", "Target_matches", "miRNA_matches", "miRNA_unmatches",
]

# RNAhybrid's scoring block, in report order. Spelled once here because the enhance reports
# render it as its own block and OUTPUT_COLUMNS orders the CSV from it.
ENERGY_COLUMNS = ["Energy", corrected_name("Energy"), "P_value"]

OUTPUT_COLUMNS = [
    "miRNA", "Gene", *ENERGY_COLUMNS, "Gene_length", "Position",
    "miRNA_unmatches", "miRNA_matches", "Target_matches", "Target_unmatches",
]


#----- Counts rendered target nucleotides preceding the first base pair in RNAhybrid's
#      alignment window. Only the unmatched line can carry them: every column left of the
#      first pair is blank in the matched line by that index's own definition -----#
def _overhang_before_first_pair(unmatched, matched):
    first_pair = next((i for i, c in enumerate(matched) if c != " "), None)
    if first_pair is None:
        return 0
    return sum(1 for u in unmatched[:first_pair] if u != " ")


#----- Parses RNAhybrid's colon-separated table, normalizes Position to a 0-1 fraction of gene length, sorts by Energy -----#
def tidy_rnahybrid(input_path, output_path):
    df = pd.read_csv(
        input_path,
        sep=":",
        header=None,
        names=COLUMNS,
        dtype={"Gene": str, "miRNA": str, "Energy": float, "P_value": float, "Gene_length": int},
        keep_default_na=False,
    )

    if (df["Gene_length"] <= 0).any():
        bad = df.loc[df["Gene_length"] <= 0, "Gene"].tolist()
        raise ValueError(f"Non-positive Gene_length in RNAhybrid output for: {bad}")

    # RNAhybrid reports the first column of its alignment window — an unpaired 5' dangling
    # target nucleotide in ~99.9% of hits — while IntaRNA's start1 is the first PAIRED base.
    # Advance past the overhang so both arms anchor alike; comparing them raw made
    # intersect.py's Site_offset_nt read a systematic +1 nt on genuinely identical sites.
    overhang = [
        _overhang_before_first_pair(u, m)
        for u, m in zip(df["Target_unmatches"], df["Target_matches"])
    ]
    df["Position"]    = (df["Position"].astype(int) + overhang) / df["Gene_length"]

    # Longer genes score better by chance, so Energy alone ranks a pair's candidate genes
    # substantially by length. Added as a column, never as a gate — see _length.py for the
    # measurements and for why RNAhybrid's own -p is the wrong correction.
    df = add_length_corrected(df, "Energy")

    # Ranks only — no gate, no per-pair cap, unlike tidy_intarna: RNAhybrid's -e and -b
    # already did both at the tool. Sort stays on raw Energy (the correction is monotone
    # within a gene, so it never reorders one gene's own sites).
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