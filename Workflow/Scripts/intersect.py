import sys
import pandas as pd

from _common import ensure_parent
from _length import corrected_name
from annotate import ANNOTATION_COLUMNS


# Columns both arms already carry, taken once from the RNAhybrid side and dropped from IntaRNA
# before the merge. The four annotation columns come from annotate.py on both arms (imported,
# so adding one there cannot leave a _x/_y pair here); Gene_length does NOT — IntaRNA's comes
# from the target FASTA, RNAhybrid's from its own output column 2.
SHARED_ANNOTATION = [*ANNOTATION_COLUMNS, "Gene_length"]

# Final consensus column order (minimal-rename: only Position is disambiguated).
CONSENSUS_COLUMNS = [
    "miRNA", "Gene",
    *SHARED_ANNOTATION,
    "Energy", corrected_name("Energy"), "P_value", "Position_rnahybrid",
    "miRNA_unmatches", "miRNA_matches", "Target_matches", "Target_unmatches",
    "E", "E_hybrid", corrected_name("E_hybrid"), "ED1", "ED2", "Pu1", "Pu2",
    "Start1", "End1", "Start2", "End2", "Position_intarna", "Site_offset_nt",
    "subseqDP", "hybridDP",
    "seedStart1", "seedEnd1", "seedE", "seedStart2", "seedEnd2",
]


#----- Keeps one row per `keys` group: the minimum `rank_col` (NaNs dropped first), tie-broken
#      by `tiebreak_col`. The dropna is load-bearing — an all-NaN group must not survive via
#      keep="first" -----#
def select_best_hit(df, keys, rank_col, tiebreak_col):
    ranked = df.dropna(subset=[rank_col]).sort_values([rank_col, tiebreak_col], kind="mergesort")
    return ranked.drop_duplicates(subset=keys, keep="first")


#----- Inner-joins the two arms' representative hits per (miRNA, Gene), each chosen by its own gated energy; shared annotation kept from RNAhybrid, Position disambiguated -----#
def intersect_annotations(rnahybrid_csv, intarna_csv, output_csv):
    keys = ["miRNA", "Gene"]
    rnahybrid = pd.read_csv(rnahybrid_csv)
    intarna = pd.read_csv(intarna_csv)

    rh_best = select_best_hit(rnahybrid, keys, "Energy", "P_value")
    in_best = select_best_hit(intarna, keys, "E_hybrid", "E")

    rh_best = rh_best.rename(columns={"Position": "Position_rnahybrid"})
    in_best = in_best.rename(columns={"Position": "Position_intarna"})

    drop_from_intarna = [c for c in SHARED_ANNOTATION if c in in_best.columns]
    in_best = in_best.drop(columns=drop_from_intarna)

    merged = rh_best.merge(in_best, on=keys, how="inner")

    # The join is pair-granular — coordinates never constrain it, so a consensus row does NOT
    # assert the two arms found the same site. Measure the disagreement rather than gating on
    # it. Signed, not absolute: a one-directional offset is how the anchor mismatch surfaced
    # in the first place, and positive means IntaRNA's site lies downstream.
    if {"Position_rnahybrid", "Position_intarna", "Gene_length"} <= set(merged.columns):
        merged["Site_offset_nt"] = (
            (merged["Position_intarna"] - merged["Position_rnahybrid"])
            * merged["Gene_length"]
        ).round(1)

    absent = [c for c in CONSENSUS_COLUMNS if c not in merged.columns]
    if absent:
        print(
            "intersect: consensus columns absent from inputs (filled NA): "
            + ", ".join(absent)
            + " — trim intarna.output.columns less if you need this metadata.",
            file=sys.stderr,
        )
    merged = merged.reindex(columns=CONSENSUS_COLUMNS)

    merged.to_csv(ensure_parent(output_csv), index=False)
    print(
        f"intersect: {len(merged)} consensus pair(s) "
        f"from {len(rh_best)} RNAhybrid × {len(in_best)} IntaRNA best-hits",
        file=sys.stderr,
    )


#----- Snakemake entry point: unpacks the injected `snakemake` object and calls the pure logic above -----#
def run_from_snakemake(snakemake):
    intersect_annotations(
        rnahybrid_csv=snakemake.input.rnahybrid,
        intarna_csv=snakemake.input.intarna,
        output_csv=snakemake.output.annotated,
    )

if "snakemake" in globals():
    run_from_snakemake(snakemake)
