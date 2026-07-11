import sys
import pandas as pd

from _common import ensure_parent


# Annotation columns both arms derive identically from the target FASTA;
# taken once from the RNAhybrid side, dropped from IntaRNA before the merge.
SHARED_ANNOTATION = ["gene_name", "locus_tag", "protein_name", "protein_id", "Gene_length"]

# Final consensus column order (minimal-rename: only Position is disambiguated).
CONSENSUS_COLUMNS = [
    "miRNA", "Gene",
    *SHARED_ANNOTATION,
    "Energy", "P_value", "Position_rnahybrid",
    "miRNA_unmatches", "miRNA_matches", "Target_matches", "Target_unmatches",
    "E", "E_hybrid", "ED1", "ED2", "Pu1", "Pu2",
    "Start1", "End1", "Start2", "End2", "Position_intarna",
    "subseqDP", "hybridDP",
    "seedStart1", "seedEnd1", "seedE", "seedStart2", "seedEnd2",
]


#----- Keeps one row per `keys` group: the minimum `rank_col` (NaNs dropped first), tie-broken by `tiebreak_col` -----#
def select_best_hit(df, keys, rank_col, tiebreak_col=None):
    ranked = df.dropna(subset=[rank_col])
    sort_cols = [rank_col] + ([tiebreak_col] if tiebreak_col else [])
    ranked = ranked.sort_values(sort_cols, kind="mergesort")  # stable
    return ranked.drop_duplicates(subset=keys, keep="first")


#----- Inner-joins the two arms' best hits per (miRNA, Gene); shared annotation kept from RNAhybrid, Position disambiguated -----#
def intersect_annotations(rnahybrid_csv, intarna_csv, output_csv):
    keys = ["miRNA", "Gene"]
    rnahybrid = pd.read_csv(rnahybrid_csv)
    intarna = pd.read_csv(intarna_csv)

    rh_best = select_best_hit(rnahybrid, keys, "P_value", "Energy")
    in_best = select_best_hit(intarna, keys, "E")

    rh_best = rh_best.rename(columns={"Position": "Position_rnahybrid"})
    in_best = in_best.rename(columns={"Position": "Position_intarna"})

    drop_from_intarna = [c for c in SHARED_ANNOTATION if c in in_best.columns]
    in_best = in_best.drop(columns=drop_from_intarna)

    merged = rh_best.merge(in_best, on=keys, how="inner")
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
