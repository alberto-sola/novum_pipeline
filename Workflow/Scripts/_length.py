"""Gene-length correction shared by both arms' tidy steps.

Longer genes score better by chance, so ranking a pair's candidate genes by raw energy
ranks them substantially by length: Spearman(Gene_length, Energy) = -0.195 on RNAhybrid
and -0.144 on IntaRNA's E_hybrid, and each pair's top 1% of hits selects genes 1.39x the
pair's own median length. RNAhybrid's own `en = e / ln(m*n)` is NOT the fix — measured,
`-p` selects 0.31x the pool's median length, a worse bias in the other direction.

Measured against "top 1% selected length / pool length" (1.00 = no length preference):

    rank within length decile (ordinal)               0.984
    within-pair OLS residual on ln(L)                 0.952
    Energy - median(Energy | length decile)           0.875
    Energy + theta*ln(L*n)*ln(L), no fit              0.484

The fitted forms lose because the arm's energy gate truncates the population before this
module sees it, understating the slope; the rank-based forms are immune. So: rank within
a length bin, then map that rank onto the middle bin's energy distribution — identical
selection to the ordinal rank (a monotone transform), but the value stays in kcal/mol on
the arm's own scale. Measured at 0.989, Spearman(Gene_length, corrected) = +0.057.

Deliberately an added ranking column, never a replacement gate: it runs on hits that
already passed the arm's gate, so it answers "which of these are good for their length",
not "is this good enough to report".
"""

import numpy as np
import pandas as pd

# Genes per (bacterium x miRNA) pair, binned by length. Measured flatness is 1.003 / 0.989
# / 0.989 / 0.978 at 5 / 10 / 20 / 50 bins — the choice is not delicate anywhere in that
# range. Ten keeps ~230 genes per bin at this project's panel sizes, enough for a stable
# within-bin rank without smearing the length range each bin covers.
LENGTH_BINS = 10

# A pair with fewer hits than this cannot fill the bins densely enough for a within-bin
# rank to mean anything, and is also a pair whose genes span little length range. Such
# pairs keep their uncorrected value rather than a noisy one — the column stays populated
# so downstream sorts never meet a NaN. All 228 pairs in the shipped data clear it.
MIN_HITS = LENGTH_BINS * 5

# One tidy output is one (sample x variant), so the bacterium is already fixed and the
# miRNA is all that is left to group by: each (bacterium x miRNA) pair ranks its own
# candidate genes. Gene_length is the column name both arms' tidy steps emit.
GROUP_COL = "miRNA"
LENGTH_COL = "Gene_length"


#----- Quantile-normalizes one pair's energies to its middle length bin -----#
def _correct_group(energy, length):
    if len(energy) < MIN_HITS:
        return energy

    binned = pd.qcut(length, LENGTH_BINS, labels=False, duplicates="drop")
    reference = energy[binned == binned.max() // 2]
    if len(reference) < 5:
        return energy

    # Rank inside the hit's own length bin, then read that rank off the reference bin.
    # Monotone within any single gene (one length, one bin), so the per-pair cap and
    # intersect.py's per-(miRNA, Gene) best-hit choice are provably unaffected.
    # np.quantile takes the reference unsorted — it computes order statistics itself.
    percentile = pd.Series(energy).groupby(binned).rank(pct=True, method="average")
    return np.quantile(reference, percentile.to_numpy())


#----- The added column's name. Spelled once because both tidy steps' OUTPUT_COLUMNS and
#      intersect's CONSENSUS_COLUMNS select on it, and each of those drops a name it does not
#      recognise silently (filtered out, or reindexed to all-NA) -----#
def corrected_name(energy_col):
    return f"{energy_col}_corrected"


#----- Adds `<energy_col>_corrected`, computed independently within each miRNA -----#
def add_length_corrected(df, energy_col):
    energy = df[energy_col].to_numpy(dtype=float)
    length = df[LENGTH_COL].to_numpy()
    corrected = np.full(len(df), np.nan)

    # groupby(...).indices hands back positions, so slice the two arrays and scatter the
    # result positionally — no label round-trip, and correct even though tidy_intarna
    # calls this on a frame whose index the sort has scrambled. An explicit loop rather
    # than groupby.apply: apply infers its return shape from what the callback hands back,
    # and reshapes a same-length result into a frame on some pandas versions.
    for _, index in df.groupby(GROUP_COL, sort=False).indices.items():
        corrected[index] = _correct_group(energy[index], length[index])

    return df.assign(**{corrected_name(energy_col): corrected})
