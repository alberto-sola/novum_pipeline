from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # Workflow/Scripts

import pandas as pd

import tidy_rnahybrid  # importable only once the trailing call is guarded


# RNAhybrid compact fields, in order:
#   Gene:Gene_length:miRNA:miRNA_length:Energy:P_value:Position:
#   Target_unmatches:Target_matches:miRNA_matches:miRNA_unmatches
def _line(gene, glen, pos, target_unmatched, target_matched):
    return ":".join([
        gene, str(glen), "m1", "22", "-25.0", "0.1", str(pos),
        target_unmatched, target_matched, " " * len(target_matched), " " * len(target_matched),
    ])


def _run(tmp_path, lines):
    src = tmp_path / "rnahybrid_output.tsv"
    src.write_text("\n".join(lines) + "\n")
    out = tmp_path / "tidy.csv"
    tidy_rnahybrid.tidy_rnahybrid(src, out)
    return pd.read_csv(out)


def test_position_anchors_on_the_first_base_pair_not_the_window_start(tmp_path):
    # Real record (E. coli NP_414542.1_1), verified against the genome: the window starts
    # at 1-based 28 on an unpaired 'A'; the first base pair is at 29, which is exactly
    # where IntaRNA reports start1 for the same site.
    unmatched = "A  A     ACCAUUACCA     AACGG   GG   G    "
    matched   = " CC CCAUC          CAGGU     UGC  GCU     "
    df = _run(tmp_path, [_line("g1", 66, 28, unmatched, matched)])
    assert df["Position"].iloc[0] * 66 == 29


def test_position_is_unchanged_when_the_window_opens_on_a_base_pair(tmp_path):
    # The ~0.1% of hits with no 5' overhang already share IntaRNA's anchor.
    df = _run(tmp_path, [_line("g1", 100, 40, "   ", "GCU")])
    assert df["Position"].iloc[0] * 100 == 40


def test_multi_nucleotide_overhang_is_fully_skipped(tmp_path):
    # Defensive: nothing in the format guarantees a single dangling nucleotide, so the
    # shift must count rendered nucleotides rather than assume +1.
    df = _run(tmp_path, [_line("g1", 100, 40, "AC  ", "  GU")])
    assert df["Position"].iloc[0] * 100 == 42
