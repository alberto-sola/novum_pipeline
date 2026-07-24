from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # Workflow/Scripts

import rnahybrid
import _rnahybrid_worker as worker


def test_split_query_sanitizes_filenames_keeps_real_keys(tmp_path):
    q = tmp_path / "q.fa"
    q.write_text(">hsa/miR-21 desc\nACGUACGU\n>hsa-let-7\nUGAGGUAG\n")
    paths = rnahybrid.split_query_per_miRNA(q, tmp_path / "queries")
    assert set(paths) == {"hsa/miR-21", "hsa-let-7"}          # dict keys = real tokens
    for p in paths.values():
        assert p.parent == tmp_path / "queries"               # no nested dirs
        assert "/" not in p.name and p.exists()

def test_build_rnahybrid_command_species_vs_dist():
    base = ["RNAhybrid", "-q", "q.fa", "-t", "c", "-c", "-m", "50000"]
    assert worker.build_rnahybrid_command(
        "RNAhybrid", "q.fa", "c", "50000", species="3utr_human", extra_args=["-b", "1"]
    ) == base + ["-s", "3utr_human", "-b", "1"]
    assert worker.build_rnahybrid_command(
        "RNAhybrid", "q.fa", "c", "50000", distribution="0.1,0.2", extra_args=["-b", "1"]
    ) == base + ["-d", "0.1,0.2", "-b", "1"]

def test_out_path_for_naming():
    assert worker.out_path_for("/out", "q.fa", "chunk_000001", False).name == "output_chunk_000001.tsv"
    assert worker.out_path_for("/out", "/x/hsa_miR-21.fa", "chunk_000001", True).name == "output_hsa_miR-21__chunk_000001.tsv"

def test_energy_flag_comes_from_max_hybrid_energy():
    # -e has always filtered pure hybridization MFE; the old parameter name claimed a
    # total energy the flag never saw.
    assert rnahybrid.build_optional_args(max_hybrid_energy=-18) == ["-e", "-18"]

def test_no_energy_flag_when_cutoff_is_none():
    assert "-e" not in rnahybrid.build_optional_args()
