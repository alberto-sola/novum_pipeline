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


import json

import pytest

LADDER = [76, 150, 300, 600, 900, 1500, 3000]


def _record(name, length):
    return f">{name}\n{'A' * length}\n"


# --- chunking: the broadcast path must not move ---

def test_chunking_without_anchors_keeps_todays_paths_and_contents(tmp_path):
    fasta = tmp_path / "t.fna"
    fasta.write_text(">a\nAAAA\n>b\nCCCC\n>c\nGGGG\n")
    chunks = rnahybrid.write_fasta_chunks(fasta, tmp_path / "chunks", max_lines=4)
    assert [p.name for p, _ in chunks] == ["chunk_000000", "chunk_000001"]
    assert [anchor for _, anchor in chunks] == [None, None]
    assert (tmp_path / "chunks" / "chunk_000000").read_text() == ">a\nAAAA\n>b\nCCCC\n"
    assert (tmp_path / "chunks" / "chunk_000001").read_text() == ">c\nGGGG\n"

def test_chunking_still_rejects_an_empty_target(tmp_path):
    fasta = tmp_path / "t.fna"
    fasta.write_text("")
    with pytest.raises(RuntimeError, match="No FASTA records"):
        rnahybrid.write_fasta_chunks(fasta, tmp_path / "chunks")


# --- chunking: the stratified path ---

def test_chunking_groups_records_by_cell_in_ascending_anchor_order(tmp_path):
    fasta = tmp_path / "t.fna"
    # 5000 -> open top cell (3000); 50 and 45 -> open bottom cell (76); 736 -> 900,
    # since the geometric midpoint of 600/900 is 734.8.
    fasta.write_text(_record("big", 5000) + _record("tiny", 50)
                     + _record("mid", 736) + _record("tiny2", 45))
    chunks = rnahybrid.write_fasta_chunks(fasta, tmp_path / "chunks", anchors=LADDER)
    assert [anchor for _, anchor in chunks] == [76, 900, 3000]
    bottom = chunks[0][0].read_text()
    assert bottom.startswith(">tiny\n") and ">tiny2\n" in bottom   # file order within a cell
    assert chunks[1][0].read_text().startswith(">mid\n")
    assert chunks[2][0].read_text().startswith(">big\n")

def test_chunking_skips_empty_cells_entirely(tmp_path):
    fasta = tmp_path / "t.fna"
    fasta.write_text(_record("only", 736))
    chunks = rnahybrid.write_fasta_chunks(fasta, tmp_path / "chunks", anchors=LADDER)
    assert [anchor for _, anchor in chunks] == [900]

def test_chunking_bounds_lines_within_a_cell_and_numbers_globally(tmp_path):
    fasta = tmp_path / "t.fna"
    # Four 2-line records in the 76 cell, two in the 900 cell, max_lines=4: the bottom
    # cell splits and the numbering carries across the boundary.
    fasta.write_text("".join(_record(f"s{i}", 50) for i in range(4))
                     + "".join(_record(f"m{i}", 736) for i in range(2)))
    chunks = rnahybrid.write_fasta_chunks(fasta, tmp_path / "chunks",
                                          anchors=LADDER, max_lines=4)
    assert [(p.name, a) for p, a in chunks] == [
        ("chunk_000000", 76), ("chunk_000001", 76), ("chunk_000002", 900),
    ]

def test_chunking_with_anchors_still_rejects_an_empty_target(tmp_path):
    fasta = tmp_path / "t.fna"
    fasta.write_text("")
    with pytest.raises(RuntimeError, match="No FASTA records"):
        rnahybrid.write_fasta_chunks(fasta, tmp_path / "chunks", anchors=LADDER)

def test_chunking_with_anchors_never_splits_a_record_over_max_lines(tmp_path):
    fasta = tmp_path / "t.fna"
    fasta.write_text(_record("solo", 50))   # 2 lines > max_lines=1, must still get its own chunk
    chunks = rnahybrid.write_fasta_chunks(fasta, tmp_path / "chunks",
                                          anchors=LADDER, max_lines=1)
    assert [anchor for _, anchor in chunks] == [76]
    assert chunks[0][0].read_text() == _record("solo", 50)

def test_sorted_glob_order_equals_emission_order(tmp_path):
    # merge_output_files needs sorted(glob) to order numerically, which holds only while
    # the counter is global and zero-padded.
    fasta = tmp_path / "t.fna"
    fasta.write_text("".join(_record(f"s{i}", 50) for i in range(12))
                     + "".join(_record(f"m{i}", 736) for i in range(12)))
    chunks = rnahybrid.write_fasta_chunks(fasta, tmp_path / "chunks",
                                          anchors=LADDER, max_lines=4)
    assert [p.name for p, _ in chunks] == sorted(p.name for p, _ in chunks)


# --- loading the distributions ---

def _write_calibration(tmp_path, payload):
    path = tmp_path / "rnacalibrate.json"
    path.write_text(json.dumps(payload))
    return path

FLAT = {"calibration": {"per_query": [
    {"query": "miR-1", "xi": 2.5, "theta": 0.25},
]}}

LADDERED = {"calibration": {
    "anchors": [{"anchor": 76}, {"anchor": 900}],
    "per_query": [{"query": "miR-1", "xi": 2.5, "theta": 0.25, "strata": [
        {"anchor": 76, "xi": 2.55, "theta": 0.2532},
        {"anchor": 900, "xi": 2.48, "theta": 0.1673},
    ]}],
}}

def test_load_distributions_returns_nothing_on_the_uncalibrated_path():
    assert rnahybrid.load_per_query_distributions(None) == (None, None)

def test_load_distributions_keys_a_flat_json_on_a_none_anchor(tmp_path):
    anchors, dist = rnahybrid.load_per_query_distributions(_write_calibration(tmp_path, FLAT))
    assert anchors is None
    assert dist == {("miR-1", None): "2.500000,0.250000"}

def test_load_distributions_keys_a_ladder_json_on_query_and_anchor(tmp_path):
    anchors, dist = rnahybrid.load_per_query_distributions(
        _write_calibration(tmp_path, LADDERED))
    assert anchors == (76, 900)
    assert dist == {("miR-1", 76): "2.550000,0.253200",
                    ("miR-1", 900): "2.480000,0.167300"}

def test_load_distributions_raises_when_a_ladder_json_has_no_strata(tmp_path):
    broken = {"calibration": {"anchors": [{"anchor": 76}],
                              "per_query": [{"query": "miR-1", "xi": 2.5, "theta": 0.25}]}}
    with pytest.raises(RuntimeError, match="strata"):
        rnahybrid.load_per_query_distributions(_write_calibration(tmp_path, broken))

def test_load_distributions_raises_on_an_empty_calibration(tmp_path):
    with pytest.raises(RuntimeError, match="No per-query distributions"):
        rnahybrid.load_per_query_distributions(
            _write_calibration(tmp_path, {"calibration": {"per_query": []}}))


# --- the job spec ---

def test_job_spec_pairs_each_chunk_with_its_own_cells_xi_theta(tmp_path):
    spec = rnahybrid.build_job_spec_tsv(
        query_paths_by_name={"miR-1": tmp_path / "miR-1.fa"},
        dist_map={("miR-1", 76): "A,B", ("miR-1", 900): "C,D"},
        chunks=[(tmp_path / "chunk_000000", 76), (tmp_path / "chunk_000001", 900)],
        tsv_path=tmp_path / "spec.tsv",
    )
    rows = [line.split("\t") for line in spec.read_text().splitlines()]
    assert [row[1] for row in rows] == ["A,B", "C,D"]
    assert [Path(row[2]).name for row in rows] == ["chunk_000000", "chunk_000001"]

def test_job_spec_raises_naming_the_missing_query_and_anchor(tmp_path):
    with pytest.raises(RuntimeError, match=r"miR-1 @ 900"):
        rnahybrid.build_job_spec_tsv(
            query_paths_by_name={"miR-1": tmp_path / "miR-1.fa"},
            dist_map={("miR-1", 76): "A,B"},
            chunks=[(tmp_path / "chunk_000000", 76), (tmp_path / "chunk_000001", 900)],
            tsv_path=tmp_path / "spec.tsv",
        )
