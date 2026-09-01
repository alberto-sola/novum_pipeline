from __future__ import annotations

import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import tempfile
from pathlib import Path

from prepare_inputs import select_emitted, _write_report

def test_emits_only_complete_taxa():
    taxa = ["a", "b", "c"]
    genome_paths = {"a": "Data/Raw/genomes/a.fna", "b": None, "c": "Data/Raw/genomes/c.fna"}
    mirna_resolved = {"a": [("hsa-miR-1", "U")], "b": [("hsa-miR-2", "U")], "c": []}
    q, t, skipped = select_emitted(taxa, genome_paths, mirna_resolved, "Data/Raw/RNAs", hits={"a": ("A", "1"), "c": ("C", "1")})
    assert q == {"a": "Data/Raw/RNAs/a.fa"}
    assert t == {"a": "Data/Raw/genomes/a.fna"}
    assert skipped == {"b": "no_assembly", "c": "no_mirnas"}


def test_report_has_rna_status_after_genome_path():
    path = os.path.join(tempfile.mkdtemp(), "r.tsv")
    _write_report(path,
                  groups={"ecoli": ["hsa-miR-1"]},
                  genome_paths={"ecoli": "g/GCF_1.1_ASM1_cds_rna_from_genomic.fna"},
                  hits_by_taxon={"ecoli": ("GCF_1.1", "ASM1")},
                  per_mirna={("ecoli", "hsa-miR-1"): ("ok", "hsa-miR-1")},
                  rna_status={"ecoli": "merged"})
    header, row = Path(path).read_text().splitlines()
    assert header.split("\t") == [
        "TAXON", "GENOME_STATUS", "ACCESSION", "ASSEMBLY_NAME", "GENOME_PATH",
        "RNA_STATUS", "MIRNA_INPUT", "MIRNA_STATUS", "MIRNA_RESOLVED"]
    assert row.split("\t")[5] == "merged"


def test_report_rna_status_defaults_to_dash_for_unresolved_taxon():
    path = os.path.join(tempfile.mkdtemp(), "r.tsv")
    _write_report(path,
                  groups={"nope": ["hsa-miR-1"]},
                  genome_paths={},
                  hits_by_taxon={},
                  per_mirna={("nope", "hsa-miR-1"): ("ok", "hsa-miR-1")},
                  rna_status={})
    row = Path(path).read_text().splitlines()[1].split("\t")
    assert row[1] == "no_assembly" and row[5] == "-"


def test_failed_accession_pin_is_labelled_bad_accession():
    # Both fail to resolve, but for different reasons: NCBI has no genome for the
    # organism vs. the accession does not exist. Collapsing them into one status
    # makes a typo look like a biology problem.
    _, _, skipped = select_emitted(
        ["escherichia_coli", "GCF_999999999.9"],
        {"escherichia_coli": None, "GCF_999999999.9": None},
        {"escherichia_coli": [("hsa-miR-1", "U")], "GCF_999999999.9": [("hsa-miR-1", "U")]},
        "Data/Raw/RNAs", hits={})
    assert skipped == {"escherichia_coli": "no_assembly",
                       "GCF_999999999.9": "bad_accession"}


def test_resolved_pin_that_fails_to_download_is_not_called_bad_accession():
    # The accession was fine; the download was not. Calling this bad_accession
    # sends the user hunting for a typo that does not exist. Reachable with an
    # unannotated GCA pin, which resolves but carries no CDS to fetch.
    _, _, skipped = select_emitted(
        ["GCA_000284435.1"], {"GCA_000284435.1": None},
        {"GCA_000284435.1": [("hsa-miR-1", "U")]}, "Data/Raw/RNAs",
        hits={"GCA_000284435.1": ("GCA_000284435.1", "ASM28443v1")})
    assert skipped == {"GCA_000284435.1": "download_failed"}


def test_network_error_on_a_pin_is_not_reported_as_a_typo():
    # bad_accession says "you mistyped this"; a dropped connection says nothing of
    # the sort. Reporting the first for the second sends the user to check an
    # accession that is perfectly correct.
    _, _, skipped = select_emitted(
        ["GCF_000284435.1"], {"GCF_000284435.1": None},
        {"GCF_000284435.1": [("hsa-miR-1", "U")]}, "Data/Raw/RNAs",
        hits={}, lookup_errors={"GCF_000284435.1"})
    assert skipped == {"GCF_000284435.1": "lookup_error"}


def test_report_distinguishes_a_bad_pin_from_an_unknown_organism():
    path = os.path.join(tempfile.mkdtemp(), "r.tsv")
    _write_report(path,
                  groups={"GCF_999999999.9": ["hsa-miR-1"], "nosuchbug": ["hsa-miR-1"]},
                  genome_paths={}, hits_by_taxon={},
                  per_mirna={("GCF_999999999.9", "hsa-miR-1"): ("ok", "hsa-miR-1"),
                             ("nosuchbug", "hsa-miR-1"): ("ok", "hsa-miR-1")},
                  rna_status={})
    rows = [ln.split("\t") for ln in Path(path).read_text().splitlines()[1:]]
    assert rows[0][1] == "bad_accession"
    assert rows[1][1] == "no_assembly"


def test_pinned_accession_becomes_the_sample_key():
    # Characterization: the taxon key is the sample identity, so a pin propagates
    # into the query FASTA name and (downstream) the Data/Results tree.
    q, t, skipped = select_emitted(
        ["GCF_000284435.1"],
        {"GCF_000284435.1": "g/GCF_000284435.1_ASM28443v1_cds_rna_from_genomic.fna"},
        {"GCF_000284435.1": [("hsa-miR-1226-5p", "U")]},
        "Data/Raw/RNAs", hits={"GCF_000284435.1": ("GCF_000284435.1", "ASM28443v1")})
    assert q == {"GCF_000284435.1": "Data/Raw/RNAs/GCF_000284435.1.fa"}
    assert skipped == {}


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for t in tests:
        t()
        print(f"ok  {t.__name__}")
    print(f"\n{len(tests)} passed")
