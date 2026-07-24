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
    q, t, skipped = select_emitted(taxa, genome_paths, mirna_resolved, "Data/Raw/RNAs")
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


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for t in tests:
        t()
        print(f"ok  {t.__name__}")
    print(f"\n{len(tests)} passed")
