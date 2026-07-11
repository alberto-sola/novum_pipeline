from __future__ import annotations

import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from prepare_inputs import select_emitted

def test_emits_only_complete_taxa():
    taxa = ["a", "b", "c"]
    genome_paths = {"a": "Data/Raw/genomes/a.fna", "b": None, "c": "Data/Raw/genomes/c.fna"}
    mirna_resolved = {"a": [("hsa-miR-1", "U")], "b": [("hsa-miR-2", "U")], "c": []}
    q, t, skipped = select_emitted(taxa, genome_paths, mirna_resolved, "Data/Raw/RNAs")
    assert q == {"a": "Data/Raw/RNAs/a.fa"}
    assert t == {"a": "Data/Raw/genomes/a.fna"}
    assert skipped == {"b": "no_assembly", "c": "no_mirnas"}

if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for t in tests:
        t()
        print(f"ok  {t.__name__}")
    print(f"\n{len(tests)} passed")
