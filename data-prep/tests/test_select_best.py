from __future__ import annotations

import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fetch_genomes import select_best

def _c(acc, cat, level, date, asm):
    return {"accession": acc, "assembly_info": {
        "refseq_category": cat, "assembly_level": level,
        "release_date": date, "assembly_name": asm}}

def test_reference_beats_better_level_and_date():
    ref = _c("GCF_1.1", "reference genome", "Contig", "2010-01-01", "ASM1")
    rep = _c("GCF_2.1", "representative genome", "Complete Genome", "2020-01-01", "ASM2")
    nun = _c("GCF_3.1", "", "Complete Genome", "2024-01-01", "ASM3")
    assert select_best([rep, nun, ref]) == ("GCF_1.1", "ASM1")

def test_no_category_picks_completeness_then_newest():
    older = _c("GCF_4.1", "", "Complete Genome", "2012-01-01", "ASM4")
    newer = _c("GCF_3.1", "", "Complete Genome", "2024-01-01", "ASM3")
    scaff = _c("GCF_5.1", "", "Scaffold", "2025-01-01", "ASM5")
    assert select_best([older, scaff, newer]) == ("GCF_3.1", "ASM3")

def test_empty_returns_none():
    assert select_best([]) is None


# --- shared auto-discovery footer ---
if __name__ == "__main__":
    import traceback
    passed = 0
    failed = 0
    for name, obj in list(globals().items()):
        if name.startswith("test_") and callable(obj):
            try:
                obj()
                print(f"  PASS  {name}")
                passed += 1
            except Exception:
                print(f"  FAIL  {name}")
                traceback.print_exc()
                failed += 1
    print(f"\n{passed} passed", end="")
    if failed:
        print(f", {failed} failed")
        sys.exit(1)
    else:
        print()
