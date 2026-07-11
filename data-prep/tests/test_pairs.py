from __future__ import annotations

import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pairs import normalize_taxon, taxon_to_query, parse_pairs

def test_normalize_and_query_roundtrip():
    assert normalize_taxon("Veillonella_parvula") == "veillonella_parvula"
    assert normalize_taxon("  Fusobacterium nucleatum ") == "fusobacterium_nucleatum"
    assert taxon_to_query("veillonella_parvula") == "veillonella parvula"

def test_separators_tab_colon_spaces():
    text = (
        "Veillonella_parvula\thsa-miR-200b-3p\n"
        "fusobacterium_nucleatum:hsa-miR-200c-3p\n"
        "Dysosmobacter_welbionis    hsa-miR-200c-3p\n"
    )
    groups, malformed = parse_pairs(text)
    assert malformed == []
    assert groups["veillonella_parvula"] == ["hsa-miR-200b-3p"]
    assert groups["fusobacterium_nucleatum"] == ["hsa-miR-200c-3p"]
    assert groups["dysosmobacter_welbionis"] == ["hsa-miR-200c-3p"]

def test_grouping_dedup_order_preserved():
    text = (
        "fusobacterium_nucleatum:hsa-miR-200c-3p\n"
        "fusobacterium_nucleatum\thsa-mir-515-5p\n"
        "Fusobacterium_nucleatum:hsa-miR-200c-3p\n"   # dup (and case-variant taxon)
    )
    groups, _ = parse_pairs(text)
    assert groups["fusobacterium_nucleatum"] == ["hsa-miR-200c-3p", "hsa-mir-515-5p"]

def test_comments_blanks_and_malformed():
    text = "# header\n\nVeillonella_parvula\thsa-miR-200b-3p\njust_one_token\n"
    groups, malformed = parse_pairs(text)
    assert list(groups) == ["veillonella_parvula"]
    assert malformed == [(4, "just_one_token")]

def test_multiword_taxon_space_separated():
    groups, malformed = parse_pairs("Escherichia coli hsa-miR-21\n")
    assert malformed == []
    assert groups == {"escherichia_coli": ["hsa-miR-21"]}

def test_multiword_taxon_colon_and_tab():
    groups, malformed = parse_pairs(
        "Escherichia coli:hsa-miR-21\nFusobacterium nucleatum\thsa-miR-155\n"
    )
    assert malformed == []
    assert groups["escherichia_coli"] == ["hsa-miR-21"]
    assert groups["fusobacterium_nucleatum"] == ["hsa-miR-155"]

if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for t in tests:
        t()
        print(f"ok  {t.__name__}")
    print(f"\n{len(tests)} passed")
