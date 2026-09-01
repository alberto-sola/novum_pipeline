from __future__ import annotations

import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pairs import normalize_taxon, taxon_to_query, parse_pairs, is_accession

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

def test_is_accession_recognizes_refseq_and_genbank_only():
    assert is_accession("GCF_000284435.1")
    assert is_accession("GCA_000284435.1")
    assert is_accession("  gcf_000284435.1  ")      # tolerant of case and padding
    # Nothing an organism name or a taxid can produce may match.
    assert not is_accession("escherichia_coli")
    assert not is_accession("candidatus_arthromitus")
    assert not is_accession("511145")
    assert not is_accession("GCF_000284435")        # no version suffix
    assert not is_accession("GCF_00028443.1")       # 8 digits, not 9
    assert not is_accession("GCF_000284435.1.fa")   # trailing junk

def test_accession_survives_the_normalize_query_roundtrip():
    # The lossy half of the round-trip is space<->underscore; an accession must
    # skip it entirely, or GCF_000284435.1 reaches NCBI as "gcf 000284435.1".
    key = normalize_taxon("  gcf_000284435.1 ")
    assert key == "GCF_000284435.1"
    assert taxon_to_query(key) == "GCF_000284435.1"

def test_names_still_lowercase_and_underscore():
    # The accession branch must not disturb the existing name path.
    assert normalize_taxon("Escherichia coli") == "escherichia_coli"
    assert taxon_to_query("escherichia_coli") == "escherichia coli"

def test_parse_pairs_keys_accession_lines_by_accession():
    text = ("GCF_000284435.1\thsa-miR-1226-5p\n"
            "gca_000284435.1:hsa-miR-21-5p\n"
            "Escherichia coli hsa-miR-21-5p\n")
    groups, malformed = parse_pairs(text)
    assert malformed == []
    assert groups["GCF_000284435.1"] == ["hsa-miR-1226-5p"]
    assert groups["GCA_000284435.1"] == ["hsa-miR-21-5p"]
    assert groups["escherichia_coli"] == ["hsa-miR-21-5p"]

def test_accession_and_name_for_one_organism_stay_distinct_keys():
    # Pinning is opt-in per line: a pinned line must not silently merge into
    # the unpinned taxon's bucket (or vice versa).
    groups, _ = parse_pairs("escherichia coli\thsa-miR-21-5p\nGCF_000005845.2\thsa-miR-21-5p\n")
    assert sorted(groups) == ["GCF_000005845.2", "escherichia_coli"]

def test_trailing_comment_is_stripped_not_parsed_as_fields():
    # A pinned line needs a comment saying which strain an opaque GCF_ key is —
    # that is the whole point of the accession form. Without this the last-separator
    # split eats the comment: the miRNA becomes "10953" and the taxon becomes the
    # rest of the line, reported as neither malformed nor pinned.
    text = ("GCF_000153625.3\thsa-mir-515-5p    # F. nucleatum ATCC 10953\n"
            "escherichia coli\thsa-miR-21-5p  # K-12\n")
    groups, malformed = parse_pairs(text)
    assert malformed == []
    assert groups["GCF_000153625.3"] == ["hsa-mir-515-5p"]
    assert groups["escherichia_coli"] == ["hsa-miR-21-5p"]

def test_line_that_is_only_a_comment_after_a_pair_is_still_malformed():
    # Stripping the comment must not turn a commented-out pair into a silent skip
    # of something else: a line with nothing but a comment is a comment, and a line
    # left with one token after stripping is still malformed.
    groups, malformed = parse_pairs("   # just a note\njust_one_token  # trailing\n")
    assert groups == {}
    assert malformed == [(2, "just_one_token  # trailing")]

if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for t in tests:
        t()
        print(f"ok  {t.__name__}")
    print(f"\n{len(tests)} passed")
