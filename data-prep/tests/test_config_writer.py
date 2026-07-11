from __future__ import annotations

import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from config_writer import replace_block

SAMPLE = """\
#----- miRNAs -----#
queries:
  escherichia: Data/Raw/RNAs/validated_escherichia.fa

#----- Genomes -----#
targets:
  escherichia: Data/Raw/genomes/old.fna

threads: 16
intarna:
  seed:
    query_range: null
""".splitlines()

def test_replace_populated_block_preserves_everything_else():
    out = replace_block(SAMPLE, "queries", {
        "veillonella_parvula": "Data/Raw/RNAs/veillonella_parvula.fa",
        "fusobacterium_nucleatum": "Data/Raw/RNAs/fusobacterium_nucleatum.fa",
    })
    text = "\n".join(out)
    assert "  veillonella_parvula: Data/Raw/RNAs/veillonella_parvula.fa" in text
    assert "  fusobacterium_nucleatum: Data/Raw/RNAs/fusobacterium_nucleatum.fa" in text
    assert "validated_escherichia.fa" not in text          # old entry gone
    assert "#----- miRNAs -----#" in text                  # comment above kept
    assert "#----- Genomes -----#" in text                 # terminator kept
    assert "threads: 16" in text                           # unrelated param kept
    assert "    query_range: null" in text                 # nested key untouched

def test_does_not_match_nested_query_range():
    # 'query_range:' is nested and must NOT be treated as the 'queries' block.
    out = replace_block(SAMPLE, "targets", {"x": "y.fna"})
    text = "\n".join(out)
    assert "  x: y.fna" in text
    assert "    query_range: null" in text
    assert "old.fna" not in text

def test_empty_block_inserts_without_removing():
    lines = ["queries:", "", "#----- next -----#", "threads: 1"]
    out = replace_block(lines, "queries", {"a": "p.fa"})
    assert out == ["queries:", "  a: p.fa", "", "#----- next -----#", "threads: 1"]

def test_missing_key_raises():
    try:
        replace_block(["threads: 1"], "queries", {"a": "b"})
    except ValueError as e:
        assert "queries" in str(e)
    else:
        raise AssertionError("expected ValueError")

if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for t in tests:
        t()
        print(f"ok  {t.__name__}")
    print(f"\n{len(tests)} passed")
