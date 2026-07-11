from __future__ import annotations

import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import os
from pathlib import Path
from mirnas import index_mature, resolve_mirna, write_query_fasta, ensure_mature_fa, MatureUnavailable

FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures", "mature_small.fa")

def test_index_lowercased_keys_canonical_values():
    idx = index_mature(FIX)
    assert idx["hsa-mir-515-5p"] == ("hsa-miR-515-5p", "UUCUCCAAAAGAAAGCACUUUCUG")
    assert "mmu-mir-292a-3p" in idx

def test_resolve_explicit_prefix_case_insensitive_canonical_out():
    idx = index_mature(FIX)
    assert resolve_mirna("hsa-mir-515-5p", idx) == ("hsa-miR-515-5p", "UUCUCCAAAAGAAAGCACUUUCUG", False)

def test_resolve_bare_defaults_to_hsa():
    idx = index_mature(FIX)
    name, _seq, used = resolve_mirna("miR-200b-3p", idx)
    assert name == "hsa-miR-200b-3p" and used is True

def test_resolve_bare_falls_back_to_mmu():
    idx = index_mature(FIX)
    name, _seq, used = resolve_mirna("mir-292a-3p", idx)   # no hsa- variant exists
    assert name == "mmu-miR-292a-3p" and used is True

def test_resolve_unknown_returns_none_and_no_reprefix():
    idx = index_mature(FIX)
    assert resolve_mirna("hsa-miR-9999-5p", idx) is None   # explicit prefix, not found
    assert resolve_mirna("zzz-miR-1", idx) is None

def test_write_query_fasta():
    import tempfile
    d = tempfile.mkdtemp()
    out = write_query_fasta("demo_taxon", [("hsa-miR-200b-3p", "UAAU"), ("hsa-miR-515-5p", "UUCU")], d)
    assert Path(out).name == "demo_taxon.fa"
    assert Path(out).read_text() == ">hsa-miR-200b-3p\nUAAU\n>hsa-miR-515-5p\nUUCU\n"

def test_ensure_mature_fa_uses_cache_without_touching_network():
    import tempfile
    d = tempfile.mkdtemp()
    cached = Path(d) / "mature.fa"
    cached.write_text(">hsa-miR-1\nUUU\n")
    class BoomSession:                       # any network call is a bug here
        def get(self, *a, **k):
            raise AssertionError("network hit despite a cached mature.fa")
    out = ensure_mature_fa(d, BoomSession())
    assert Path(out) == cached and Path(out).read_text() == ">hsa-miR-1\nUUU\n"

def test_ensure_mature_fa_raises_actionable_error_on_network_failure():
    import tempfile, requests
    d = tempfile.mkdtemp()
    class DeadSession:                       # mimics the mirBase 443 outage
        def get(self, *a, **k):
            raise requests.exceptions.ConnectTimeout("connect timed out")
    try:
        ensure_mature_fa(d, DeadSession())
    except MatureUnavailable as e:
        msg = str(e)
        assert str(Path(d) / "mature.fa") in msg   # tells the user where to drop the file
    else:
        raise AssertionError("expected MatureUnavailable, got no error")

if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for t in tests:
        t()
        print(f"ok  {t.__name__}")
    print(f"\n{len(tests)} passed")
