"""mirnas.py — miRBase mature.fa: cache/download, index, name resolution, FASTA out."""
from __future__ import annotations

import gzip
from pathlib import Path

# miRBase current mature set. gz is auto-detected by magic bytes, so either the
# plain or the .gz URL works; adjust here if miRBase moves the path.
MIRBASE_URL = "https://www.mirbase.org/download/mature.fa"

# A name is "bare" (no species prefix) when it starts with a miRNA stem.
_STEMS = ("mir-", "let-", "lin-")


def index_mature(path) -> dict[str, tuple[str, str]]:
    index: dict[str, tuple[str, str]] = {}
    name: str | None = None
    seq_parts: list[str] = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                index[name.lower()] = (name, "".join(seq_parts))
            name = line[1:].split()[0]
            seq_parts = []
        elif (s := line.strip()):
            seq_parts.append(s)
    if name is not None:
        index[name.lower()] = (name, "".join(seq_parts))
    return index


def resolve_mirna(raw: str, index: dict[str, tuple[str, str]],
                  prefix_chain: tuple[str, ...] = ("hsa", "mmu")):
    key = raw.strip().lower()
    if key in index:
        canonical, seq = index[key]
        return canonical, seq, False
    if key.startswith(_STEMS):
        for prefix in prefix_chain:
            hit = index.get(f"{prefix}-{key}")
            if hit:
                return hit[0], hit[1], True
    return None


#----- The query-FASTA path for a taxon: <rnas_dir>/<taxon>.fa (single owner) -----#
def query_fasta_path(rnas_dir, taxon_key: str) -> Path:
    return Path(rnas_dir) / f"{taxon_key}.fa"


def write_query_fasta(taxon_key: str, resolved, rnas_dir) -> Path:
    out = query_fasta_path(rnas_dir, taxon_key)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        for canonical, seq in resolved:
            fh.write(f">{canonical}\n{seq}\n")
    return out


def ensure_mature_fa(mirnas_dir, session, force: bool = False, url: str = MIRBASE_URL) -> Path:
    out = Path(mirnas_dir) / "mature.fa"
    if out.exists() and not force:
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    resp = session.get(url, timeout=120)
    resp.raise_for_status()
    data = resp.content
    if url.endswith(".gz") or data[:2] == b"\x1f\x8b":     # gzip magic
        data = gzip.decompress(data)
    out.write_bytes(data)
    return out
