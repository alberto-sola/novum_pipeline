"""pairs.py — tolerant parser for taxon-first (taxon, miRNA) pair input."""
from __future__ import annotations

import re

# Split a line on the last run of tab / colon / spaces. miRNA and taxon tokens
# contain none of these, so the last run is always the field boundary.
_SEP = re.compile(r"[\t: ]+")

_ACCESSION = re.compile(r"GC[AF]_\d{9}\.\d+\Z", re.I)


#----- True for an NCBI assembly accession — the taxon field's "pin this exact assembly" form -----#
def is_accession(s: str) -> bool:
    return bool(_ACCESSION.match(s.strip()))


#----- The taxon key: an accession keeps case and underscore, a name is lowercased and underscored -----#
def normalize_taxon(raw: str) -> str:
    raw = raw.strip()
    # An accession is a literal, not a name: uppercase is its canonical form, and
    # its underscore must survive taxon_to_query below.
    if is_accession(raw):
        return raw.upper()
    return raw.lower().replace(" ", "_")


#----- Inverse of normalize_taxon for the NCBI lookup; identity for an accession -----#
def taxon_to_query(key: str) -> str:
    # Identity for an accession — the underscore->space mapping that de-normalizes
    # a name would split "GCF_000284435.1" into two words before NCBI ever sees it.
    if is_accession(key):
        return key
    return key.replace("_", " ")


#----- Groups "taxon<sep>miRNA" lines by taxon key in first-seen order; returns (groups, malformed) -----#
def parse_pairs(text: str) -> tuple[dict[str, list[str]], list[tuple[int, str]]]:
    groups: dict[str, list[str]] = {}
    malformed: list[tuple[int, str]] = []
    for lineno, raw in enumerate(text.splitlines(), 1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        matches = list(_SEP.finditer(line))
        if not matches:
            malformed.append((lineno, raw))
            continue
        last = matches[-1]
        taxon_raw, mirna = line[: last.start()], line[last.end():].strip()
        if not taxon_raw.strip() or not mirna:
            malformed.append((lineno, raw))
            continue
        taxon = normalize_taxon(taxon_raw)
        bucket = groups.setdefault(taxon, [])
        if mirna not in bucket:
            bucket.append(mirna)
    return groups, malformed
