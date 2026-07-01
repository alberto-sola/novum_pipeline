"""pairs.py — tolerant parser for taxon-first (taxon, miRNA) pair input."""
from __future__ import annotations

import re

# Split a line on the first run of tab / colon / spaces. miRNA and taxon tokens
# contain none of these, so the first run is always the field boundary.
_SEP = re.compile(r"[\t: ]+")


def normalize_taxon(raw: str) -> str:
    return raw.strip().lower().replace(" ", "_")


def taxon_to_query(key: str) -> str:
    return key.replace("_", " ")


def parse_pairs(text: str) -> tuple[dict[str, list[str]], list[tuple[int, str]]]:
    groups: dict[str, list[str]] = {}
    malformed: list[tuple[int, str]] = []
    for lineno, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = _SEP.split(line, maxsplit=1)
        if len(parts) != 2 or not parts[0].strip() or not parts[1].strip():
            malformed.append((lineno, raw))
            continue
        taxon = normalize_taxon(parts[0])
        mirna = parts[1].strip()
        bucket = groups.setdefault(taxon, [])
        if mirna not in bucket:
            bucket.append(mirna)
    return groups, malformed
