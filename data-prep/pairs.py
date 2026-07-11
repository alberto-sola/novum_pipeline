"""pairs.py — tolerant parser for taxon-first (taxon, miRNA) pair input."""
from __future__ import annotations

import re

# Split a line on the last run of tab / colon / spaces. miRNA and taxon tokens
# contain none of these, so the last run is always the field boundary.
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
