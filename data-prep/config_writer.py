"""config_writer.py — surgically rewrite queries:/targets: blocks in config.yaml,
preserving every other line (comments, params, formatting)."""
from __future__ import annotations

import os
import re
import tempfile
from pathlib import Path


#----- Swaps one top-level key's entries for `mapping`, carrying the block's comments through -----#
def replace_block(lines: list[str], key: str, mapping: dict[str, str]) -> list[str]:
    key_re = re.compile(rf"^{re.escape(key)}:\s*$")          # top-level, anchored
    i = next((idx for idx, ln in enumerate(lines) if key_re.match(ln)), None)
    if i is None:
        raise ValueError(f"expected top-level '{key}:' in config")

    # The block runs to the next top-level KEY — a column-0 line that is neither
    # blank nor a comment. Only a key may end it: a comment written at column 0
    # inside the block used to terminate it early, so the entries below survived
    # next to the new ones. That is a duplicate YAML mapping key, which loads
    # without complaint (last wins), silently running the pipeline on the old path.
    end = i + 1
    for j in range(i + 1, len(lines)):
        ln = lines[j]
        if ln.strip() == "" or ln[:1] in (" ", "\t") or ln.lstrip().startswith("#"):
            end = j + 1
            continue
        break

    # Entries are replaced wholesale, but blanks and comments are this module's
    # whole reason to exist, so they are kept and re-emitted below the new entries
    # — including any that sit between the last entry and the next key, which is
    # where they were already.
    kept = [ln for ln in lines[i + 1 : end]
            if ln.strip() == "" or ln.lstrip().startswith("#")]
    new_entries = [f"  {k}: {v}" for k, v in mapping.items()]
    return lines[: i + 1] + new_entries + kept + lines[end:]


#----- Rewrites queries: and targets: via a temp file, so a crash never leaves a partial config -----#
def update_path_blocks(config_path, queries: dict[str, str], targets: dict[str, str]) -> None:
    path = Path(config_path)
    lines = path.read_text().splitlines()
    lines = replace_block(lines, "queries", queries)
    lines = replace_block(lines, "targets", targets)
    text = "\n".join(lines) + "\n"

    fd, tmp = tempfile.mkstemp(dir=str(path.parent), suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(text)
        os.replace(tmp, path)               # atomic; never a partial config
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise
