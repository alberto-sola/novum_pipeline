"""config_writer.py — surgically rewrite queries:/targets: blocks in config.yaml,
preserving every other line (comments, params, formatting)."""
from __future__ import annotations

import os
import re
import tempfile
from pathlib import Path


def replace_block(lines: list[str], key: str, mapping: dict[str, str]) -> list[str]:
    key_re = re.compile(rf"^{re.escape(key)}:\s*$")          # top-level, anchored
    i = next((idx for idx, ln in enumerate(lines) if key_re.match(ln)), None)
    if i is None:
        raise ValueError(f"expected top-level '{key}:' in config")

    # Block body = the run of INDENTED non-blank lines after the key. Blank lines
    # are scanned through but not consumed; a column-0 non-blank line terminates.
    end = i + 1
    for j in range(i + 1, len(lines)):
        ln = lines[j]
        if ln.strip() == "":
            continue
        if ln[:1] in (" ", "\t"):
            end = j + 1
            continue
        break

    new_entries = [f"  {k}: {v}" for k, v in mapping.items()]
    return lines[: i + 1] + new_entries + lines[end:]


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
