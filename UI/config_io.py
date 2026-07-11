"""Load/save Config/config.yaml as structured data for the UI bridge.

Pure and dependency-light (stdlib + PyYAML only) so it is unit-testable without
importing pywebview/PyQt. `launcher.py` delegates its API methods here."""
from __future__ import annotations

import os
import tempfile
from pathlib import Path

import yaml


#----- Parse the on-disk config into a plain dict; never raises to the UI -----#
def load_config(path):
    p = Path(path)
    try:
        text = p.read_text()
    except OSError:
        return {}
    try:
        data = yaml.safe_load(text)
    except yaml.YAMLError as exc:
        return {"__error__": str(exc)}
    return data if isinstance(data, dict) else {}


#----- Emit `config` as YAML with safe quoting, written atomically (temp + os.replace) -----#
def save_config(path, config):
    p = Path(path)
    text = yaml.safe_dump(config, sort_keys=False, default_flow_style=False, allow_unicode=True)
    p.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=str(p.parent), suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(text)
        os.replace(tmp, p)                  # atomic; never a partial config
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise
    return {"ok": True, "path": str(p)}
