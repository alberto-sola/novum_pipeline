"""PyWebView launcher for the Novum Pipeline config editor.

Loads ``UI/Novum Pipeline.html`` in a native window and bridges it to the
local filesystem via a small ``js_api`` exposed as ``window.pywebview.api``
on the JavaScript side.

Setup
-----
- Activate the snakemake conda env first so ``snakemake`` is on PATH::

      source ~/miniconda3/bin/activate snakemake-modern
      pip install pywebview pyqt5 pyqtwebengine

- We force the Qt backend (``webview.start(gui="qt")``) so no system-level
  GTK / WebKit2 packages are required — PyQt5 + PyQtWebEngine pull a
  self-contained Chromium into the conda env.

- Then::

      python UI/launcher.py

Notes
-----
- The launcher resolves the repo root as ``Path(__file__).parent.parent``,
  so it works from any working directory.
- ``run_pipeline`` does not guard against double-clicks; Snakemake's own
  ``.snakemake/locks/`` directory prevents overlapping runs.
- Snakemake's stdout/stderr are inherited by the launcher's terminal —
  that is where to watch the run.
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

# Qt high-DPI: must be set before pywebview imports PyQt5/QtWebEngine.
# QtWebEngine on WSLg ignores --force-device-scale-factor, so page scaling is
# handled in CSS (`html { zoom: ... }` in UI/src/tokens.css) instead.
os.environ.setdefault("QT_AUTO_SCREEN_SCALE_FACTOR", "1")
os.environ.setdefault("QT_ENABLE_HIGHDPI_SCALING", "1")

import webview


ROOT = Path(__file__).resolve().parent.parent
QUERIES_DIR = ROOT / "Data" / "Raw" / "RNAs"
TARGETS_DIR = ROOT / "Data" / "Raw" / "genomes"
CONFIG_PATH = ROOT / "Config" / "config.yaml"
HTML_PATH = ROOT / "UI" / "Novum Pipeline.html"


def _list_dir(directory: Path, prefix: str) -> list[dict]:
    if not directory.is_dir():
        return []
    entries = []
    for entry in sorted(directory.iterdir(), key=lambda p: p.name):
        if entry.is_file():
            entries.append({"name": entry.name, "path": f"{prefix}/{entry.name}"})
    return entries


class API:
    def list_queries(self) -> list[dict]:
        return _list_dir(QUERIES_DIR, "Data/Raw/RNAs")

    def list_targets(self) -> list[dict]:
        return _list_dir(TARGETS_DIR, "Data/Raw/genomes")

    def save_config(self, yaml_text: str) -> dict:
        try:
            CONFIG_PATH.parent.mkdir(parents=True, exist_ok=True)
            CONFIG_PATH.write_text(yaml_text)
            return {"ok": True, "path": str(CONFIG_PATH)}
        except OSError as exc:
            return {"ok": False, "error": str(exc)}

    def run_pipeline(self, yaml_text: str) -> dict:
        save_result = self.save_config(yaml_text)
        if not save_result["ok"]:
            return save_result
        try:
            proc = subprocess.Popen(
                ["snakemake", "--use-conda", "--cores", "all"],
                cwd=str(ROOT),
            )
        except FileNotFoundError as exc:
            return {"ok": False, "error": f"snakemake not found on PATH ({exc})"}
        return {"ok": True, "pid": proc.pid}


def main() -> None:
    webview.create_window(
        "Novum Pipeline",
        str(HTML_PATH),
        js_api=API(),
        width=1280,
        height=820,
        min_size=(960, 640),
        maximized=True,
    )
    # Force the Qt backend on Linux: the GTK backend needs system PyGObject
    # (`python3-gi`), which conda envs don't see. PyQt5 + QtWebEngineWidgets
    # ship via pip into the env and Just Work.
    webview.start(gui="qt")


if __name__ == "__main__":
    main()
