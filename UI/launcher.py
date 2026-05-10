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
- ``run_pipeline`` records the live ``Popen`` so the UI can poll
  ``pipeline_status`` and call ``cancel_pipeline``. Snakemake's own
  ``.snakemake/locks/`` directory still prevents overlapping runs.
- Snakemake's stdout/stderr are tee'd to ``Data/Results/.pipeline.log``
  and to the launcher's terminal, so the UI can surface the tail and the
  user can still watch the live output.
"""
from __future__ import annotations

import os
import signal
import subprocess
import time
from pathlib import Path
from threading import Lock

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
LOG_PATH = ROOT / "Data" / "Results" / ".pipeline.log"


def _list_dir(directory: Path, prefix: str) -> list[dict]:
    if not directory.is_dir():
        return []
    entries = []
    for entry in sorted(directory.iterdir(), key=lambda p: p.name):
        if entry.is_file():
            entries.append({"name": entry.name, "path": f"{prefix}/{entry.name}"})
    return entries


class API:
    def __init__(self) -> None:
        self._lock = Lock()
        self._proc: subprocess.Popen | None = None
        self._log_fh = None
        self._started_at: float | None = None
        self._finished_at: float | None = None
        self._returncode: int | None = None
        self._cancelled = False

    # -- file pickers used by the React form --------------------------------

    def list_queries(self) -> list[dict]:
        return _list_dir(QUERIES_DIR, "Data/Raw/RNAs")

    def list_targets(self) -> list[dict]:
        return _list_dir(TARGETS_DIR, "Data/Raw/genomes")

    # -- config save --------------------------------------------------------

    def save_config(self, yaml_text: str) -> dict:
        try:
            CONFIG_PATH.parent.mkdir(parents=True, exist_ok=True)
            CONFIG_PATH.write_text(yaml_text)
            return {"ok": True, "path": str(CONFIG_PATH)}
        except OSError as exc:
            return {"ok": False, "error": str(exc)}

    # -- pipeline lifecycle -------------------------------------------------

    def _is_running(self) -> bool:
        return self._proc is not None and self._proc.poll() is None

    def _close_log(self) -> None:
        if self._log_fh is not None:
            try:
                self._log_fh.flush()
                self._log_fh.close()
            finally:
                self._log_fh = None

    def run_pipeline(self, yaml_text: str) -> dict:
        with self._lock:
            if self._is_running():
                return {
                    "ok": False,
                    "error": "Pipeline already running",
                    "pid": self._proc.pid,
                }
            save_result = self.save_config(yaml_text)
            if not save_result["ok"]:
                return save_result
            try:
                LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
                self._log_fh = LOG_PATH.open("w", buffering=1)
                self._proc = subprocess.Popen(
                    ["snakemake", "--use-conda", "--cores", "all"],
                    cwd=str(ROOT),
                    stdout=self._log_fh,
                    stderr=subprocess.STDOUT,
                )
            except FileNotFoundError as exc:
                self._close_log()
                return {"ok": False, "error": f"snakemake not found on PATH ({exc})"}
            self._started_at = time.time()
            self._finished_at = None
            self._returncode = None
            self._cancelled = False
            return {
                "ok": True,
                "pid": self._proc.pid,
                "log_path": str(LOG_PATH),
            }

    def pipeline_status(self) -> dict:
        with self._lock:
            if self._proc is None:
                return {"state": "idle"}

            if self._proc.poll() is None:
                return {
                    "state": "running",
                    "pid": self._proc.pid,
                    "elapsed": time.time() - (self._started_at or time.time()),
                }

            if self._returncode is None:
                self._returncode = self._proc.returncode
                self._finished_at = time.time()
                self._close_log()

            if self._cancelled:
                state = "cancelled"
            elif self._returncode == 0:
                state = "succeeded"
            else:
                state = "failed"

            duration = (self._finished_at or 0) - (self._started_at or 0)
            return {
                "state": state,
                "pid": self._proc.pid,
                "returncode": self._returncode,
                "elapsed": duration,
                "log_path": str(LOG_PATH),
            }

    def cancel_pipeline(self) -> dict:
        with self._lock:
            if not self._is_running():
                return {"ok": False, "error": "No pipeline running"}
            self._cancelled = True
            try:
                self._proc.send_signal(signal.SIGTERM)
            except OSError as exc:
                return {"ok": False, "error": str(exc)}
            return {"ok": True, "pid": self._proc.pid}

    def acknowledge_pipeline(self) -> dict:
        """Reset the post-run latch so the UI returns to idle.

        The terminal status (succeeded/failed/cancelled) is sticky until the
        UI explicitly acknowledges it — that way a fast finish cannot be
        missed between two polls.
        """
        with self._lock:
            if self._is_running():
                return {"ok": False, "error": "Pipeline still running"}
            self._proc = None
            self._started_at = None
            self._finished_at = None
            self._returncode = None
            self._cancelled = False
            return {"ok": True}

    def read_log_tail(self, n_lines: int = 40) -> dict:
        # Read-whole-file tail: a single snakemake run's log is bounded to
        # kilobytes–low-MB, so seek-from-end isn't worth the complexity.
        try:
            data = LOG_PATH.read_bytes().decode("utf-8", errors="replace")
        except FileNotFoundError:
            return {"ok": True, "lines": []}
        except OSError as exc:
            return {"ok": False, "error": str(exc)}
        return {"ok": True, "lines": data.splitlines()[-n_lines:]}


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
