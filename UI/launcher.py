"""PyWebView launcher for the Novum Pipeline config editor.

Loads ``UI/Novum Pipeline.html`` in a native window and bridges it to the local
filesystem via an ``API`` instance exposed as ``window.pywebview.api`` on the JS
side (file pickers, config save, and the pipeline run/poll/cancel lifecycle).

We force the Qt backend (``webview.start(gui="qt")``) so no system GTK/WebKit2 is
needed — PyQt5 + PyQtWebEngine pull a self-contained Chromium into the conda env.
See the README for env setup; run with ``python UI/launcher.py``.

``run_pipeline`` records the live ``Popen`` so the UI can poll ``pipeline_status``
and call ``cancel_pipeline``; Snakemake's ``.snakemake/locks/`` still guards against
overlapping runs. Subprocess output is written to ``Data/Results/.pipeline.log``,
which ``pipeline_status`` tails back into the UI on a non-clean exit.
"""
from __future__ import annotations

import os
import signal
import subprocess
import time
from pathlib import Path
from threading import Lock, Timer

# Qt high-DPI: must be set before pywebview imports PyQt5/QtWebEngine.
# QtWebEngine on WSLg ignores --force-device-scale-factor, so page scaling is
# handled in CSS (`html { zoom: ... }` in UI/src/tokens.css) instead.
os.environ.setdefault("QT_AUTO_SCREEN_SCALE_FACTOR", "1")
os.environ.setdefault("QT_ENABLE_HIGHDPI_SCALING", "1")
# WSLg has no GPU passthrough; skip Chromium's GPU init to silence the
# transient `GpuChannelMsg_CreateCommandBuffer` error before fallback.
os.environ.setdefault("QTWEBENGINE_CHROMIUM_FLAGS", "--disable-gpu --log-level=3")

import webview

from config_io import load_config as _load_config, save_config as _save_config


ROOT = Path(__file__).resolve().parent.parent
QUERIES_DIR = ROOT / "Data" / "Raw" / "RNAs"
TARGETS_DIR = ROOT / "Data" / "Raw" / "genomes"
CONFIG_PATH = ROOT / "Config" / "config.yaml"
HTML_PATH = ROOT / "UI" / "Novum Pipeline.html"
LOG_PATH = ROOT / "Data" / "Results" / ".pipeline.log"


#----- Lists files in `directory` as {name, path} dicts the React file pickers consume -----#
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
        self._log_fh = None
        self._reset_run_state()

    #----- Resets the run/latch fields to the idle state (shared by init and acknowledge) -----#
    def _reset_run_state(self) -> None:
        self._proc: subprocess.Popen | None = None
        self._started_at: float | None = None
        self._finished_at: float | None = None
        self._returncode: int | None = None
        self._cancelled = False

    #----- file pickers used by the React form -----#

    def list_queries(self) -> list[dict]:
        return _list_dir(QUERIES_DIR, "Data/Raw/RNAs")

    def list_targets(self) -> list[dict]:
        return _list_dir(TARGETS_DIR, "Data/Raw/genomes")

    #----- config save -----#

    def load_config(self) -> dict:
        return _load_config(CONFIG_PATH)

    def save_config(self, config: dict) -> dict:
        try:
            return _save_config(CONFIG_PATH, config)
        except OSError as exc:
            return {"ok": False, "error": str(exc)}

    #----- pipeline lifecycle -----#

    def _is_running(self) -> bool:
        return self._proc is not None and self._proc.poll() is None

    def _close_log(self) -> None:
        if self._log_fh is not None:
            try:
                self._log_fh.flush()
                self._log_fh.close()
            finally:
                self._log_fh = None

    @staticmethod
    def _read_log_tail(max_lines: int = 80, max_bytes: int = 65536) -> str:
        """Last ``max_lines`` lines of the run log, so a failure can be read
        in-app instead of sending the user back to the launcher terminal.

        Reads only the final ``max_bytes`` so a long run's log is never slurped
        whole; the leading partial line is dropped by the line slice."""
        try:
            with LOG_PATH.open("rb") as fh:
                fh.seek(0, os.SEEK_END)
                fh.seek(max(0, fh.tell() - max_bytes))
                tail = fh.read().decode("utf-8", "replace")
        except OSError:
            return ""
        return "\n".join(tail.splitlines()[-max_lines:])

    def run_pipeline(self, config: dict) -> dict:
        with self._lock:
            if self._is_running():
                return {
                    "ok": False,
                    "error": "Pipeline already running",
                    "pid": self._proc.pid,
                }
            save_result = self.save_config(config)
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
            result = {
                "state": state,
                "pid": self._proc.pid,
                "returncode": self._returncode,
                "elapsed": duration,
                "log_path": str(LOG_PATH),
            }
            # Surface the trace in-app for non-clean exits; success stays terse.
            if state in ("failed", "cancelled"):
                result["log_tail"] = self._read_log_tail()
            return result

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
            self._reset_run_state()
            return {"ok": True}


#----- Builds the PyWebView window, wires the min-width snap-back, and starts the Qt loop -----#
def main() -> None:
    # Width has a floor (cards break below ~960); height is left unconstrained
    # so the user can collapse the window vertically as far as they want.
    min_w = 960
    window = webview.create_window(
        "Novum Pipeline",
        str(HTML_PATH),
        js_api=API(),
        width=1280,
        height=820,
        min_size=(min_w, 1),
    )

    #----- Applies Qt's advisory min-width hint once the native widget exists -----#
    # The hint alone doesn't hold on WSLg (the compositor lets the user drag past it),
    # so the `resized` snap-back below does the real enforcement.
    def _set_min_width_hint(*_args) -> None:
        widget = getattr(window, "native", None)
        if widget is None:
            print("[launcher] min-width: window.native is None; constraint not applied")
            return
        top = widget.window() if hasattr(widget, "window") and callable(widget.window) else widget
        if hasattr(top, "setMinimumWidth"):
            top.setMinimumWidth(min_w)
            print(f"[launcher] min-width {min_w} hint set on {type(top).__name__}")
        else:
            print(f"[launcher] min-width: {type(top).__name__} has no setMinimumWidth")

    #----- Snaps sub-floor widths back to min_w, debounced so the correction fires once -----#
    # Firing on every resize during a drag makes the window drift in the pull direction,
    # since each resize() re-anchors to the just-moved top-left; waiting for the drag to
    # settle gives one clean correction.
    timer: Timer | None = None

    def _snap_back(width, _height, *_args) -> None:
        nonlocal timer
        # Cancel any pending snap: the drag may have crossed back above the floor.
        if timer is not None:
            timer.cancel()
            timer = None
        if width >= min_w:
            return

        def _fire() -> None:
            # Read the live size at fire time so we restore the user's current
            # height, not the (possibly stale) height when the timer was armed.
            try:
                w = window.width
                h = window.height
            except Exception as exc:
                print(f"[launcher] min-width snap-back: live size read failed: {exc}")
                return
            if w >= min_w:
                return
            try:
                window.resize(min_w, h)
            except Exception as exc:
                print(f"[launcher] min-width snap-back failed: {exc}")

        timer = Timer(0.12, _fire)
        timer.daemon = True
        timer.start()

    window.events.shown += _set_min_width_hint
    window.events.resized += _snap_back

    # Force the Qt backend on Linux: the GTK backend needs system PyGObject
    # (`python3-gi`), which conda envs don't see. PyQt5 + QtWebEngineWidgets
    # ship via pip into the env and Just Work.
    # Flip ``debug=True`` for the rare React debugging session — it both opens
    # DevTools at startup and adds the right-click → Inspect Element entry.
    webview.start(gui="qt", debug=False)


if __name__ == "__main__":
    main()
