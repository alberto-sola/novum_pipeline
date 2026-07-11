from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # Workflow/Scripts

import intarna  # importable because the trailing call is guarded


def _outnumber(cmd):
    for tok in cmd:
        if tok.startswith("--outNumber="):
            return int(tok.split("=", 1)[1])
    raise AssertionError("no --outNumber flag emitted")


def test_outnumber_clamped_and_warns(capsys):
    cmd = []
    intarna._add_output_flags(cmd, {"overlap": "B"}, max_suboptimal_hits=5000)
    assert _outnumber(cmd) == intarna.INTARNA_MAX_OUTNUMBER  # 1000
    assert "outNumber" in capsys.readouterr().err            # warned

def test_outnumber_passthrough_under_cap(capsys):
    cmd = []
    intarna._add_output_flags(cmd, {"overlap": "B"}, max_suboptimal_hits=500)
    assert _outnumber(cmd) == 500
    assert capsys.readouterr().err == ""                     # no warning

def test_outnumber_null_defaults_to_cap():
    cmd = []
    intarna._add_output_flags(cmd, {"overlap": "B"}, max_suboptimal_hits=None)
    assert _outnumber(cmd) == intarna.INTARNA_MAX_OUTNUMBER
