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


def test_helix_flags_emitted_for_each_key():
    cmd = []
    intarna._add_helix_flags(cmd, {
        "min_bp": 3,
        "max_bp": 8,
        "max_internal_loop": 2,
        "min_unpaired_probability": 0.05,
        "max_energy": -1.5,
        "full_energy": True,
    })
    assert cmd == [
        "--helixMinBP=3",
        "--helixMaxBP=8",
        "--helixMaxIL=2",
        "--helixMinPu=0.05",
        "--helixMaxE=-1.5",
        "--helixFullE",
    ]

def test_helix_flags_skip_nulls_and_false():
    cmd = []
    intarna._add_helix_flags(cmd, {
        "min_bp": None,
        "max_bp": 4,
        "max_internal_loop": None,
        "min_unpaired_probability": None,
        "max_energy": None,
        "full_energy": False,
    })
    assert cmd == ["--helixMaxBP=4"]

def test_helix_flags_empty_block_emits_nothing():
    cmd = []
    intarna._add_helix_flags(cmd, {})
    assert cmd == []
    intarna._add_helix_flags(cmd, None)
    assert cmd == []

def test_helix_zero_is_a_real_value_not_a_skip():
    # 0 is meaningful for max_internal_loop (pure stacks) and max_energy.
    # `is not None` must not collapse it the way a falsy check would.
    cmd = []
    intarna._add_helix_flags(cmd, {"max_internal_loop": 0, "max_energy": 0})
    assert cmd == ["--helixMaxIL=0", "--helixMaxE=0"]


def _flag(cmd, prefix):
    for tok in cmd:
        if tok.startswith(prefix):
            return tok
    return None


def _build(cfg, **kw):
    kw.setdefault("acc", "N")
    kw.setdefault("query", "q.fa")
    kw.setdefault("target", "t.fa")
    kw.setdefault("out_path", "out.csv")
    kw.setdefault("threads", 1)
    kw.setdefault("max_suboptimal_hits", 1)
    kw.setdefault("max_total_energy", None)
    kw.setdefault("derived_seed", None)
    return intarna.build_command(cfg, **kw)


def test_energy_set_defaults_to_turner04():
    cmd = _build({})
    assert _flag(cmd, "--energyVRNA=") == "--energyVRNA=Turner04"

def test_energy_set_uses_configured_value():
    cmd = _build({"energy_set": "Andronescu07"})
    assert _flag(cmd, "--energyVRNA=") == "--energyVRNA=Andronescu07"

def test_energy_set_explicit_null_falls_back_to_default():
    cmd = _build({"energy_set": None})
    assert _flag(cmd, "--energyVRNA=") == "--energyVRNA=Turner04"

def test_model_defaults_to_X_not_S():
    # IntaRNA's real default is X (seed-extension). The old 'S' default silently
    # selected the legacy Busch-2008 algorithm whenever `model:` was absent.
    cmd = _build({})
    assert _flag(cmd, "--model=") == "--model=X"

def test_model_explicit_null_falls_back_to_X():
    cmd = _build({"model": None})
    assert _flag(cmd, "--model=") == "--model=X"

def test_explicit_null_numeric_falls_back_to_default():
    # The fallbacks must equal IntaRNA's own defaults, so an unset toggle in the UI
    # (which writes null here) is a no-op rather than a silent re-parameterization.
    cmd = _build({"max_loop_size": None, "max_interaction_length": None})
    assert _flag(cmd, "--intLoopMax=") == "--intLoopMax=10"
    assert _flag(cmd, "--intLenMax=") == "--intLenMax=0"

def test_zero_is_preserved_not_treated_as_unset():
    # 0 is a deliberate value here; an `or default` idiom would silently replace it.
    cmd = _build({"max_loop_size": 0})
    assert _flag(cmd, "--intLoopMax=") == "--intLoopMax=0"

def test_build_command_emits_helix_flags():
    cmd = _build({"helix": {"max_bp": 6, "full_energy": True}})
    assert "--helixMaxBP=6" in cmd
    assert "--helixFullE" in cmd

def test_build_command_ignores_legacy_extra_args():
    # extra_args is removed from the schema. A leftover key in a stale config must
    # be inert, not appended: it could otherwise smuggle --qAcc=N past --acc=C and
    # mislabel results written into w_accessibility/.
    cmd = _build({"extra_args": ["--qAcc=N", "--personality=IntaRNAduplex"]})
    assert "--qAcc=N" not in cmd
    assert "--personality=IntaRNAduplex" not in cmd
