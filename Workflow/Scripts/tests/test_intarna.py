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

def test_overlap_explicit_null_falls_back_to_default():
    # An explicit `overlap:` null reaches here as None. A plain .get(k, "B") would pass it
    # through and stringify it into the flag as the literal "--outOverlap=None", which
    # IntaRNA rejects — this is exactly what opt() exists to prevent.
    cmd = []
    intarna._add_output_flags(cmd, {"overlap": None}, max_suboptimal_hits=1)
    assert "--outOverlap=B" in cmd


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


def test_accessibility_shared_flags_still_emitted():
    cmd = []
    intarna._add_accessibility_flags(cmd, {"window": 150, "max_bp_span": 100})
    assert cmd == ["--accW=150", "--accL=100"]

def test_accessibility_per_side_flags_emitted():
    # IntaRNA's own README gives this asymmetric example: global query, local target.
    cmd = []
    intarna._add_accessibility_flags(cmd, {
        "window": None,
        "max_bp_span": None,
        "query_window": 0,
        "query_max_bp_span": 0,
        "target_window": 150,
        "target_max_bp_span": 100,
    })
    assert cmd == ["--qAccW=0", "--qAccL=0", "--tAccW=150", "--tAccL=100"]

def test_accessibility_zero_is_a_real_value_not_a_skip():
    # window 0 = whole sequence (global). A truthiness check would silently drop it
    # and leave IntaRNA on its 150 nt local default — the opposite setting.
    cmd = []
    intarna._add_accessibility_flags(cmd, {"window": 0, "max_bp_span": 0})
    assert cmd == ["--accW=0", "--accL=0"]

def test_accessibility_nulls_are_omitted_entirely():
    # Omission is what hands the decision back to IntaRNA's own defaults (150/100).
    cmd = []
    intarna._add_accessibility_flags(cmd, {
        "window": None, "max_bp_span": None,
        "query_window": None, "query_max_bp_span": None,
        "target_window": None, "target_max_bp_span": None,
        "forbid_lonely_pairs": False, "forbid_gu_at_ends": False,
    })
    assert cmd == []

def test_accessibility_shared_and_per_side_coexist_when_equal():
    # IntaRNA tolerates equal values; only a MISMATCH is rejected, and that is
    # _intarna_config's job, not the emitter's.
    cmd = []
    intarna._add_accessibility_flags(cmd, {"window": 40, "target_window": 40})
    assert cmd == ["--accW=40", "--tAccW=40"]

def test_accessibility_constraint_pills_follow_the_windows():
    cmd = []
    intarna._add_accessibility_flags(cmd, {
        "target_window": 40,
        "forbid_lonely_pairs": True,
        "forbid_gu_at_ends": True,
    })
    assert cmd == ["--tAccW=40", "--accNoLP", "--accNoGUend"]

def test_build_command_emits_per_side_accessibility_flags():
    cmd = _build({"accessibility": {"query_window": 0, "target_window": 150}})
    assert "--qAccW=0" in cmd
    assert "--tAccW=150" in cmd


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
    kw.setdefault("out_max_energy", None)
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

def test_outmaxe_emitted_when_set():
    # Only reached on the wo_accessibility arm, where E == E_hybrid so --outMaxE (a TOTAL
    # energy filter) expresses the hybridization gate exactly.
    cmd = _build({}, out_max_energy=-12.9)
    assert "--outMaxE=-12.9" in cmd

def test_outmaxe_omitted_when_none():
    # w_accessibility passes None: under acc=C, --outMaxE would filter E_hybrid+ED1+ED2,
    # a ~2x stricter bar. tidy_intarna gates on E_hybrid downstream instead.
    cmd = _build({}, out_max_energy=None)
    assert _flag(cmd, "--outMaxE=") is None


SHIPPED_COLUMNS = ("id1,id2,start1,end1,start2,end2,subseqDP,hybridDP,"
                   "E,E_hybrid,ED1,ED2,Pu1,Pu2,seedStart1,seedEnd1,seedE,seedStart2,seedEnd2")


def _csv_cols(cmd):
    return _flag(cmd, "--outCsvCols=").split("=", 1)[1].split(",")


def test_seed_columns_withheld_under_noseed():
    # --noSeed makes IntaRNA write the literal NAN in all five seed columns on every row —
    # ~12% of a panel-scale CSV carrying nothing. Withhold them rather than parse them back.
    cmd = _build({"output": {"columns": SHIPPED_COLUMNS}})        # _build's seed default: None
    assert "--noSeed" in cmd
    assert not [c for c in _csv_cols(cmd) if c.startswith("seed")]
    assert _csv_cols(cmd)[:4] == ["id1", "id2", "start1", "end1"]  # order otherwise intact


def test_seed_columns_kept_when_a_seed_is_enforced():
    cmd = _build({"output": {"columns": SHIPPED_COLUMNS}},
                 derived_seed={"length": 6, "query_range": "2-7"})
    assert _csv_cols(cmd) == SHIPPED_COLUMNS.split(",")
