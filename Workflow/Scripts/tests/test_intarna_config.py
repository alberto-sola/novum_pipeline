from __future__ import annotations
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # Workflow/Scripts

import _intarna_config as cfgmod

validate = cfgmod.validate_intarna_config


def test_module_is_dependency_free():
    # The Snakemake DRIVER imports this at DAG-build time, where none of the arms'
    # scientific deps exist. A third-party import here breaks every snakemake invocation.
    source = (Path(cfgmod.__file__)).read_text()
    imports = [l.strip() for l in source.splitlines()
               if l.startswith("import ") or l.startswith("from ")]
    assert imports == [], f"_intarna_config must stay stdlib-only, found: {imports}"


# --- opt(): explicit null must fall back, but 0 must survive ---

def test_opt_falls_back_on_missing_and_null():
    assert cfgmod.opt({}, "k", 7) == 7
    assert cfgmod.opt({"k": None}, "k", 7) == 7

def test_opt_preserves_zero():
    assert cfgmod.opt({"k": 0}, "k", 7) == 0


# --- derive_seed_from_string ---

def test_derive_seed_width_is_inclusive():
    assert cfgmod.derive_seed_from_string("2,7") == {"length": 6, "query_range": "2-7"}

def test_derive_seed_none_passes_through():
    assert cfgmod.derive_seed_from_string(None) is None


# --- (a) seed must fit inside one helix under model B ---

def test_validate_raises_when_seed_exceeds_default_helix_max_bp():
    # seed 1,12 -> seedBP=12 > IntaRNA's default helixMaxBP of 10
    with pytest.raises(ValueError, match="helix"):
        validate({"model": "B"}, "1,12")

def test_validate_raises_when_seed_exceeds_explicit_helix_max_bp():
    # The likelier path: the stock seed 2,7 (seedBP = 7-2+1 = 6) against a lowered cap.
    with pytest.raises(ValueError, match="helix"):
        validate({"model": "B", "helix": {"max_bp": 2}}, "2,7")

def test_validate_allows_seed_that_fits_the_helix():
    assert validate({"model": "B", "helix": {"max_bp": 10}}, "2,7") == []

def test_validate_seed_helix_check_is_scoped_to_model_B():
    # helixMaxBP is inert under X/S/P — verified: identical energies with a maximally
    # restrictive helix block. So a seed too wide for max_bp must NOT raise here, the way
    # it does under model B. It does still earn the (c) warning, which is correct advice.
    warnings = validate({"model": "X", "helix": {"max_bp": 2}}, "2,7")
    assert len(warnings) == 1
    assert "will be ignored" in warnings[0]

def test_validate_no_seed_cannot_conflict():
    assert validate({"model": "B", "helix": {"max_bp": 2}}, None) == []

def test_validate_rejects_a_malformed_seed_under_every_model():
    # The seed string is only parsed here, so a bad one must fail at DAG build whatever
    # the model — not just on the model-B path that happens to need seedBP.
    for model in ("X", "S", "B", "P"):
        with pytest.raises(ValueError, match="seed"):
            validate({"model": model}, "banana")
        with pytest.raises(ValueError, match="seed"):
            validate({"model": model}, "1,30")  # seedBP=30, outside [2,20]


# --- (b) energy_set enum ---

def test_validate_raises_on_unknown_energy_set():
    with pytest.raises(ValueError, match="energy_set"):
        validate({"energy_set": "Turner4"}, None)

@pytest.mark.parametrize("name", ["Turner99", "Turner04", "Andronescu07"])
def test_validate_accepts_every_named_energy_set(name):
    assert validate({"energy_set": name}, None) == []

def test_validate_accepts_missing_and_null_energy_set():
    assert validate({}, None) == []
    assert validate({"energy_set": None}, None) == []


# --- (c) helix set but model is not B ---

def test_validate_warns_when_helix_set_under_wrong_model():
    warnings = validate({"model": "X", "helix": {"min_bp": 3}}, None)
    assert len(warnings) == 1
    assert "helix" in warnings[0] and "'X'" in warnings[0]

def test_validate_helix_warning_counts_zero_as_set():
    # 0 is a real value for max_internal_loop; it must not be treated as "unset".
    assert len(validate({"model": "X", "helix": {"max_internal_loop": 0}}, None)) == 1

def test_validate_helix_warning_ignores_false_full_energy():
    assert validate({"model": "X", "helix": {"full_energy": False}}, None) == []

def test_validate_silent_for_empty_helix_block():
    assert validate({"model": "X", "helix": {}}, None) == []
    assert validate({"model": "X"}, None) == []

def test_validate_no_helix_warning_under_model_B():
    assert validate({"model": "B", "helix": {"min_bp": 3}}, None) == []


# --- removed top-level keys ---

def test_reject_removed_keys_raises_on_max_total_energy():
    with pytest.raises(ValueError, match="max_total_energy has been removed"):
        cfgmod.reject_removed_keys({"max_total_energy": -18})

def test_reject_removed_keys_raises_even_when_the_value_is_null():
    # Nulling the key is the likeliest "I disabled it" migration attempt; silence there
    # would leave the user believing they had turned the cutoff off on both arms.
    with pytest.raises(ValueError, match="max_total_energy"):
        cfgmod.reject_removed_keys({"max_total_energy": None})

def test_reject_removed_keys_message_names_both_replacements():
    with pytest.raises(ValueError) as exc:
        cfgmod.reject_removed_keys({"max_total_energy": -18})
    msg = str(exc.value)
    assert "rnahybrid:" in msg and "intarna:" in msg
    assert "-18" in msg and "-12.9" in msg

def test_reject_removed_keys_silent_on_clean_config():
    assert cfgmod.reject_removed_keys({"threads": 16}) is None
    assert cfgmod.reject_removed_keys({}) is None
    assert cfgmod.reject_removed_keys(None) is None


# --- (d) the two arms must gate at comparable bars ---

def test_offset_constant_is_the_measured_value():
    assert cfgmod.INTARNA_RNAHYBRID_ENERGY_OFFSET == 5.10

def test_energy_warning_silent_when_the_bars_are_comparable():
    assert validate({"max_hybrid_energy": -12.9}, None,
                    rnahybrid_max_hybrid_energy=-18) == []

def test_energy_warning_fires_on_non_comparable_bars():
    # -18 on BOTH arms is the pre-refactor bug: it demands ~2x the literature threshold
    # of IntaRNA. That must be loud, not silently accepted.
    warnings = validate({"max_hybrid_energy": -18.0}, None,
                        rnahybrid_max_hybrid_energy=-18)
    assert len(warnings) == 1
    assert "non-comparable" in warnings[0]
    assert "-12.90" in warnings[0]          # names the value the user probably wanted

def test_energy_warning_silent_when_both_cutoffs_are_null():
    # Both unset is a deliberate "no gating anywhere" — nothing to warn about.
    assert validate({"max_hybrid_energy": None}, None) == []

def test_energy_warning_fires_when_only_intarna_cutoff_is_null():
    # One arm gated, the other not, is a BIGGER mismatch than two numeric cutoffs
    # merely drifting apart — and exactly where a half-finished hand migration off
    # max_total_energy lands. It must not be silent.
    warnings = validate({"max_hybrid_energy": None}, None,
                        rnahybrid_max_hybrid_energy=-18)
    assert len(warnings) == 1
    assert "intarna.max_hybrid_energy is null" in warnings[0]
    assert "-18" in warnings[0]

def test_energy_warning_fires_when_only_rnahybrid_cutoff_is_null():
    warnings = validate({"max_hybrid_energy": -12.9}, None)
    assert len(warnings) == 1
    assert "rnahybrid.max_hybrid_energy is null" in warnings[0]
    assert "-12.9" in warnings[0]

def test_energy_warning_tolerance_is_one_kcal():
    # A usability threshold, not a confidence bound — the median offset is pinned far
    # tighter than this, so 1 kcal/mol of drift is a user choice, not measurement noise.
    assert validate({"max_hybrid_energy": -13.8}, None,
                    rnahybrid_max_hybrid_energy=-18) == []           # 0.9 off — quiet
    assert len(validate({"max_hybrid_energy": -14.0}, None,
                        rnahybrid_max_hybrid_energy=-18)) == 1       # 1.1 off — warns

def test_energy_warning_coexists_with_the_helix_warning():
    warnings = validate({"model": "X", "helix": {"min_bp": 3},
                         "max_hybrid_energy": -18.0}, None,
                        rnahybrid_max_hybrid_energy=-18)
    assert len(warnings) == 2


# --- (e) w_accessibility withholds --outMaxE, so only the per-pair cap bounds the CSV ---

def test_unbounded_output_warning_fires_on_accessibility_with_null_cap():
    warnings = validate({}, None, max_suboptimal_hits=None, accessibility_on=True)
    assert len(warnings) == 1
    assert str(cfgmod.INTARNA_MAX_OUTNUMBER) in warnings[0]

def test_unbounded_output_warning_says_the_depth_knob_cannot_bound_it():
    # accessibility_search_depth is a FLOOR on --outNumber: with a null cap the max()
    # that would apply it is skipped entirely, so advising it here would be wrong.
    warnings = validate({"accessibility_search_depth": 20}, None,
                        max_suboptimal_hits=None, accessibility_on=True)
    assert "accessibility_search_depth cannot" in warnings[0]

def test_unbounded_output_warning_silent_once_the_cap_is_set():
    assert validate({}, None, max_suboptimal_hits=1, accessibility_on=True) == []

def test_unbounded_output_warning_silent_without_the_accessibility_arm():
    # wo_accessibility carries --outMaxE, so the intermediate is already bounded there.
    assert validate({}, None, max_suboptimal_hits=None, accessibility_on=False) == []

def test_unbounded_output_warning_defaults_off_for_callers_that_omit_it():
    assert validate({}, None) == []


# --- (f) IntaRNA refuses a shared accessibility window beside a differing per-side one ---

@pytest.mark.parametrize("acc,shared_key,side_key", [
    ({"window": 0, "query_window": 150}, "window", "query_window"),
    ({"window": 0, "target_window": 150}, "window", "target_window"),
    ({"max_bp_span": 0, "query_max_bp_span": 100}, "max_bp_span", "query_max_bp_span"),
    ({"max_bp_span": 0, "target_max_bp_span": 100}, "max_bp_span", "target_max_bp_span"),
])
def test_accessibility_conflict_is_fatal(acc, shared_key, side_key):
    # Measured on IntaRNA 3.4.1: `--accW=0 --accL=0 --tAccW=40 --tAccL=40` exits with
    # "# ERROR : --accW and --tAccW are set to different values". The flags do not
    # layer and neither wins, so this must never reach the tool.
    with pytest.raises(ValueError) as excinfo:
        validate({"accessibility": acc}, None)
    message = str(excinfo.value)
    assert shared_key in message
    assert side_key in message

def test_accessibility_equal_values_are_allowed():
    # IntaRNA itself tolerates the redundant-but-consistent form.
    assert validate({"accessibility": {"window": 40, "target_window": 40}}, None) == []

def test_accessibility_per_side_alone_is_the_normal_override():
    assert validate({"accessibility": {
        "query_window": 0, "query_max_bp_span": 0,
        "target_window": 150, "target_max_bp_span": 100,
    }}, None) == []

def test_accessibility_shared_alone_is_allowed():
    assert validate({"accessibility": {"window": 0, "max_bp_span": 0}}, None) == []

def test_accessibility_zero_versus_null_is_not_a_conflict():
    # `window: 0` with an unset per-side key is an override-free config, not a mismatch.
    assert validate({"accessibility": {"window": 0, "query_window": None}}, None) == []

def test_accessibility_block_absent_is_allowed():
    assert validate({}, None) == []
