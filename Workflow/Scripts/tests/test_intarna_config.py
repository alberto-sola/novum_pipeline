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
