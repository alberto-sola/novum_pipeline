"""Guards the UI's hand-maintained mirrors of the Python/YAML config surface.

Three invariants:

1. INTARNA_BLOCKS key order is the order safe_dump(sort_keys=False) writes back to
   Config/config.yaml, so the two must not drift — and every key rendered from a block must
   exist in INITIAL_CONFIG or hydrate throws.
2. The JS constants whose comments claim to mirror _intarna_config.py really do. There is no
   shared runtime (plain <script> tags, no Node build), so the duplication is unavoidable;
   its drifting silently is not.
3. The rnacalibrate block round-trips through the REAL hydrateConfig -> configToObject pair.
   F1 broke this: the emit line, the hydrate line, and INITIAL_CONFIG's default can each go
   missing independently, and only executing the real JS catches all three.

Most of this file bracket-matches data.js/config.js rather than running a JS engine — they
are plain `window.X = {...}` assignments and the repo has no Node *build* dependency. The
round-trip tests are the exception: they shell out to `node` (see _js_roundtrip.js), because
comparing source text cannot catch a behavioural break.
"""

from __future__ import annotations
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parent.parent.parent  # repo root
sys.path.insert(0, str(ROOT / "Workflow" / "Scripts"))

import _intarna_config  # stdlib-only by contract, so importable from any env

NODE = shutil.which("node")
JS_ROUNDTRIP_HARNESS = Path(__file__).resolve().parent / "_js_roundtrip.js"

# Every intarna sub-block that is both a UI block and a nested mapping in config.yaml.
# `top` is excluded here — its keys live at the intarna level, so it has its own pair of
# tests below rather than a nested mapping to compare against.
BLOCKS = ["helix", "seed", "accessibility", "output"]

# The intarna-level keys config.js writes by hand, outside INTARNA_BLOCKS.top.
INTARNA_TOP_LITERALS = ("accessibility_variant", "prediction_mode", "model", "energy_set")


def _data_js_source():
    return (ROOT / "UI" / "src" / "data.js").read_text()


def _config_js_source():
    return (ROOT / "UI" / "src" / "config.js").read_text()


#----- Text between the bracket at `pos` and its match, by nesting depth -----#
def _body_at(source, pos, open_ch, close_ch):
    depth, i = 1, pos + 1
    while depth:
        if source[i] == open_ch:
            depth += 1
        elif source[i] == close_ch:
            depth -= 1
        i += 1
    return source[pos + 1:i - 1]


#----- Body of a top-level `window.<name> = { … }` / `= [ … ]` assignment -----#
def _assignment_body(name, open_ch="[", close_ch="]"):
    source = _data_js_source()
    match = re.search(rf"window\.{re.escape(name)}\s*=\s*{re.escape(open_ch)}", source)
    assert match, f"could not locate window.{name} in data.js"
    return _body_at(source, match.end() - 1, open_ch, close_ch)


#----- Body of a `<name>: { … }` / `[ … ]` entry nested in an already-extracted body. Scoped rather than searched from the file start because INITIAL_CONFIG's top-level `seed:` would otherwise shadow `intarna.seed` -----#
def _nested_body(container, name, open_ch, close_ch):
    match = re.search(rf"\b{re.escape(name)}:\s*{re.escape(open_ch)}", container)
    assert match, f"could not locate {name}"
    return _body_at(container, match.end() - 1, open_ch, close_ch)


def _intarna_initial_config():
    return _nested_body(_assignment_body("INITIAL_CONFIG", "{", "}"), "intarna", "{", "}")


def _block_entries(name):
    body = _nested_body(_assignment_body("INTARNA_BLOCKS", "{", "}"), name, "[", "]")
    # Kind letter is \w, not [npc]: an entry with a drifted kind must still be picked up, or
    # the order assertion could pass against a block that has actually diverged.
    return re.findall(r'\["([a-z_]+)",\s*"(\w+)"\]', body)


def _block_keys(name):
    return [key for key, _kind in _block_entries(name)]


def _initial_config_keys(name):
    body = _nested_body(_intarna_initial_config(), name, "{", "}")
    return re.findall(r"(?m)^\s*(\w+):", body)


def _config_yaml_keys(name):
    cfg = yaml.safe_load((ROOT / "Config" / "config.yaml").read_text())
    return list(cfg["intarna"][name].keys())


@pytest.mark.parametrize("block", BLOCKS)
def test_ui_block_order_matches_config_yaml(block):
    assert _block_keys(block) == _config_yaml_keys(block)


@pytest.mark.parametrize("block", BLOCKS)
def test_initial_config_has_every_editable_block_key(block):
    # Guards UI/src/config.js (`state[k].value`): a key in INTARNA_BLOCKS but missing from
    # INITIAL_CONFIG throws during hydrate, producing a blank UI on launch. Kind "c" is
    # exempt — a constant config.js deliberately never hydrates back.
    editable = [key for key, kind in _block_entries(block) if kind != "c"]
    missing = [k for k in editable if k not in _initial_config_keys(block)]
    assert not missing, (
        f"INITIAL_CONFIG.intarna.{block} is missing keys present in "
        f"INTARNA_BLOCKS.{block}: {missing}"
    )


def test_ui_top_block_matches_config_yaml_scalars():
    # configToObject rebuilds obj.intarna wholesale from INTARNA_TOP_LITERALS plus
    # INTARNA_BLOCKS.top, and save_config dumps that object over the file — so an
    # intarna-level key in neither is DELETED the first time the editor saves, silently,
    # because absent and null mean the same thing everywhere downstream.
    cfg = yaml.safe_load((ROOT / "Config" / "config.yaml").read_text())
    scalars = [key for key, value in cfg["intarna"].items()
               if not isinstance(value, dict) and key not in INTARNA_TOP_LITERALS]
    assert _block_keys("top") == scalars


def test_ui_rnacalibrate_block_matches_config_yaml_scalars():
    # The rnacalibrate analogue of test_ui_top_block_matches_config_yaml_scalars above;
    # config.js hand-lists `obj.rnacalibrate` (there is no RNACALIBRATE_BLOCKS table).
    # This is exactly how `rng_seed` went missing (F1), after which the next UI-launched
    # run went unpinned with no error or warning.
    source = _config_js_source()
    match = re.search(r"obj\.rnacalibrate\s*=\s*\{", source)
    assert match, "could not locate obj.rnacalibrate in config.js"
    body = _body_at(source, match.end() - 1, "{", "}")
    js_keys = re.findall(r"(?m)^\s*(\w+):", body)

    cfg = yaml.safe_load((ROOT / "Config" / "config.yaml").read_text())
    assert js_keys == list(cfg["rnacalibrate"].keys())


def _run_js_roundtrip(raw):
    # Hydrates then re-emits `raw` through the real hydrateConfig -> configToObject pair
    # (see _js_roundtrip.js). Returns the harness's {"hydrated": ..., "obj": ...} dict.
    result = subprocess.run(
        [NODE, str(JS_ROUNDTRIP_HARNESS),
         str(ROOT / "UI" / "src" / "data.js"), str(ROOT / "UI" / "src" / "config.js")],
        input=json.dumps({"raw": raw}), capture_output=True, text=True,
    )
    assert result.returncode == 0, (
        f"node round-trip harness exited {result.returncode}; stderr:\n{result.stderr}"
    )
    return json.loads(result.stdout)


@pytest.mark.skipif(NODE is None, reason="node not found on PATH")
@pytest.mark.parametrize("rng_seed", [1, 0, None])
def test_rnacalibrate_block_round_trips_through_the_real_js(rng_seed):
    # Catches the three breaks the emit-direction test above cannot: a missing hydrate line
    # (F1 verbatim — the disk value falls back to the INITIAL_CONFIG default and re-saves as
    # it), emitting the raw {set, value} object unwrapped, and the key vanishing from
    # INITIAL_CONFIG (hydrateConfig throws, blanking the UI).
    # Asserting the whole block costs nothing but only extends cover to keys whose shipped
    # value DIFFERS from its default — calibration_variant and rng_seed, verified by
    # mutation. Where they coincide (k, max_target_length, randomize_targets) a dropped
    # hydrate line stays invisible; closing that needs the per-key table config.js already
    # uses for intarna.
    # All three params matter: null means "don't pin" not "unset", 0 guards falsy-collapse,
    # and null alone cannot detect a deleted hydrate line (its fallback emits null too).
    cfg = yaml.safe_load((ROOT / "Config" / "config.yaml").read_text())
    cfg["rnacalibrate"]["rng_seed"] = rng_seed
    out = _run_js_roundtrip(cfg)
    assert out["obj"]["rnacalibrate"] == cfg["rnacalibrate"]


def test_initial_config_has_every_top_block_key():
    # Same hydrate guard as the per-block test, for the keys that sit at the intarna level.
    # Four spaces is that level exactly: nested blocks indent their own keys further.
    keys = re.findall(r"(?m)^\s{4}(\w+):", _intarna_initial_config())
    missing = [key for key, kind in _block_entries("top") if kind != "c" and key not in keys]
    assert not missing, f"INITIAL_CONFIG.intarna is missing: {missing}"


def test_per_side_accessibility_keys_are_present():
    keys = _config_yaml_keys("accessibility")
    for key in ("query_window", "query_max_bp_span",
                "target_window", "target_max_bp_span"):
        assert key in keys, f"{key} missing from Config/config.yaml"


def test_js_accessibility_conflicts_mirror_python():
    # data.js's copy drives the UI's live "IntaRNA will abort" note. A pair added on the
    # Python side alone leaves the UI silent about a config that is still fatal.
    body = _assignment_body("ACCESSIBILITY_WINDOW_CONFLICTS")
    js_pairs = [tuple(p) for p in re.findall(r'\["(\w+)",\s*"(\w+)"\]', body)]
    assert js_pairs == list(_intarna_config.ACCESSIBILITY_WINDOW_CONFLICTS)


def test_js_accessibility_sides_cover_the_python_flag_table():
    # The --qAccW/--qAccL/--tAccW/--tAccL labels the per-side fields display must name the
    # flags intarna.py actually emits for those keys.
    body = _assignment_body("ACCESSIBILITY_SIDES")
    flag_for = dict(_intarna_config.ACCESSIBILITY_WINDOW_FLAGS)
    sides = re.findall(
        r'side:\s*"(\w+)".*?wFlag:\s*"(--\w+)",\s*lFlag:\s*"(--\w+)"', body
    )
    assert sides, "could not parse ACCESSIBILITY_SIDES"
    for side, w_flag, l_flag in sides:
        assert flag_for[f"{side}_window"] == w_flag
        assert flag_for[f"{side}_max_bp_span"] == l_flag


def test_js_energy_sets_mirror_python():
    # The dropdown is [{value, label}]; only `value` reaches config.yaml and fatal (b).
    body = _assignment_body("INTARNA_ENERGY_SETS")
    js_sets = tuple(re.findall(r'value:\s*"(\w+)"', body))
    assert js_sets == _intarna_config.INTARNA_ENERGY_SETS


def test_js_default_helix_max_bp_mirrors_python():
    match = re.search(r"INTARNA_DEFAULT_HELIX_MAX_BP\s*=\s*(\d+)", _data_js_source())
    assert match, "could not locate INTARNA_DEFAULT_HELIX_MAX_BP in data.js"
    assert int(match.group(1)) == _intarna_config.INTARNA_DEFAULT_HELIX_MAX_BP
