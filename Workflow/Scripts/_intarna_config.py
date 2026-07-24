"""IntaRNA config schema: defaults, seed derivation, and DAG-build-time validation.

Imported by BOTH `intarna.py` (at rule-execution time, inside the IntaRNA conda env) and
the Snakefile (at DAG-build time, in the Snakemake driver env). That second consumer is
why this module must stay **dependency-free** — stdlib only, no pandas, no third-party
imports — exactly like `_common.py`. An import added here runs in the driver env, which
ships none of the arms' scientific dependencies.

Keeping this separate from `intarna.py` means the driver never imports a script that ends
in `run_from_snakemake(snakemake)`, so nothing here depends on that call staying guarded.
"""

INTARNA_ENERGY_SETS = ("Turner99", "Turner04", "Andronescu07")
INTARNA_DEFAULT_ENERGY_SET = "Turner04"

# IntaRNA's own defaults, mirrored here because validate_intarna_config reproduces two
# of its constraints.
INTARNA_DEFAULT_MODEL = "X"
INTARNA_DEFAULT_HELIX_MAX_BP = 10

# --outNumber's hard ceiling: IntaRNA accepts [0,1000], where 0 means "report nothing", so
# 1000 is how `max_suboptimal_hits: null` ("all hits") is spelled. intarna.py clamps to it;
# warning (e) quotes it.
INTARNA_MAX_OUTNUMBER = 1000

# Ours, not IntaRNA's: the --outNumber floor applied on w_accessibility only, where acc=C
# ranks by total E while the gate is on E_hybrid.
INTARNA_DEFAULT_ACCESSIBILITY_SEARCH_DEPTH = 20

# Measured median of (IntaRNA E_hybrid - RNAhybrid mfe) over matched, unfiltered acc=N
# hits on 600 E. coli CDS (n=1122): additive, regression slope +0.0019 against mfe, sd
# 1.06. Independently reproduced at 5.12. IntaRNA's E_hybrid carries duplex initiation
# (~+4.1 under Turner04) plus terminal-AU/dangling-end terms; RNAhybrid's mfe carries
# neither. See docs/superpowers/specs/2026-07-23-cross-arm-energy-gate-design.md.
INTARNA_RNAHYBRID_ENERGY_OFFSET = 5.10

# Warning (d) fires past this much drift from the expected relationship. Chosen as a
# usability threshold, not a confidence bound: the median offset is known far more tightly
# than this (SE ~0.03 at n=1122), so anything past 1 kcal/mol is a deliberate choice by the
# user rather than measurement noise.
ENERGY_OFFSET_TOLERANCE = 1.0

# Top-level keys that nothing reads any more, mapped to the message that explains the
# hand migration. Keeping such a key inert-but-present is the exact failure this refactor
# exists to end, so its presence is fatal rather than ignored.
REMOVED_TOP_LEVEL_KEYS = {
    "max_total_energy": (
        "max_total_energy has been removed. It filtered pure hybridization energy on the "
        "RNAhybrid arm but TOTAL interaction energy (E_hybrid + ED1 + ED2) on the IntaRNA "
        "arm, so one number meant two different bars — roughly a 2x stricter one for "
        "IntaRNA. There is no single correct automatic migration; replace it with the "
        "per-arm keys:\n"
        "    rnahybrid:\n"
        "      max_hybrid_energy: -18      # -e ; the literature MFE threshold\n"
        "    intarna:\n"
        "      max_hybrid_energy: -12.9    # -18 + 5.10 measured convention offset"
    ),
}


#----- Raises on a top-level key nothing reads any more. Called from the Snakefile at DAG-build time, before any rule runs -----#
def reject_removed_keys(config):
    for key, message in REMOVED_TOP_LEVEL_KEYS.items():
        if key in (config or {}):
            raise ValueError(message)


#----- cfg.get() that survives an explicit null: YAML `key:` parses to None, which a plain .get(key, default) would pass through and stringify into the flag as "None" -----#
def opt(cfg, key, default):
    value = cfg.get(key)
    return default if value is None else value


#----- Parses the shared top-level seed "x,y" into IntaRNA's seedBP+seedQRange shape; raises if the width is outside [2,20] -----#
def derive_seed_from_string(seed_str):
    if seed_str is None:
        return None
    try:
        parts = str(seed_str).split(",")
        if len(parts) != 2:
            raise ValueError
        x, y = int(parts[0].strip()), int(parts[1].strip())
    except ValueError:
        raise ValueError(f"top-level seed must be 'x,y' integers, got: {seed_str!r}")
    bp = y - x + 1
    if not (2 <= bp <= 20):
        raise ValueError(
            f"top-level seed '{seed_str}' gives seedBP={bp}, which is outside IntaRNA's [2,20] range"
        )
    return {"length": bp, "query_range": f"{x}-{y}"}


#----- (d) The two arms' cutoffs must sit ~one convention offset apart, or the consensus silently intersects non-comparable bars -----#
def _energy_comparability_warnings(cfg, rnahybrid_max_hybrid_energy):
    intarna_cutoff = cfg.get("max_hybrid_energy")

    # Both unset: a deliberate "no gating anywhere" — nothing to warn about.
    if intarna_cutoff is None and rnahybrid_max_hybrid_energy is None:
        return []

    # Exactly one unset is the more asymmetric case — a bigger mismatch than two
    # numeric cutoffs merely drifting apart, and exactly where a half-finished hand
    # migration off max_total_energy lands.
    if intarna_cutoff is None or rnahybrid_max_hybrid_energy is None:
        null_key, null_arm, set_key, set_value = (
            ("intarna", "IntaRNA", "rnahybrid", rnahybrid_max_hybrid_energy)
            if intarna_cutoff is None
            else ("rnahybrid", "RNAhybrid", "intarna", intarna_cutoff)
        )
        return [
            f"{null_key}.max_hybrid_energy is null ({null_arm} arm ungated on "
            f"hybridization energy) while {set_key}.max_hybrid_energy is {set_value}. "
            f"The consensus intersects a gated arm with an ungated one."
        ]

    expected = rnahybrid_max_hybrid_energy + INTARNA_RNAHYBRID_ENERGY_OFFSET
    drift = abs(intarna_cutoff - expected)
    if drift <= ENERGY_OFFSET_TOLERANCE:
        return []
    return [
        f"intarna.max_hybrid_energy is {intarna_cutoff} while rnahybrid.max_hybrid_energy "
        f"is {rnahybrid_max_hybrid_energy}, which is {expected:.2f} on IntaRNA's E_hybrid "
        f"scale (measured offset {INTARNA_RNAHYBRID_ENERGY_OFFSET:+.2f} kcal/mol). The two "
        f"arms are gating at non-comparable bars, {drift:.2f} kcal/mol apart, which biases "
        f"the consensus toward whichever arm is looser."
    ]


#----- (e) w_accessibility withholds --outMaxE, so the per-pair cap is the only thing bounding the intermediate CSV. accessibility_search_depth cannot stand in: it is a FLOOR on --outNumber, and a null cap skips the max() that would apply it -----#
def _unbounded_output_warnings(max_suboptimal_hits, accessibility_on):
    if not accessibility_on or max_suboptimal_hits is not None:
        return []
    return [
        "intarna: w_accessibility withholds --outMaxE (its energy gate lives downstream "
        "in tidy_intarna.py), so with max_suboptimal_hits: null IntaRNA writes up to "
        f"{INTARNA_MAX_OUTNUMBER} hits/pair, bounded only by its own --outMaxE default "
        "of 0. Set max_suboptimal_hits to bound intarna_output.csv; "
        "intarna.accessibility_search_depth cannot — it is a floor on --outNumber, which "
        "a null cap ignores entirely."
    ]


#----- DAG-build-time config validation: raises ValueError on a config IntaRNA would die on, so the run fails before the hour-long RNAhybrid arm starts; returns advisory warnings -----#
def validate_intarna_config(cfg, seed_str, rnahybrid_max_hybrid_energy=None,
                            max_suboptimal_hits=None, accessibility_on=False):
    cfg = cfg or {}
    helix = cfg.get("helix") or {}
    model = opt(cfg, "model", INTARNA_DEFAULT_MODEL)

    # IntaRNA reads an unknown energy set as a parameter-FILE path, fails to load it, and
    # exits 255 asking the user to file a bug report. Reject it here instead.
    energy_set = opt(cfg, "energy_set", INTARNA_DEFAULT_ENERGY_SET)
    if energy_set not in INTARNA_ENERGY_SETS:
        raise ValueError(
            f"intarna.energy_set must be one of {list(INTARNA_ENERGY_SETS)}, got {energy_set!r}"
        )

    # The only validation of the top-level seed string, so it must run for every model.
    derived_seed = derive_seed_from_string(seed_str)

    warnings = _energy_comparability_warnings(cfg, rnahybrid_max_hybrid_energy)
    warnings += _unbounded_output_warnings(max_suboptimal_hits, accessibility_on)

    if model != "B":
        # `0` is a legitimate value (max_internal_loop: 0 = pure stacks), so compare
        # against None/False rather than testing truthiness.
        inert = sorted(k for k, v in helix.items() if v is not None and v is not False)
        if inert:
            warnings.append(
                f"intarna.helix {inert} is set but intarna.model is {model!r}: helix "
                f"parameters apply only to model B and will be ignored."
            )
        return warnings

    # Under model B the seed must fit inside a single helix or IntaRNA hard-errors —
    # reachable both by widening the seed and by lowering the cap beneath it.
    if derived_seed is not None:
        max_bp = opt(helix, "max_bp", INTARNA_DEFAULT_HELIX_MAX_BP)
        if derived_seed["length"] > max_bp:
            raise ValueError(
                f"seed {seed_str!r} needs {derived_seed['length']} base pairs, but "
                f"intarna.helix.max_bp is {max_bp}. Under intarna.model: B the seed must "
                f"fit inside one helix. Raise intarna.helix.max_bp to at least "
                f"{derived_seed['length']}, or shorten the top-level seed."
            )

    return warnings
