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


#----- DAG-build-time config validation: raises ValueError on a config IntaRNA would die on, so the run fails before the hour-long RNAhybrid arm starts; returns advisory warnings -----#
def validate_intarna_config(cfg, seed_str):
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

    if model != "B":
        # `0` is a legitimate value (max_internal_loop: 0 = pure stacks), so compare
        # against None/False rather than testing truthiness.
        inert = sorted(k for k, v in helix.items() if v is not None and v is not False)
        if not inert:
            return []
        return [
            f"intarna.helix {inert} is set but intarna.model is {model!r}: helix "
            f"parameters apply only to model B and will be ignored."
        ]

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

    return []
