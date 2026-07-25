import subprocess
import sys

from _common import which_required, ensure_parent
from _intarna_config import (
    ACCESSIBILITY_WINDOW_FLAGS,
    INTARNA_DEFAULT_ENERGY_SET,
    INTARNA_DEFAULT_MODEL,
    INTARNA_MAX_OUTNUMBER,
    derive_seed_from_string,
    opt,
)


#----- Emits `--flag=value` per (config key, flag) that is set. Guard is `is not None`, never
#      truthiness: `0` is a real value here (accessibility window 0 = whole sequence) -----#
def _add_value_flags(cmd, block, table):
    if not block:
        return
    for key, flag in table:
        value = block.get(key)
        if value is not None:
            cmd.append(f"{flag}={value}")


#----- Emits the bare `--flag` for each (config key, flag) the block sets truthy -----#
def _add_bool_flags(cmd, block, table):
    if not block:
        return
    for key, flag in table:
        if block.get(key):
            cmd.append(flag)


# Flag tables, in Config/config.yaml order. Every group is emitted unconditionally — IntaRNA
# ignores the inapplicable ones cleanly (helix outside model B, accessibility under --acc=N),
# and validate_intarna_config reports the mismatch instead of the emitter withholding flags.
_SEED_VALUE_FLAGS = (
    ("max_energy", "--seedMaxE"),
    ("max_hybrid_energy", "--seedMaxEhybrid"),
    ("min_unpaired_probability", "--seedMinPu"),
    ("target_range", "--seedTRange"),
)
_SEED_BOOL_FLAGS = (
    ("forbid_gu", "--seedNoGU"),
    ("forbid_gu_at_ends", "--seedNoGUend"),
    ("report_best_only", "--outBestSeedOnly"),
)
_ACCESSIBILITY_BOOL_FLAGS = (
    ("forbid_lonely_pairs", "--accNoLP"),
    ("forbid_gu_at_ends", "--accNoGUend"),
)
_HELIX_VALUE_FLAGS = (
    ("min_bp", "--helixMinBP"),
    ("max_bp", "--helixMaxBP"),
    ("max_internal_loop", "--helixMaxIL"),
    ("min_unpaired_probability", "--helixMinPu"),
    ("max_energy", "--helixMaxE"),
)
_HELIX_BOOL_FLAGS = (("full_energy", "--helixFullE"),)
_OUTPUT_VALUE_FLAGS = (
    ("max_delta_energy", "--outDeltaE"),
    ("min_unpaired_probability", "--outMinPu"),
)
_OUTPUT_BOOL_FLAGS = (
    ("forbid_lonely_pairs", "--outNoLP"),
    ("forbid_gu_at_ends", "--outNoGUend"),
)


#----- Seed enforcement is inherited from the top-level `seed`: absent → --noSeed. When present, --seedBP/--seedQRange are derived and the sub-keys are independent overrides -----#
def _add_seed_flags(cmd, seed, derived_seed=None):
    if derived_seed is None:
        cmd.append("--noSeed")
        return

    cmd.append(f"--seedBP={derived_seed['length']}")
    cmd.append(f"--seedQRange={derived_seed['query_range']}")
    # --seedMaxUP: explicit override wins, else RNAhybrid-parity 0
    cmd.append(f"--seedMaxUP={opt(seed, 'max_unpaired_bases', 0)}")

    _add_value_flags(cmd, seed, _SEED_VALUE_FLAGS)
    _add_bool_flags(cmd, seed, _SEED_BOOL_FLAGS)


#----- Accessibility window/span flags plus --accNoLP/--accNoGUend. The per-side flags override
#      the shared pair but do NOT layer: IntaRNA hard-errors on a shared value beside a differing
#      per-side one, so that conflict is fatal in validate_intarna_config, not resolved here -----#
def _add_accessibility_flags(cmd, acc):
    _add_value_flags(cmd, acc, ACCESSIBILITY_WINDOW_FLAGS)
    _add_bool_flags(cmd, acc, _ACCESSIBILITY_BOOL_FLAGS)


#----- Appends --helix* flags from the helix block -----#
def _add_helix_flags(cmd, helix):
    _add_value_flags(cmd, helix, _HELIX_VALUE_FLAGS)
    _add_bool_flags(cmd, helix, _HELIX_BOOL_FLAGS)


#----- --outNumber accepts [0,1000]; a null cap means "all hits", spelled as the ceiling -----#
def _resolve_outnumber(max_suboptimal_hits):
    if max_suboptimal_hits is None:
        return INTARNA_MAX_OUTNUMBER
    if max_suboptimal_hits > INTARNA_MAX_OUTNUMBER:
        print(
            f"intarna: max_suboptimal_hits={max_suboptimal_hits} exceeds IntaRNA's "
            f"--outNumber ceiling; clamping to {INTARNA_MAX_OUTNUMBER}.",
            file=sys.stderr,
        )
        return INTARNA_MAX_OUTNUMBER
    return max_suboptimal_hits


#----- --out* filtering flags. `out_max_energy` is IntaRNA's TOTAL-energy bound (--outMaxE), supplied by the Snakefile only where it equals the hybridization gate (acc=N) -----#
def _add_output_flags(cmd, out, max_suboptimal_hits=None, out_max_energy=None):
    if out_max_energy is not None:
        cmd.append(f"--outMaxE={out_max_energy}")
    _add_value_flags(cmd, out, _OUTPUT_VALUE_FLAGS)
    cmd.append(f"--outNumber={_resolve_outnumber(max_suboptimal_hits)}")
    cmd.append(f"--outOverlap={opt(out, 'overlap', 'B')}")
    _add_bool_flags(cmd, out, _OUTPUT_BOOL_FLAGS)
    # Truthiness, not `is not None`: an empty `columns:` must fall back to IntaRNA's own
    # column set rather than emit a bare `--outCsvCols=`.
    if out.get("columns"):
        cmd.append(f"--outCsvCols={out['columns']}")


#----- Assembles the IntaRNA command line for one (query, target); `acc` is the resolved --acc mode (N=none, C=constrained) from the Snakefile -----#
def build_command(cfg, acc, query, target, out_path, threads, max_suboptimal_hits, out_max_energy, derived_seed):
    cmd = [
        which_required("IntaRNA"),
        "-q", query,
        "-t", target,
        "--outMode=C",
        "--out", out_path,
        f"--threads={threads}",
        f"--acc={acc}",
    ]

    # Interaction model (applies to both variants)
    cmd += [
        f"--mode={opt(cfg, 'prediction_mode', 'H')}",
        f"--model={opt(cfg, 'model', INTARNA_DEFAULT_MODEL)}",
        f"--energyVRNA={opt(cfg, 'energy_set', INTARNA_DEFAULT_ENERGY_SET)}",
        f"--intLenMax={opt(cfg, 'max_interaction_length', 0)}",
        f"--intLoopMax={opt(cfg, 'max_loop_size', 10)}",
    ]

    # Seed (inherited from the top-level seed; null → --noSeed)
    _add_seed_flags(cmd, cfg.get("seed") or {}, derived_seed=derived_seed)

    # Helix (only honoured under --model=B, but harmless to pass always)
    _add_helix_flags(cmd, cfg.get("helix") or {})

    # Accessibility (only meaningful when --acc=C, but harmless to pass always)
    _add_accessibility_flags(cmd, cfg.get("accessibility") or {})

    _add_output_flags(cmd, cfg.get("output") or {}, max_suboptimal_hits=max_suboptimal_hits, out_max_energy=out_max_energy)

    return cmd


#----- Top-level driver: builds the command and runs IntaRNA once per (sample × accessibility variant) -----#
def run_intarna(query, target, out_path, acc, cfg, threads, max_suboptimal_hits, out_max_energy, derived_seed):
    ensure_parent(out_path)
    cmd = build_command(cfg, acc, query, target, out_path, threads, max_suboptimal_hits, out_max_energy, derived_seed)
    subprocess.run(cmd, check=True)


#----- Snakemake entry point: unpacks the injected `snakemake` object and calls the pure logic above -----#
def run_from_snakemake(snakemake):
    derived_seed = derive_seed_from_string(snakemake.params.seed)
    run_intarna(
        query=snakemake.input.query,
        target=snakemake.input.target,
        out_path=snakemake.output.csv,
        acc=snakemake.params.acc,
        cfg=dict(snakemake.params.intarna),
        threads=snakemake.threads,
        max_suboptimal_hits=snakemake.params.max_suboptimal_hits,
        out_max_energy=snakemake.params.out_max_energy,
        derived_seed=derived_seed,
    )

if "snakemake" in globals():
    run_from_snakemake(snakemake)
