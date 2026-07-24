import subprocess
import sys

from _common import which_required, ensure_parent
from _intarna_config import (
    INTARNA_DEFAULT_ENERGY_SET,
    INTARNA_DEFAULT_MODEL,
    INTARNA_MAX_OUTNUMBER,
    derive_seed_from_string,
    opt,
)


#----- Seed enforcement is inherited from the top-level `seed`: absent → --noSeed. When present, --seedBP/--seedQRange are derived; the remaining sub-keys are independent overrides -----#
def _add_seed_flags(cmd, seed, derived_seed=None):
    if derived_seed is None:
        cmd.append("--noSeed")
        return

    cmd.append(f"--seedBP={derived_seed['length']}")
    cmd.append(f"--seedQRange={derived_seed['query_range']}")

    # --seedMaxUP: explicit override wins, else RNAhybrid-parity 0
    if seed.get("max_unpaired_bases") is not None:
        cmd.append(f"--seedMaxUP={seed['max_unpaired_bases']}")
    else:
        cmd.append("--seedMaxUP=0")

    if seed.get("max_energy") is not None:
        cmd.append(f"--seedMaxE={seed['max_energy']}")
    if seed.get("max_hybrid_energy") is not None:
        cmd.append(f"--seedMaxEhybrid={seed['max_hybrid_energy']}")
    if seed.get("min_unpaired_probability") is not None:
        cmd.append(f"--seedMinPu={seed['min_unpaired_probability']}")
    if seed.get("forbid_gu"):
        cmd.append("--seedNoGU")
    if seed.get("forbid_gu_at_ends"):
        cmd.append("--seedNoGUend")
    if seed.get("target_range") is not None:
        cmd.append(f"--seedTRange={seed['target_range']}")
    if seed.get("report_best_only"):
        cmd.append("--outBestSeedOnly")


#----- Appends --accW/--accL/--accNoLP/--accNoGUend from the accessibility block, skipping any null/false entries -----#
def _add_accessibility_flags(cmd, acc):
    if not acc:
        return
    if acc.get("window") is not None:
        cmd.append(f"--accW={acc['window']}")
    if acc.get("max_bp_span") is not None:
        cmd.append(f"--accL={acc['max_bp_span']}")
    if acc.get("forbid_lonely_pairs"):
        cmd.append("--accNoLP")
    if acc.get("forbid_gu_at_ends"):
        cmd.append("--accNoGUend")


#----- Appends --helix* flags from the helix block, skipping null/false entries. Emitted unconditionally: IntaRNA ignores them under any model but B (verified), and the mismatch is reported by validate_intarna_config rather than by withholding flags -----#
def _add_helix_flags(cmd, helix):
    if not helix:
        return
    if helix.get("min_bp") is not None:
        cmd.append(f"--helixMinBP={helix['min_bp']}")
    if helix.get("max_bp") is not None:
        cmd.append(f"--helixMaxBP={helix['max_bp']}")
    if helix.get("max_internal_loop") is not None:
        cmd.append(f"--helixMaxIL={helix['max_internal_loop']}")
    if helix.get("min_unpaired_probability") is not None:
        cmd.append(f"--helixMinPu={helix['min_unpaired_probability']}")
    if helix.get("max_energy") is not None:
        cmd.append(f"--helixMaxE={helix['max_energy']}")
    if helix.get("full_energy"):
        cmd.append("--helixFullE")


#----- Appends --out* filtering flags. `out_max_energy` is IntaRNA's TOTAL-energy bound (--outMaxE), supplied by the Snakefile only where it equals the hybridization gate (acc=N); `max_suboptimal_hits` is the shared knob, the rest come from intarna.output -----#
def _add_output_flags(cmd, out, max_suboptimal_hits=None, out_max_energy=None):
    if out_max_energy is not None:
        cmd.append(f"--outMaxE={out_max_energy}")
    if out.get("max_delta_energy") is not None:
        cmd.append(f"--outDeltaE={out['max_delta_energy']}")
    if max_suboptimal_hits is None:
        n_hits = INTARNA_MAX_OUTNUMBER
    elif max_suboptimal_hits > INTARNA_MAX_OUTNUMBER:
        print(
            f"intarna: max_suboptimal_hits={max_suboptimal_hits} exceeds IntaRNA's "
            f"--outNumber ceiling; clamping to {INTARNA_MAX_OUTNUMBER}.",
            file=sys.stderr,
        )
        n_hits = INTARNA_MAX_OUTNUMBER
    else:
        n_hits = max_suboptimal_hits
    cmd.append(f"--outNumber={n_hits}")
    cmd.append(f"--outOverlap={out.get('overlap', 'B')}")
    if out.get("min_unpaired_probability") is not None:
        cmd.append(f"--outMinPu={out['min_unpaired_probability']}")
    if out.get("forbid_lonely_pairs"):
        cmd.append("--outNoLP")
    if out.get("forbid_gu_at_ends"):
        cmd.append("--outNoGUend")
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
