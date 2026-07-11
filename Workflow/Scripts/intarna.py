import subprocess
import sys

from _common import which_required, ensure_parent


INTARNA_MAX_OUTNUMBER = 1000

#----- Parses the shared top-level seed "x,y" into IntaRNA's seedBP+seedQRange shape; raises if the width is outside [2,20] -----#
def _derive_seed_from_string(seed_str):
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


#----- Appends --out* filtering flags; `max_suboptimal_hits` and `max_total_energy` come from the shared top-level config, the rest from intarna.output -----#
def _add_output_flags(cmd, out, max_suboptimal_hits=None, max_total_energy=None):
    if max_total_energy is not None:
        cmd.append(f"--outMaxE={max_total_energy}")
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
def build_command(cfg, acc, query, target, out_path, threads, max_suboptimal_hits, max_total_energy, derived_seed):
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
        f"--mode={cfg.get('prediction_mode', 'H')}",
        f"--model={cfg.get('model', 'S')}",
        f"--intLenMax={cfg.get('max_interaction_length', 0)}",
        f"--intLoopMax={cfg.get('max_loop_size', 10)}",
    ]

    # Seed (inherited from the top-level seed; null → --noSeed)
    _add_seed_flags(cmd, cfg.get("seed") or {}, derived_seed=derived_seed)

    # Accessibility (only meaningful when --acc=C, but harmless to pass always)
    _add_accessibility_flags(cmd, cfg.get("accessibility") or {})

    _add_output_flags(cmd, cfg.get("output") or {}, max_suboptimal_hits=max_suboptimal_hits, max_total_energy=max_total_energy)
    cmd += list(cfg.get("extra_args") or [])

    return cmd


#----- Top-level driver: builds the command and runs IntaRNA once per (sample × accessibility variant) -----#
def run_intarna(query, target, out_path, acc, cfg, threads, max_suboptimal_hits, max_total_energy, derived_seed):
    ensure_parent(out_path)
    cmd = build_command(cfg, acc, query, target, out_path, threads, max_suboptimal_hits, max_total_energy, derived_seed)
    subprocess.run(cmd, check=True)


#----- Snakemake entry point: unpacks the injected `snakemake` object and calls the pure logic above -----#
def run_from_snakemake(snakemake):
    derived_seed = _derive_seed_from_string(snakemake.params.seed)
    run_intarna(
        query=snakemake.input.query,
        target=snakemake.input.target,
        out_path=snakemake.output.csv,
        acc=snakemake.params.acc,
        cfg=dict(snakemake.params.intarna),
        threads=snakemake.threads,
        max_suboptimal_hits=snakemake.params.max_suboptimal_hits,
        max_total_energy=snakemake.params.max_total_energy,
        derived_seed=derived_seed,
    )

if "snakemake" in globals():
    run_from_snakemake(snakemake)
