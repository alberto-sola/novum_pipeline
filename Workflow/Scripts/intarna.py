import subprocess
from pathlib import Path


WO_ACCESSIBILITY = "wo_accessibility"


def _derive_seed_from_string(seed_str):
    """Parse top-level seed "x,y" and return derivation dict for --seedBP/--seedQRange/--seedMaxUP=0."""
    if seed_str is None:
        return None
    try:
        parts = str(seed_str).split(",")
        if len(parts) != 2:
            raise ValueError
        x, y = int(parts[0].strip()), int(parts[1].strip())
    except (ValueError, IndexError):
        raise ValueError(f"top-level seed must be 'x,y' integers, got: {seed_str!r}")
    bp = y - x + 1
    if not (2 <= bp <= 20):
        raise ValueError(
            f"top-level seed '{seed_str}' gives seedBP={bp}, which is outside IntaRNA's [2,20] range"
        )
    return {"bp": bp, "q_range": f"{x}-{y}"}


def _add_seed_flags(cmd, seed, derived_seed=None):
    if seed.get("disable"):
        cmd.append("--noSeed")
        return

    # --seedBP
    if seed.get("bp") is not None:
        cmd.append(f"--seedBP={seed['bp']}")
    elif derived_seed is not None:
        cmd.append(f"--seedBP={derived_seed['bp']}")

    # --seedQRange
    if seed.get("q_range") is not None:
        cmd.append(f"--seedQRange={seed['q_range']}")
    elif derived_seed is not None:
        cmd.append(f"--seedQRange={derived_seed['q_range']}")

    # --seedMaxUP
    if seed.get("max_UP") is not None:
        cmd.append(f"--seedMaxUP={seed['max_UP']}")
    elif derived_seed is not None:
        cmd.append("--seedMaxUP=0")

    # Remaining sub-keys
    if seed.get("max_E") is not None:
        cmd.append(f"--seedMaxE={seed['max_E']}")
    if seed.get("max_E_hybrid") is not None:
        cmd.append(f"--seedMaxEhybrid={seed['max_E_hybrid']}")
    if seed.get("min_Pu") is not None:
        cmd.append(f"--seedMinPu={seed['min_Pu']}")
    if seed.get("no_GU"):
        cmd.append("--seedNoGU")
    if seed.get("no_GU_end"):
        cmd.append("--seedNoGUend")
    if seed.get("t_range") is not None:
        cmd.append(f"--seedTRange={seed['t_range']}")
    if seed.get("out_best_only"):
        cmd.append("--outBestSeedOnly")


def _add_accessibility_flags(cmd, acc):
    if not acc:
        return
    if acc.get("window") is not None:
        cmd.append(f"--accW={acc['window']}")
    if acc.get("max_bp_span") is not None:
        cmd.append(f"--accL={acc['max_bp_span']}")
    if acc.get("no_lp"):
        cmd.append("--accNoLP")
    if acc.get("no_gu_end"):
        cmd.append("--accNoGUend")


def _add_output_flags(cmd, out, hits=None, max_energy=None):
    if max_energy is not None:
        cmd.append(f"--outMaxE={max_energy}")
    if out.get("delta_E") is not None:
        cmd.append(f"--outDeltaE={out['delta_E']}")
    if hits is not None:
        cmd.append(f"--outNumber={hits}")
    cmd.append(f"--outOverlap={out.get('overlap', 'B')}")
    if out.get("min_Pu") is not None:
        cmd.append(f"--outMinPu={out['min_Pu']}")
    if out.get("no_lp"):
        cmd.append("--outNoLP")
    if out.get("no_gu_end"):
        cmd.append("--outNoGUend")
    if out.get("csv_cols"):
        cmd.append(f"--outCsvCols={out['csv_cols']}")


def build_command(cfg, variant, query, target, out_path, threads, hits, max_energy, derived_seed):
    cmd = [
        "IntaRNA",
        "-q", query,
        "-t", target,
        "--outMode=C",
        "--out", out_path,
        f"--threads={threads}",
    ]

    # Accessibility
    if variant == WO_ACCESSIBILITY:
        cmd.append("--acc=N")
    else:
        cmd.append("--acc=C")

    # Interaction model (applies to both variants)
    cmd += [
        f"--mode={cfg.get('mode', 'H')}",
        f"--model={cfg.get('model', 'S')}",
        f"--intLenMax={cfg.get('int_len_max', 0)}",
        f"--intLoopMax={cfg.get('int_loop_max', 10)}",
    ]

    # Seed (applies to both variants; disable: true → --noSeed)
    _add_seed_flags(cmd, cfg.get("seed") or {}, derived_seed=derived_seed)

    # Accessibility params (only meaningful when --acc=C, but harmless to pass always)
    _add_accessibility_flags(cmd, cfg.get("accessibility_params") or {})

    _add_output_flags(cmd, cfg.get("output") or {}, hits=hits, max_energy=max_energy)
    cmd += list(cfg.get("extra_args") or [])

    return cmd


def run_intarna(query, target, out_path, variant, cfg, threads, hits, max_energy, derived_seed):
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    cmd = build_command(cfg, variant, query, target, out_path, threads, hits, max_energy, derived_seed)
    subprocess.run(cmd, check=True)


def run_from_snakemake(snakemake):
    derived_seed = _derive_seed_from_string(snakemake.params.seed)
    run_intarna(
        query=snakemake.input.query,
        target=snakemake.input.target,
        out_path=snakemake.output.csv,
        variant=snakemake.wildcards.variant,
        cfg=dict(snakemake.params.intarna),
        threads=snakemake.threads,
        hits=snakemake.params.hits,
        max_energy=snakemake.params.max_energy,
        derived_seed=derived_seed,
    )

run_from_snakemake(snakemake)
