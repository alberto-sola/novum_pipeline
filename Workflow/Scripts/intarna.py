import subprocess
from pathlib import Path


WO_ACCESSIBILITY = "wo_accessibility"


def _add_seed_flags(cmd, seed):
    if not seed:
        return
    cmd.append(f"--seedBP={seed.get('bp', 7)}")
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
    if seed.get("q_range") is not None:
        cmd.append(f"--seedQRange={seed['q_range']}")
    if seed.get("t_range") is not None:
        cmd.append(f"--seedTRange={seed['t_range']}")
    if seed.get("max_UP") is not None:
        cmd.append(f"--seedMaxUP={seed['max_UP']}")
    if seed.get("out_best_only"):
        cmd.append("--outBestSeedOnly")


def _add_accessibility_flags(cmd, acc):
    if not acc:
        return
    if acc.get("w") is not None:
        cmd.append(f"--accW={acc['w']}")
    if acc.get("l") is not None:
        cmd.append(f"--accL={acc['l']}")
    if acc.get("no_lp"):
        cmd.append("--accNoLP")
    if acc.get("no_gu_end"):
        cmd.append("--accNoGUend")


def _add_output_flags(cmd, out):
    if not out:
        return
    if out.get("max_E") is not None:
        cmd.append(f"--outMaxE={out['max_E']}")
    if out.get("delta_E") is not None:
        cmd.append(f"--outDeltaE={out['delta_E']}")
    cmd.append(f"--outNumber={out.get('number', 1)}")
    cmd.append(f"--outOverlap={out.get('overlap', 'B')}")
    if out.get("min_Pu") is not None:
        cmd.append(f"--outMinPu={out['min_Pu']}")
    if out.get("no_lp"):
        cmd.append("--outNoLP")
    if out.get("no_gu_end"):
        cmd.append("--outNoGUend")
    if out.get("csv_cols"):
        cmd.append(f"--outCsvCols={out['csv_cols']}")


def build_command(cfg, variant, query, target, out_path, threads):
    cmd = [
        "IntaRNA",
        "-q", query,
        "-t", target,
        "--outMode=C",
        "--out", out_path,
        f"--threads={threads}",
    ]

    if variant == WO_ACCESSIBILITY:
        # §10.1 RNAhybrid-emulation profile: no seed, no accessibility, exact mode.
        cmd += ["--noSeed", "--acc=N", "--intLoopMax=30", "--mode=M"]
    else:
        cmd += [
            "--acc=C",
            f"--mode={cfg.get('mode', 'H')}",
            f"--model={cfg.get('model', 'X')}",
            f"--intLenMax={cfg.get('int_len_max', 0)}",
            f"--intLoopMax={cfg.get('int_loop_max', 10)}",
        ]
        _add_seed_flags(cmd, cfg.get("seed") or {})
        _add_accessibility_flags(cmd, cfg.get("accessibility") or {})

    _add_output_flags(cmd, cfg.get("output") or {})
    cmd += list(cfg.get("extra_args") or [])

    return cmd


def run_intarna(query, target, out_path, variant, cfg, threads):
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    cmd = build_command(cfg, variant, query, target, out_path, threads)
    subprocess.run(cmd, check=True)


def run_from_snakemake(snakemake):
    run_intarna(
        query=snakemake.input.query,
        target=snakemake.input.target,
        out_path=snakemake.output.csv,
        variant=snakemake.wildcards.variant,
        cfg=dict(snakemake.params.intarna),
        threads=snakemake.threads,
    )

run_from_snakemake(snakemake)
