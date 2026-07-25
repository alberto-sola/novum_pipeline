import json
import math
import subprocess

from _common import which_required, iter_fasta_records, ensure_parent


#----- Mean and population stdev of the target FASTA's record lengths (the input to RNAcalibrate's `-l`) -----#
def compute_target_length_stats(target_file):
    lengths = [len(seq) for _header, seq in iter_fasta_records(target_file) if seq]
    if not lengths:
        raise RuntimeError(f"No FASTA records found in target file: {target_file}")

    mean_length = sum(lengths) / len(lengths)
    if len(lengths) == 1:
        std_length = 0.0
    else:
        variance = sum((length - mean_length) ** 2 for length in lengths) / len(lengths)
        std_length = math.sqrt(variance)

    return {"count": len(lengths), "mean": mean_length, "std": std_length}


#----- Rounds the (mean, std) pair to ints and formats them as RNAcalibrate's "<mean>,<std>" -l argument -----#
def build_length_arg(stats):
    mean_length = int(round(stats["mean"]))
    std_length = int(round(stats["std"]))
    return {
        "mean": mean_length,
        "std": std_length,
        "value": f"{mean_length},{std_length}",
    }


#----- Assembles the RNAcalibrate command line; optional flags (max_internal_loop, max_bulge_loop, seed, randomize) are appended only if set -----#
def build_command(executable, query, target, k, max_target_length, length_arg, randomize_targets=False, max_internal_loop=None, max_bulge_loop=None, seed=None):
    command = [
        executable,
        "-k",
        str(k),
        "-q",
        str(query),
        "-t",
        str(target),
        "-m",
        str(max_target_length),
        "-l",
        length_arg["value"],
    ]
    if max_internal_loop is not None:
        command.extend(["-u", str(max_internal_loop)])
    if max_bulge_loop is not None:
        command.extend(["-v", str(max_bulge_loop)])
    if seed is not None:
        command.extend(["-f", str(seed)])
    if randomize_targets:
        command.append("-s")

    return command


#----- Parses the four-column RNAcalibrate stdout into one xi/theta record per query miRNA -----#
def parse_rnacalibrate_output(stdout):
    per_query = []

    for raw_line in stdout.splitlines():
        line = raw_line.strip()
        if not line:
            continue

        fields = line.split()
        if len(fields) != 4:
            raise RuntimeError(
                "Unexpected RNAcalibrate output: expected 4 columns "
                f"(query, sample_size, xi, theta), got {len(fields)} in line: {raw_line!r}"
            )

        xi = float(fields[2])
        theta = float(fields[3])
        # NaN sneaks in when RNAcalibrate's sample is degenerate; refuse to
        # let it poison the downstream RNAhybrid -d argument.
        if math.isnan(xi) or math.isnan(theta):
            raise RuntimeError(f"RNAcalibrate produced NaN parameters: {raw_line}")

        per_query.append({
            "query": fields[0],
            "sample_size": int(fields[1]),
            "xi": xi,
            "theta": theta,
        })

    if not per_query:
        raise RuntimeError("RNAcalibrate did not produce any calibration rows.")

    return {"per_query": per_query}


#----- Top-level driver: stat the target FASTA, run RNAcalibrate, persist command + result as JSON -----#
def run_rnacalibrate(query, target, output_file, k, max_target_length, randomize_targets=False, max_internal_loop=None, max_bulge_loop=None, seed=None):
    executable = which_required("RNAcalibrate", "rnacalibrate")
    stats = compute_target_length_stats(target)
    length_arg = build_length_arg(stats)
    command = build_command(
        executable=executable,
        query=query,
        target=target,
        k=k,
        max_target_length=max_target_length,
        length_arg=length_arg,
        randomize_targets=randomize_targets,
        max_internal_loop=max_internal_loop,
        max_bulge_loop=max_bulge_loop,
        seed=seed,
    )

    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    parsed_output = parse_rnacalibrate_output(completed.stdout)

    ensure_parent(output_file).write_text(json.dumps(
        {
            "command": command,
            "target_length_stats": stats,
            "target_length_argument": length_arg,
            "calibration": parsed_output,
        }, indent=2) + "\n"
    )


#----- Snakemake entry point: unpacks the injected `snakemake` object and calls the pure logic above -----#
def run_from_snakemake(snakemake):
    run_rnacalibrate(
        query=snakemake.input.query,
        target=snakemake.input.target,
        output_file=snakemake.output.calibration,
        k=snakemake.params.k,
        max_target_length=snakemake.params.max_target_length,
        randomize_targets=snakemake.params.randomize_targets,
        max_internal_loop=snakemake.params.max_internal_loop,
        max_bulge_loop=snakemake.params.max_bulge_loop,
        seed=snakemake.params.seed,
    )

if "snakemake" in globals():
    run_from_snakemake(snakemake)