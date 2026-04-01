import json
import math
import shutil
import subprocess
from pathlib import Path


def ensure_dependency():
    executable = shutil.which("RNAcalibrate")
    if executable is None:
        raise RuntimeError("Required executable not found in PATH: RNAcalibrate")
    return executable


#----- Counts nucleotides after each FASTA ">" annotation -----#
def iter_fasta_lengths(path):
    current_length = 0

    with open(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_length:
                    yield current_length
                current_length = 0
                continue

            current_length += len(line)

    if current_length:
        yield current_length


#----- Calculates the mean length of the annotations and the standard deviation... -----#
def compute_target_length_stats(target_file):
    lengths = list(iter_fasta_lengths(target_file))
    if not lengths:
        raise RuntimeError(f"No FASTA records found in target file: {target_file}")

    # Total sum of nucleotides divided by how many lengths (annotations) there are
    mean_length = sum(lengths) / len(lengths)
    if len(lengths) == 1:
        std_length = 0.0
    else:
        variance = sum((length - mean_length) ** 2 for length in lengths) / len(lengths)
        std_length = math.sqrt(variance)

    return {
        "count": len(lengths),
        "mean": mean_length,
        "std": std_length,
    }


def build_command(executable, query, target, k, max_target_length, stats, randomize_targets=False, u=None, v=None, seed=None):
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
    ]

    # RNAcalibrate appears to be sensitive to argument order here:
    # with some builds, `-s` must appear before `-l` to avoid NaN fits.
    if randomize_targets:
        command.append("-s")

    command.extend(
        [
            "-l",
            f"{stats['mean']:.6f},{stats['std']:.6f}",
        ]
    )

    if u is not None:
        command.extend(["-u", str(u)])
    if v is not None:
        command.extend(["-v", str(v)])
    if seed is not None:
        command.extend(["-f", str(seed)])

    return command


#----- -----#
def parse_rnacalibrate_output(stdout):
    per_query = []

    for raw_line in stdout.splitlines():
        line = raw_line.strip()
        if not line:
            continue

        fields = line.split()
        if len(fields) < 4:
            raise RuntimeError(f"Unexpected RNAcalibrate output line: {raw_line}")

        query_name = fields[0]
        sample_size = int(fields[1])
        xi = float(fields[2])
        theta = float(fields[3])
        per_query.append(
            {
                "query": query_name,
                "sample_size": sample_size,
                "xi": xi,
                "theta": theta,
            }
        )

    if not per_query:
        raise RuntimeError("RNAcalibrate did not produce any calibration rows.")

    mean_xi = sum(item["xi"] for item in per_query) / len(per_query)
    mean_theta = sum(item["theta"] for item in per_query) / len(per_query)

    return {
        "per_query": per_query,
        "mean_xi": mean_xi,
        "mean_theta": mean_theta,
        "distribution": f"{mean_xi:.6f},{mean_theta:.6f}",
    }


def run_rnacalibrate(query, target, output_file, k, max_target_length, randomize_targets=False, u=None, v=None, seed=None):
    executable = ensure_dependency()
    stats = compute_target_length_stats(target)
    command = build_command(
        executable=executable,
        query=query,
        target=target,
        k=k,
        max_target_length=max_target_length,
        stats=stats,
        randomize_targets=randomize_targets,
        u=u,
        v=v,
        seed=seed,
    )

    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    parsed_output = parse_rnacalibrate_output(completed.stdout)

    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(
            {
                "command": command,
                "target_length_stats": stats,
                "calibration": parsed_output,
                "raw_stdout": completed.stdout,
            },
            indent=2,
        )
        + "\n"
    )


def run_from_snakemake(snakemake):
    run_rnacalibrate(
        query=snakemake.input.query,
        target=snakemake.input.target,
        output_file=snakemake.output.calibration,
        k=snakemake.params.k,
        max_target_length=snakemake.params.max_target_length,
        randomize_targets=snakemake.params.randomize_targets,
        u=snakemake.params.u,
        v=snakemake.params.v,
        seed=snakemake.params.seed,
    )

run_from_snakemake(snakemake)
