import json
import math
import shutil
import subprocess
from pathlib import Path


#----- Launches a message error if RNAcalibrate is not installed -----#
def ensure_dependency():
    executable = shutil.which("RNAcalibrate") or shutil.which("rnacalibrate")
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


#----- Calculates mean of length and standard deviation of the annotations in the target FASTA file -----#
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


#----- Rounds mean and std values and return them in the <x>,<y> format -----#
def build_length_arg(stats):
    mean_length = int(round(stats["mean"]))
    std_length = int(round(stats["std"]))

    return {
        "mean": mean_length,
        "std": std_length,
        "value": f"{mean_length},{std_length}",
    }


#----- Builds the RNAcalibrate bash command -----#
def build_command(executable, query, target, k, max_target_length, length_arg, randomize_targets=False, u=None, v=None, seed=None):
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
    if u is not None:
        command.extend(["-u", str(u)])
    if v is not None:
        command.extend(["-v", str(v)])
    if seed is not None:
        command.extend(["-f", str(seed)])
    if randomize_targets:
        command.append("-s")

    return command


#----- Parses the calibration values from the RNAcalibrate bash output and returns a dictionary of per-query values -----#
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

        query_name = fields[0]
        sample_size = int(fields[1])
        xi = float(fields[2])
        theta = float(fields[3])

        # Launches an error message is NaN is output
        if math.isnan(xi) or math.isnan(theta):
            raise RuntimeError(f"RNAcalibrate produced NaN parameters: {raw_line}")

        per_query.append(
            {
                "query": query_name,
                "sample_size": sample_size,
                "xi": xi,
                "theta": theta,
            }
        )
    # Launches an error message if RNAcalibrate output is empty
    if not per_query:
        raise RuntimeError("RNAcalibrate did not produce any calibration rows.")

    return {"per_query": per_query}


#----- Runs RNAcalibrate and craft a JSON file -----#
def run_rnacalibrate(query, target, output_file, k, max_target_length, randomize_targets=False, u=None, v=None, seed=None):
    executable = ensure_dependency()
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
        u=u,
        v=v,
        seed=seed,
    )

    # Here RNAcalibrate is run and the output is stored
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    parsed_output = parse_rnacalibrate_output(completed.stdout)

    # Craft a JSON file containing metadata of RNAcalibrate command + per-query distribution values
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(
        {
            "command": command,
            "target_length_stats": stats,
            "target_length_argument": length_arg,
            "calibration": parsed_output,
        }, indent=2) + "\n"
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