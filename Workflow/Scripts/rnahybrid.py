from pathlib import Path
import shutil
import subprocess


def write_fasta_chunks(target_file, chunk_prefix, max_lines=800):
    """Split a FASTA file into chunk files while preserving full records."""
    chunk_paths = []
    current_chunk_index = 0
    current_chunk_line_total = 0
    current_record = []
    current_chunk_records = []

    def flush_chunk():
        nonlocal current_chunk_index, current_chunk_line_total, current_chunk_records
        if not current_chunk_records:
            return

        first_letter = chr(ord("a") + ((current_chunk_index // 26) % 26))
        second_letter = chr(ord("a") + (current_chunk_index % 26))
        chunk_path = Path(f"{chunk_prefix}{first_letter}{second_letter}")
        chunk_path.write_text("".join(current_chunk_records))
        chunk_paths.append(chunk_path)
        current_chunk_index += 1
        current_chunk_line_total = 0
        current_chunk_records = []

    def add_record(record_lines):
        nonlocal current_chunk_line_total, current_chunk_records
        record_line_total = len(record_lines)
        if current_chunk_records and current_chunk_line_total + record_line_total > max_lines:
            flush_chunk()
        current_chunk_records.extend(record_lines)
        current_chunk_line_total += record_line_total

    with open(target_file) as handle:
        for line in handle:
            if line.startswith(">"):
                if current_record:
                    add_record(current_record)
                current_record = [line]
            else:
                current_record.append(line)

    if current_record:
        add_record(current_record)
    flush_chunk()

    if not chunk_paths:
        raise RuntimeError(f"No FASTA records found in target file: {target_file}")

    return chunk_paths


def ensure_dependencies():
    parallel_executable = shutil.which("parallel")
    if parallel_executable is None:
        raise RuntimeError("Required executable not found in PATH: parallel")

    rnahybrid_executable = shutil.which("RNAhybrid") or shutil.which("rnahybrid")
    if rnahybrid_executable is None:
        raise RuntimeError("Required executable not found in PATH: RNAhybrid")

    return {
        "parallel": parallel_executable,
        "rnahybrid": rnahybrid_executable,
    }


def cleanup_temp_files(chunk_prefix, output_pattern):
    for path in Path(".").glob(output_pattern):
        path.unlink()
    for path in Path(".").glob(f"{chunk_prefix}*"):
        path.unlink()


def build_optional_args(hits=None, u=None, v=None, energy=None, pvalue=None, seed=None):
    optional_args = []

    if hits is not None:
        optional_args.extend(["-b", str(hits)])
    if u is not None:
        optional_args.extend(["-u", str(u)])
    if v is not None:
        optional_args.extend(["-v", str(v)])
    if energy is not None:
        optional_args.extend(["-e", str(energy)])
    if pvalue is not None:
        optional_args.extend(["-p", str(pvalue)])
    if seed is not None:
        optional_args.extend(["-f", str(seed)])

    return optional_args


def build_parallel_command(query, chunk_paths, species, optional_args, threads, parallel_executable, rnahybrid_executable):
    job_template = " ".join(
        [
            rnahybrid_executable,
            "-q",
            "{1}",
            "-t",
            "{2}",
            "-s",
            "{3}",
            "-c",
            "-m",
            "500000",
            *optional_args,
            ">",
            "output_{2/.}.tsv",
        ]
    )

    return [
        parallel_executable,
        f"-j{max(1, threads)}",
        "--load=100%",
        job_template,
        ":::",
        query,
        ":::",
        *[str(chunk_path) for chunk_path in chunk_paths],
        ":::",
        species,
    ]


def merge_output_files(output_pattern, merged_output_path):
    output_files = sorted(Path(".").glob(output_pattern))
    if not output_files:
        raise RuntimeError("No output files were produced by RNAhybrid.")

    with open(merged_output_path, "w") as merged_output:
        for output_file in output_files:
            merged_output.write(output_file.read_text())


def run_rnahybrid(query, target, species, output_file, threads=1, hits=None, u=None, v=None, energy=None, pvalue=None, seed=None):
    chunk_prefix = "chunk_"
    output_pattern = "output_chunk_*.tsv"

    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    executables = ensure_dependencies()
    cleanup_temp_files(chunk_prefix, output_pattern)

    try:
        chunk_paths = write_fasta_chunks(target, chunk_prefix)
        optional_args = build_optional_args(
            hits=hits,
            u=u,
            v=v,
            energy=energy,
            pvalue=pvalue,
            seed=seed,
        )
        command = build_parallel_command(
            query=query,
            chunk_paths=chunk_paths,
            species=species,
            optional_args=optional_args,
            threads=threads,
            parallel_executable=executables["parallel"],
            rnahybrid_executable=executables["rnahybrid"],
        )

        subprocess.run(command, check=True)
        merge_output_files(output_pattern, output_path)
    finally:
        cleanup_temp_files(chunk_prefix, output_pattern)


def run_from_snakemake(snakemake):
    run_rnahybrid(
        query=snakemake.input.query,
        target=snakemake.input.target,
        output_file=snakemake.output.compact,
        threads=snakemake.threads,
        species=snakemake.params.species,
        hits=snakemake.params.hits,
        u=snakemake.params.u,
        v=snakemake.params.v,
        energy=snakemake.params.energy,
        pvalue=snakemake.params.pvalue,
        seed=snakemake.params.seed,
    )


run_from_snakemake(snakemake)
