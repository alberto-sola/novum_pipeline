from pathlib import Path
import json
import shutil
import subprocess
import sys
import tempfile

from _common import which_required, ensure_parent


#----- Splits the target FASTA into chunk files of bounded line count, never breaking a record across chunks -----#
def write_fasta_chunks(target_file, chunk_dir, chunk_prefix="chunk_", max_lines=800):
    chunk_dir = Path(chunk_dir)
    chunk_dir.mkdir(parents=True, exist_ok=True)

    chunk_paths = []
    current_chunk_index = 0
    current_chunk_line_total = 0
    current_record = []
    current_chunk_records = []

    def flush_chunk():
        nonlocal current_chunk_index, current_chunk_line_total, current_chunk_records
        if not current_chunk_records:
            return

        # Zero-padded numeric suffix so we do not silently wrap around on large targets.
        chunk_path = chunk_dir / f"{chunk_prefix}{current_chunk_index:06d}"
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


#----- Explodes a multi-record query FASTA into one file per miRNA, keyed by the first header token -----#
def split_query_per_miRNA(query_file, query_dir):
    query_dir = Path(query_dir)
    query_dir.mkdir(parents=True, exist_ok=True)

    paths_by_name = {}
    current_name = None
    current_record = []

    def flush():
        nonlocal current_name, current_record
        if current_name is None or not current_record:
            return
        if current_name in paths_by_name:
            raise RuntimeError(f"Duplicate query ID in FASTA: {current_name}")
        path = query_dir / f"{current_name}.fa"
        path.write_text("".join(current_record))
        paths_by_name[current_name] = path
        current_name = None
        current_record = []

    with open(query_file) as handle:
        for line in handle:
            if line.startswith(">"):
                flush()
                # Match RNAcalibrate's name parsing: first whitespace-delimited token after '>'.
                current_name = line[1:].strip().split()[0]
                current_record = [line]
            else:
                current_record.append(line)

    flush()

    if not paths_by_name:
        raise RuntimeError(f"No FASTA records found in query file: {query_file}")

    return paths_by_name


#----- Resolves the RNAhybrid + GNU Parallel binaries this rule depends on -----#
def ensure_dependencies():
    return {
        "parallel": which_required("parallel"),
        "rnahybrid": which_required("RNAhybrid", "rnahybrid", label="RNAhybrid"),
    }


#----- Translates the optional RNAhybrid params into CLI flags, omitting any that are None -----#
def build_optional_args(max_suboptimal_hits=None, max_internal_loop=None, max_bulge_loop=None, max_total_energy=None, pvalue_threshold=None, seed=None, distribution=None):
    # (flag, value) in emission order; each pair is appended only when set.
    flag_values = [
        ("-b", max_suboptimal_hits),
        ("-u", max_internal_loop),
        ("-v", max_bulge_loop),
        ("-e", max_total_energy),
        ("-p", pvalue_threshold),
        ("-f", seed),
        ("-d", distribution),
    ]

    optional_args = []
    for flag, value in flag_values:
        if value is not None:
            optional_args.extend([flag, str(value)])

    return optional_args


#----- Reads per-miRNA xi/theta from the calibration JSON; returns None on the uncalibrated path -----#
def load_per_query_distributions(distribution_file):
    if distribution_file is None:
        return None

    payload = json.loads(Path(distribution_file).read_text())
    per_query = payload.get("calibration", {}).get("per_query")
    if not per_query:
        raise RuntimeError(f"No per-query distributions found in calibration file: {distribution_file}")

    return {entry["query"]: f"{entry['xi']:.6f},{entry['theta']:.6f}" for entry in per_query}


#----- Resolves the species/distribution conflict: an explicit distribution wins, species is dropped with a warning -----#
def validate_rnahybrid_args(species=None, distribution=None):
    if distribution is None and species is None:
        raise ValueError("Either 'species' or 'distribution' must be provided for RNAhybrid.")

    if distribution is not None:
        if species is not None:
            print(
                f"rnahybrid: both 'species' ({species!r}) and 'distribution' set; "
                "RNAhybrid will use the static distribution and ignore species.",
                file=sys.stderr,
            )
        return None

    return species


#----- Materializes the (miRNA × target chunk) cross-product as a TSV that GNU Parallel reads line-by-line -----#
def build_job_spec_tsv(query_paths_by_name, dist_map, chunk_paths, tsv_path):
    missing = sorted(name for name in query_paths_by_name if name not in dist_map)
    if missing:
        raise RuntimeError(
            "miRNA(s) in query FASTA have no calibration entry: " + ", ".join(missing)
        )

    lines = []
    for name, query_path in query_paths_by_name.items():
        distribution = dist_map[name]
        for chunk_path in chunk_paths:
            lines.append(
                f"{Path(query_path).as_posix()}\t{distribution}\t{Path(chunk_path).as_posix()}"
            )

    tsv_path = ensure_parent(tsv_path)
    tsv_path.write_text("\n".join(lines) + "\n")
    return tsv_path


#----- GNU Parallel command for the uncalibrated path: one shared query fanned out across every target chunk -----#
def build_parallel_command_broadcast(query, chunk_paths, output_dir, species, optional_args, threads, max_target_length, parallel_executable, rnahybrid_executable):
    species_args = ["-s", "{3}"] if species is not None else []
    # Keep the '{2/}' parallel replacement literal — build the template as a plain string.
    output_template = f"{Path(output_dir).as_posix()}/output_{{2/}}.tsv"

    job_template = " ".join([
        rnahybrid_executable,
        "-q", "{1}",
        "-t", "{2}",
        *species_args,
        "-c",
        "-m", str(max_target_length),
        *optional_args,
        ">", output_template,
    ])

    command = [
        parallel_executable, f"-j{max(1, threads)}",
        job_template,
        ":::", str(query),
        ":::", *[str(chunk_path) for chunk_path in chunk_paths],
    ]
    if species is not None:
        command.extend([":::", species])
    return command, "output_chunk_*.tsv"


#----- GNU Parallel command for the calibrated path: each TSV row is (query_path, distribution, chunk_path) -----#
def build_parallel_command_calibrated(spec_path, output_dir, optional_args, threads, max_target_length, parallel_executable, rnahybrid_executable):
    # {1/.} = query basename without extension; {3/.} = chunk basename without extension.
    output_template = f"{Path(output_dir).as_posix()}/output_{{1/.}}__{{3/.}}.tsv"

    job_template = " ".join([
        rnahybrid_executable,
        "-q", "{1}",
        "-d", "{2}",
        "-t", "{3}",
        "-c",
        "-m", str(max_target_length),
        *optional_args,
        ">", output_template,
    ])

    command = [
        parallel_executable, f"-j{max(1, threads)}",
        "--colsep", "\\t",
        job_template,
        "::::", str(spec_path),
    ]
    return command, "output_*__chunk_*.tsv"


#----- Concatenates the per-chunk RNAhybrid outputs into one TSV, streaming through the kernel to bound memory -----#
def merge_output_files(output_dir, output_pattern, merged_output_path):
    output_files = sorted(Path(output_dir).glob(output_pattern))
    if not output_files:
        raise RuntimeError("No output files were produced by RNAhybrid.")

    with open(merged_output_path, "wb") as merged_output:
        for output_file in output_files:
            with open(output_file, "rb") as src:
                shutil.copyfileobj(src, merged_output, length=1024 * 1024)


#----- Top-level driver: chunk the target, dispatch RNAhybrid via GNU Parallel, merge per-chunk outputs back together -----#
def run_rnahybrid(query, target, species, output_file, max_target_length, threads=1, max_suboptimal_hits=None, max_internal_loop=None, max_bulge_loop=None, max_total_energy=None, pvalue_threshold=None, seed=None, distribution=None, distribution_file=None):
    output_path = ensure_parent(output_file)

    executables = ensure_dependencies()

    # Per-sample tmp dir next to the final output so parallel samples never collide.
    tmp_dir = Path(tempfile.mkdtemp(prefix=".rnahybrid_", dir=output_path.parent))
    chunk_dir = tmp_dir / "chunks"
    split_output_dir = tmp_dir / "outputs"
    split_output_dir.mkdir(parents=True, exist_ok=True)

    try:
        chunk_paths = write_fasta_chunks(target, chunk_dir)
        dist_map = load_per_query_distributions(distribution_file)

        if dist_map is not None:
            # Calibrated branch: one RNAhybrid invocation per (miRNA, chunk), each with its own -d.
            # `species` and `distribution` are intentionally ignored here — calibration provides
            # the per-miRNA xi/theta that supersede both.
            query_paths_by_name = split_query_per_miRNA(query, tmp_dir / "queries")
            spec_path = build_job_spec_tsv(
                query_paths_by_name=query_paths_by_name,
                dist_map=dist_map,
                chunk_paths=chunk_paths,
                tsv_path=tmp_dir / "job_spec.tsv",
            )
            optional_args = build_optional_args(
                max_suboptimal_hits=max_suboptimal_hits, max_internal_loop=max_internal_loop, max_bulge_loop=max_bulge_loop, max_total_energy=max_total_energy, pvalue_threshold=pvalue_threshold, seed=seed,
            )
            command, output_pattern = build_parallel_command_calibrated(
                spec_path=spec_path,
                output_dir=split_output_dir,
                optional_args=optional_args,
                threads=threads,
                max_target_length=max_target_length,
                parallel_executable=executables["parallel"],
                rnahybrid_executable=executables["rnahybrid"],
            )
        else:
            # Broadcast/uncalibrated branch: one shared query across all chunks.
            species = validate_rnahybrid_args(species=species, distribution=distribution)
            optional_args = build_optional_args(
                max_suboptimal_hits=max_suboptimal_hits, max_internal_loop=max_internal_loop, max_bulge_loop=max_bulge_loop, max_total_energy=max_total_energy, pvalue_threshold=pvalue_threshold, seed=seed,
                distribution=distribution,
            )
            command, output_pattern = build_parallel_command_broadcast(
                query=query,
                chunk_paths=chunk_paths,
                output_dir=split_output_dir,
                species=species,
                optional_args=optional_args,
                threads=threads,
                max_target_length=max_target_length,
                parallel_executable=executables["parallel"],
                rnahybrid_executable=executables["rnahybrid"],
            )

        subprocess.run(command, check=True)
        merge_output_files(split_output_dir, output_pattern, output_path)
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


#----- Snakemake entry point: unpacks the injected `snakemake` object and calls the pure logic above -----#
def run_from_snakemake(snakemake):
    # `input.calibration` is the calibration JSON path for the calibrated variant
    # and an empty list for the uncalibrated variant — coerce the empty case to None.
    calibration_input = getattr(snakemake.input, "calibration", None)
    distribution_file = calibration_input if calibration_input else None

    run_rnahybrid(
        query=snakemake.input.query,
        target=snakemake.input.target,
        output_file=snakemake.output.compact,
        threads=snakemake.threads,
        species=snakemake.params.species,
        max_suboptimal_hits=snakemake.params.max_suboptimal_hits,
        max_internal_loop=snakemake.params.max_internal_loop,
        max_bulge_loop=snakemake.params.max_bulge_loop,
        max_total_energy=snakemake.params.max_total_energy,
        pvalue_threshold=snakemake.params.pvalue_threshold,
        seed=snakemake.params.seed,
        distribution=snakemake.params.distribution,
        max_target_length=snakemake.params.max_target_length,
        distribution_file=distribution_file,
    )

run_from_snakemake(snakemake)