from pathlib import Path
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

from _common import which_required, query_key, ensure_parent
from _rnahybrid_worker import OUTPUT_GLOB_BROADCAST, OUTPUT_GLOB_CALIBRATED


#----- Yields one FASTA record at a time as its ORIGINAL lines. Both splitters below need the
#      verbatim text (and its line count), which _common.iter_fasta_records discards when it
#      joins the sequence — so this is the raw-text counterpart, not a duplicate of it -----#
def _iter_raw_records(path):
    record = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">") and record:
                yield record
                record = []
            record.append(line)
    if record:
        yield record


#----- Splits the target FASTA into chunk files of bounded line count, never breaking a record across chunks -----#
def write_fasta_chunks(target_file, chunk_dir, chunk_prefix="chunk_", max_lines=800):
    chunk_dir = Path(chunk_dir)
    chunk_dir.mkdir(parents=True, exist_ok=True)

    chunk_paths = []
    buffered = []

    #----- Writes the buffered records out as the next numbered chunk file -----#
    def flush_chunk():
        if not buffered:
            return
        # Zero-padded so merge_output_files' sorted(glob) orders chunks numerically rather
        # than lexicographically (chunk_10 before chunk_9).
        chunk_path = chunk_dir / f"{chunk_prefix}{len(chunk_paths):06d}"
        chunk_path.write_text("".join(buffered))
        chunk_paths.append(chunk_path)
        buffered.clear()

    for record in _iter_raw_records(target_file):
        if buffered and len(buffered) + len(record) > max_lines:
            flush_chunk()
        buffered.extend(record)
    flush_chunk()

    if not chunk_paths:
        raise RuntimeError(f"No FASTA records found in target file: {target_file}")

    return chunk_paths


_UNSAFE_STEM = re.compile(r"[^A-Za-z0-9._-]")

#----- Maps a FASTA header token to a filesystem-safe, collision-free file stem -----#
def _safe_file_stem(name, used):
    stem = _UNSAFE_STEM.sub("_", name) or "query"
    candidate, n = stem, 1
    while candidate in used:
        candidate, n = f"{stem}_{n}", n + 1
    used.add(candidate)
    return candidate


#----- Explodes a multi-record query FASTA into one file per miRNA, keyed by the first header token -----#
def split_query_per_miRNA(query_file, query_dir):
    query_dir = Path(query_dir)
    query_dir.mkdir(parents=True, exist_ok=True)

    paths_by_name = {}
    used_stems = set()

    for record in _iter_raw_records(query_file):
        if not record[0].startswith(">"):
            continue                      # anything before the first header is not a record
        # query_key: the same first-whitespace-token rule RNAcalibrate prints in column 1.
        name = query_key(record[0][1:].strip())
        if name in paths_by_name:
            raise RuntimeError(f"Duplicate query ID in FASTA: {name}")
        path = query_dir / f"{_safe_file_stem(name, used_stems)}.fa"
        path.write_text("".join(record))
        paths_by_name[name] = path

    if not paths_by_name:
        raise RuntimeError(f"No FASTA records found in query file: {query_file}")

    return paths_by_name


#----- Returns the GNU Parallel binary, having first confirmed RNAhybrid exists. RNAhybrid is
#      resolved again per-invocation inside _rnahybrid_worker.py; the check here (path
#      discarded) is only to fail fast, before dispatching the parallel fan-out -----#
def ensure_dependencies():
    which_required("RNAhybrid", "rnahybrid")
    return which_required("parallel")


#----- Translates the optional RNAhybrid params into CLI flags, omitting any that are None -----#
def build_optional_args(max_suboptimal_hits=None, max_internal_loop=None, max_bulge_loop=None,
                        max_hybrid_energy=None, pvalue_threshold=None, seed=None,
                        distribution=None):
    # (flag, value) in emission order; each pair is appended only when set.
    flag_values = [
        ("-b", max_suboptimal_hits),
        ("-u", max_internal_loop),
        ("-v", max_bulge_loop),
        ("-e", max_hybrid_energy),
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


#----- GNU Parallel command (uncalibrated): worker per (query × chunk); paths travel as {n} args / env, never argv0 -----#
def build_parallel_command_broadcast(query, chunk_paths, species, optional_args, threads,
                                     max_target_length, parallel_executable, worker_executable):
    worker_call = " ".join([
        "python3", "{1}",
        "--query", "{2}",
        "--target", "{3}",
        "--max-target-length", str(max_target_length),
        *(["--species", "{4}"] if species is not None else []),
        "--", *optional_args,
    ])
    command = [
        parallel_executable, f"-j{max(1, threads)}",
        worker_call,
        ":::", str(worker_executable),
        ":::", str(query),
        ":::", *[str(chunk_path) for chunk_path in chunk_paths],
    ]
    if species is not None:
        command.extend([":::", species])
    return command


#----- GNU Parallel command (calibrated): each TSV row is (query, distribution, chunk); paths as {n} args / env -----#
def build_parallel_command_calibrated(spec_path, optional_args, threads, max_target_length,
                                      parallel_executable, worker_executable):
    worker_call = " ".join([
        "python3", "{1}",
        "--query", "{2}",
        "--dist", "{3}",
        "--target", "{4}",
        "--max-target-length", str(max_target_length),
        "--", *optional_args,
    ])
    command = [
        parallel_executable, f"-j{max(1, threads)}",
        "--colsep", "\\t",
        worker_call,
        ":::", str(worker_executable),
        "::::", str(spec_path),
    ]
    return command


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
def run_rnahybrid(query, target, species, output_file, max_target_length, threads=1,
                  max_suboptimal_hits=None, max_internal_loop=None, max_bulge_loop=None,
                  max_hybrid_energy=None, pvalue_threshold=None, seed=None,
                  distribution=None, distribution_file=None):
    output_path = ensure_parent(output_file)

    parallel_executable = ensure_dependencies()
    worker_executable = str(Path(__file__).resolve().parent / "_rnahybrid_worker.py")

    # Per-sample tmp dir next to the final output so parallel samples never collide.
    tmp_dir = Path(tempfile.mkdtemp(prefix=".rnahybrid_", dir=output_path.parent))
    chunk_dir = tmp_dir / "chunks"
    split_output_dir = tmp_dir / "outputs"
    split_output_dir.mkdir(parents=True, exist_ok=True)

    try:
        chunk_paths = write_fasta_chunks(target, chunk_dir)
        dist_map = load_per_query_distributions(distribution_file)
        calibrated = dist_map is not None

        # One call for both branches. On the calibrated path the static `distribution` is
        # dropped explicitly rather than by omission: calibration supplies per-miRNA xi/theta
        # via the job spec's -d column, which supersedes both it and `species`.
        optional_args = build_optional_args(
            max_suboptimal_hits=max_suboptimal_hits, max_internal_loop=max_internal_loop,
            max_bulge_loop=max_bulge_loop, max_hybrid_energy=max_hybrid_energy,
            pvalue_threshold=pvalue_threshold, seed=seed,
            distribution=None if calibrated else distribution,
        )

        if calibrated:
            # One RNAhybrid invocation per (miRNA, chunk), each with its own -d.
            query_paths_by_name = split_query_per_miRNA(query, tmp_dir / "queries")
            spec_path = build_job_spec_tsv(
                query_paths_by_name=query_paths_by_name,
                dist_map=dist_map,
                chunk_paths=chunk_paths,
                tsv_path=tmp_dir / "job_spec.tsv",
            )
            command = build_parallel_command_calibrated(
                spec_path=spec_path,
                optional_args=optional_args,
                threads=threads,
                max_target_length=max_target_length,
                parallel_executable=parallel_executable,
                worker_executable=worker_executable,
            )
        else:
            # Broadcast/uncalibrated branch: one shared query across all chunks.
            species = validate_rnahybrid_args(species=species, distribution=distribution)
            command = build_parallel_command_broadcast(
                query=query,
                chunk_paths=chunk_paths,
                species=species,
                optional_args=optional_args,
                threads=threads,
                max_target_length=max_target_length,
                parallel_executable=parallel_executable,
                worker_executable=worker_executable,
            )

        # Fixed by the branch above, so it is chosen here rather than returned by both builders.
        output_pattern = OUTPUT_GLOB_CALIBRATED if calibrated else OUTPUT_GLOB_BROADCAST

        env = {**os.environ, "RNAHYBRID_OUT_DIR": str(split_output_dir)}
        subprocess.run(command, check=True, env=env)
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
        max_hybrid_energy=snakemake.params.max_hybrid_energy,
        pvalue_threshold=snakemake.params.pvalue_threshold,
        seed=snakemake.params.seed,
        distribution=snakemake.params.distribution,
        max_target_length=snakemake.params.max_target_length,
        distribution_file=distribution_file,
    )

if "snakemake" in globals():
    run_from_snakemake(snakemake)