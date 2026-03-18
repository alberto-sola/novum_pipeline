from pathlib import Path  # Work with filesystem paths in an OS-agnostic way.
import shlex  # Safely quote shell arguments before building a shell command.
import shutil  # Check whether required executables are available in PATH.
import subprocess  # Run the assembled external command pipeline.


# FASTA chunking.
def write_fasta_chunks(target_file, chunk_prefix, max_lines=650):  # Split a FASTA file into smaller temporary files.
    chunk_paths = []  # Collect the paths of all chunk files that get written.
    current_chunk_index = 0  # Track which chunk file number we are currently building.
    current_chunk_line_total = 0  # Count how many FASTA lines are already in the current chunk.
    current_record = []  # Buffer the FASTA record currently being read from disk.
    current_chunk_records = []  # Buffer all record lines assigned to the current chunk.

    def flush_chunk():  # Persist the current chunk buffer to disk and reset the chunk state.
        nonlocal current_chunk_index, current_chunk_line_total, current_chunk_records  # Rebind outer-scope chunk state variables.
        if not current_chunk_records:  # Skip writing when the current chunk has no buffered records.
            return  # Nothing to flush.

        first_letter = chr(ord("a") + ((current_chunk_index // 26) % 26))  # Derive the first suffix letter from the chunk index.
        second_letter = chr(ord("a") + (current_chunk_index % 26))  # Derive the second suffix letter from the chunk index.
        chunk_path = Path(f"{chunk_prefix}{first_letter}{second_letter}")  # Build a filename like chunk_aa or chunk_ab.
        chunk_path.write_text("".join(current_chunk_records))  # Write the buffered FASTA records into the chunk file.
        chunk_paths.append(chunk_path)  # Remember the path so the caller can process the chunk later.
        current_chunk_index += 1  # Advance to the next chunk name.
        current_chunk_line_total = 0  # Reset the line counter for the next chunk.
        current_chunk_records = []  # Clear the in-memory chunk buffer.

    def add_record(record_lines):  # Add one parsed FASTA record to the active chunk buffer.
        nonlocal current_chunk_line_total, current_chunk_records  # Update the outer chunk buffer and line count.
        record_line_total = len(record_lines)  # Measure how many lines this FASTA record occupies.
        if current_chunk_records and current_chunk_line_total + record_line_total > max_lines:  # Flush first if this record would overflow the chunk.
            flush_chunk()  # Start a fresh chunk before adding the record.
        current_chunk_records.extend(record_lines)  # Append the full record to the current chunk buffer.
        current_chunk_line_total += record_line_total  # Increase the line count to reflect the appended record.

    with open(target_file) as handle:  # Stream through the target FASTA file line by line.
        for line in handle:  # Process each raw line from the FASTA file.
            if line.startswith(">"):  # A header line marks the beginning of a new FASTA record.
                if current_record:  # If a previous record was buffered, finish assigning it to a chunk.
                    add_record(current_record)  # Store the completed record in the current chunk.
                current_record = [line]  # Start buffering the new record with its header line.
            else:  # Non-header lines belong to the current FASTA sequence body.
                current_record.append(line)  # Add the sequence line to the active record buffer.

    if current_record:  # After the read loop, there may still be one final buffered record.
        add_record(current_record)  # Add the last record to a chunk.
    flush_chunk()  # Write any remaining buffered chunk data to disk.

    if not chunk_paths:  # An empty chunk list means the target file had no FASTA records.
        raise RuntimeError(f"No FASTA records found in target file: {target_file}")  # Fail fast with a clear error.

    return chunk_paths  # Return all generated chunk paths for downstream processing.


def ensure_dependencies():  # Verify that required external executables are available before starting work.
    for executable in ("parallel", "RNAhybrid"):
        if shutil.which(executable) is None:
            raise RuntimeError(f"Required executable not found in PATH: {executable}")


def cleanup_temp_files(chunk_prefix, output_pattern):  # Remove stale chunk and RNAhybrid output files from the working directory.
    for path in Path(".").glob(output_pattern):
        path.unlink()
    for path in Path(".").glob(f"{chunk_prefix}*"):
        path.unlink()


def build_optional_args(hits=None, u=None, v=None, energy=None, pvalue=None, seed=None):  # Translate optional crafter settings into RNAhybrid CLI flags.
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


def build_parallel_command(query, chunk_paths, species, cores, optional_args):  # Build the GNU parallel shell command that fans out RNAhybrid across target chunks.
    optional_arg_string = " ".join(shlex.quote(arg) for arg in optional_args)
    command_template = f"RNAhybrid -q {{1}} -t {{2}} -s {{3}} -c -m 500000 {optional_arg_string} > output_{{2/.}}.tsv".strip()

    return " ".join(
        [
            "parallel",
            "-j",
            shlex.quote(str(cores)),
            shlex.quote(command_template),
            ":::",
            shlex.quote(query),
            ":::",
            *[shlex.quote(str(chunk_path)) for chunk_path in chunk_paths],
            ":::",
            shlex.quote(species),
        ]
    )


def merge_output_files(output_pattern, merged_output_path):  # Merge all per-chunk RNAhybrid outputs into the final output file.
    output_files = sorted(Path(".").glob(output_pattern))
    if not output_files:
        raise RuntimeError("No output files were produced by RNAhybrid.")

    with open(merged_output_path, "w") as merged_output:
        for output_file in output_files:
            merged_output.write(output_file.read_text())


def run_crafter(query, target, species, cores, output_file, hits=None, u=None, v=None, energy=None, pvalue=None, seed=None):
    chunk_prefix = "chunk_"
    output_pattern = "output_chunk_*.tsv"

    ensure_dependencies()
    cleanup_temp_files(chunk_prefix, output_pattern)

    chunk_paths = write_fasta_chunks(target, chunk_prefix)
    optional_args = build_optional_args(
        hits=hits,
        u=u,
        v=v,
        energy=energy,
        pvalue=pvalue,
        seed=seed,
    )
    command_run = build_parallel_command(
        query=query,
        chunk_paths=chunk_paths,
        species=species,
        cores=cores,
        optional_args=optional_args,
    )

    subprocess.run(command_run, shell=True, check=True)
    merge_output_files(output_pattern, output_file)
    cleanup_temp_files(chunk_prefix, output_pattern)