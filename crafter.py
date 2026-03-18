# Command-line interface.
from argparse import ArgumentParser  # Build the CLI argument parser.
from pathlib import Path  # Work with filesystem paths in an OS-agnostic way.
import shlex  # Safely quote shell arguments before building a shell command.
import shutil  # Check whether required executables are available in PATH.
import subprocess  # Run the assembled external command pipeline.

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('-q', '--query', type=str, required=True, help='query file')
    parser.add_argument('-t', '--target', type=str, required=True, help='target file')
    parser.add_argument('-s', type=str, required=True, help='(3utr_fly|3utr_worm|3utr_human)')
    parser.add_argument('-j', '--cores', type=int, required=True, help='<number of cores>')
    parser.add_argument('output_file', help='output file')
    parser.add_argument('-b', type=int, required=False, help='<number of hits per target>')
    parser.add_argument('-u', type=int, required=False, help='<max internal loop size (per side)>')
    parser.add_argument('-v', type=int, required=False, help='<max bulge loop size>')
    parser.add_argument('-e', '--energy', type=int, required=False, help='<energy cut-off>')
    parser.add_argument('-p', '--pvalue', type=float, required=False, help='<p-value cut-off>')
    return parser.parse_args()


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


# Pipeline execution.
def main():  # Parse CLI inputs, run RNAhybrid in parallel, and merge temporary outputs.
    args = parse_args()  # Read user-supplied command-line arguments.

    for executable in ("parallel", "RNAhybrid"):  # Verify that both external dependencies are installed.
        if shutil.which(executable) is None:  # Search PATH for the current executable name.
            raise RuntimeError(f"Required executable not found in PATH: {executable}")  # Stop before doing any work if a dependency is missing.

    chunk_prefix = "chunk_"  # Prefix used for temporary FASTA chunk filenames.
    output_pattern = "output_chunk_*.tsv"  # Glob pattern matching temporary RNAhybrid outputs.
    optional_args = []  # Accumulate optional CLI flags that should be forwarded to RNAhybrid.

    for path in Path(".").glob(output_pattern):  # Remove stale output files from previous runs.
        path.unlink()  # Delete the matching temporary output file.
    for path in Path(".").glob(f"{chunk_prefix}*"):  # Remove stale FASTA chunk files from previous runs.
        path.unlink()  # Delete the matching temporary chunk file.

    if args.b is not None: optional_args.extend(["-b", str(args.b)])  # Forward the optional hit-count limit if provided.
    if args.u is not None: optional_args.extend(["-u", str(args.u)])  # Forward the optional internal-loop limit if provided.
    if args.v is not None: optional_args.extend(["-v", str(args.v)])  # Forward the optional bulge-loop limit if provided.
    if args.energy is not None: optional_args.extend(["-e", str(args.energy)])  # Forward the optional energy cutoff if provided.
    if args.pvalue is not None: optional_args.extend(["-p", str(args.pvalue)])  # Forward the optional p-value cutoff if provided.

    chunk_paths = write_fasta_chunks(args.target, chunk_prefix)  # Split the target FASTA into manageable chunk files.
    optional_arg_string = " ".join(shlex.quote(arg) for arg in optional_args)  # Turn optional arguments into a shell-safe string.
    command_template = (  # Build the per-job GNU parallel command template.
        f"RNAhybrid -q {{1}} -t {{2}} -s {{3}} -c -m 500000 {optional_arg_string} > output_{{2/.}}.tsv"  # Run RNAhybrid on one chunk and redirect to a chunk-specific TSV.
    ).strip()  # Remove extra whitespace if no optional arguments were added.
    command_run = " ".join(  # Assemble the full GNU parallel command line as a single shell string.
        [  # Build the command from individually quoted shell tokens.
            "parallel",  # Invoke GNU parallel.
            "-j",  # Specify the parallel job count flag.
            shlex.quote(str(args.cores)),  # Pass the requested number of concurrent jobs.
            shlex.quote(command_template),  # Pass the shell-quoted job template to GNU parallel.
            ":::",  # Start the first input source list for placeholder {1}.
            shlex.quote(args.query),  # Provide the same query file to every parallel job.
            ":::",  # Start the second input source list for placeholder {2}.
            *[shlex.quote(str(chunk_path)) for chunk_path in chunk_paths],  # Provide every generated target chunk as a separate job value.
            ":::",  # Start the third input source list for placeholder {3}.
            shlex.quote(args.s),  # Provide the selected species/model argument to every job.
        ]
    )

    subprocess.run(command_run, shell=True, check=True)  # Execute the GNU parallel command and raise if it fails.

    output_files = sorted(Path(".").glob(output_pattern))  # Collect per-chunk RNAhybrid outputs in stable order.
    if not output_files:  # Guard against silent failures that produced no output files.
        raise RuntimeError("No output files were produced by RNAhybrid.")  # Fail if the parallel run produced nothing.

    with open(args.output_file, "w") as merged_output:  # Open the final merged output file.
        for output_file in output_files:  # Visit each temporary RNAhybrid output file in order.
            merged_output.write(output_file.read_text())  # Append each chunk result to the final output.

    for path in Path(".").glob(output_pattern):  # Clean up all temporary RNAhybrid output files.
        path.unlink()  # Remove the temporary per-chunk RNAhybrid output file.
    for path in Path(".").glob(f"{chunk_prefix}*"):  # Clean up all temporary FASTA chunk files.
        path.unlink()  # Remove the temporary FASTA chunk file.


if __name__ == "__main__":  # Run the CLI entry point only when the file is executed as a script.
    main()  # Start the program.
