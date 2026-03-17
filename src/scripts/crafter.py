# Command-line interface.
from argparse import ArgumentParser
from pathlib import Path
import shlex
import shutil
import subprocess

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
def write_fasta_chunks(target_file, chunk_prefix, max_lines=650):
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


# Pipeline execution.
def main():
    args = parse_args()

    for executable in ("parallel", "RNAhybrid"):
        if shutil.which(executable) is None:
            raise RuntimeError(f"Required executable not found in PATH: {executable}")

    chunk_prefix = "chunk_"
    output_pattern = "output_chunk_*.tsv"
    optional_args = []

    for path in Path("../..").glob(output_pattern):
        path.unlink()
    for path in Path("../..").glob(f"{chunk_prefix}*"):
        path.unlink()

    if args.b is not None:optional_args.extend(["-b", str(args.b)])
    if args.u is not None:optional_args.extend(["-u", str(args.u)])
    if args.v is not None:optional_args.extend(["-v", str(args.v)])
    if args.energy is not None:optional_args.extend(["-e", str(args.energy)])
    if args.pvalue is not None:optional_args.extend(["-p", str(args.pvalue)])

    chunk_paths = write_fasta_chunks(args.target, chunk_prefix)
    optional_arg_string = " ".join(shlex.quote(arg) for arg in optional_args)
    command_template = (
        f"RNAhybrid -q {{1}} -t {{2}} -s {{3}} -c -m 50000 {optional_arg_string} > output_{{2/.}}.tsv"
    ).strip()
    command_run = " ".join(
        [
            "parallel",
            "-j",
            shlex.quote(str(args.cores)),
            shlex.quote(command_template),
            ":::",
            shlex.quote(args.query),
            ":::",
            *[shlex.quote(str(chunk_path)) for chunk_path in chunk_paths],
            ":::",
            shlex.quote(args.s),
        ]
    )

    subprocess.run(command_run, shell=True, check=True)

    output_files = sorted(Path("../..").glob(output_pattern))  # Collect per-chunk RNAhybrid outputs in stable order.
    if not output_files:
        raise RuntimeError("No output files were produced by RNAhybrid.")  # Fail if the parallel run produced nothing.

    with open(args.output_file, "w") as merged_output:  # Open the final merged output file.
        for output_file in output_files:
            merged_output.write(output_file.read_text())  # Append each chunk result to the final output.

    for path in Path("../..").glob(output_pattern):
        path.unlink()  # Remove temporary per-chunk RNAhybrid output files.
    for path in Path("../..").glob(f"{chunk_prefix}*"):
        path.unlink()  # Remove temporary FASTA chunk files.


if __name__ == "__main__":
    main()
