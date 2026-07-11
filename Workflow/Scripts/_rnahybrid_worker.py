"""Worker invoked by GNU Parallel: runs ONE RNAhybrid call and writes its output
to a file derived from the input basenames — no shell redirection, so paths with
spaces or shell metacharacters can never break the command line.

The output directory is passed via the RNAHYBRID_OUT_DIR environment variable
(never on the command line) so a spaced results path is not word-split by the shell."""
from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path

from _common import which_required, ensure_parent


#----- Deterministic per-job output filename from the chunk (and, when calibrated, query) basename -----#
def out_path_for(out_dir, query_path, chunk_path, calibrated):
    chunk_name = Path(chunk_path).name
    if calibrated:
        name = f"output_{Path(query_path).stem}__{chunk_name}.tsv"
    else:
        name = f"output_{chunk_name}.tsv"
    return Path(out_dir) / name


#----- Builds the single RNAhybrid command; -d (calibrated) wins over -s (species) -----#
def build_rnahybrid_command(rnahybrid, query, target, max_target_length, species=None, distribution=None, extra_args=()):
    cmd = [rnahybrid, "-q", query, "-t", target, "-c", "-m", str(max_target_length)]
    if distribution is not None:
        cmd += ["-d", distribution]
    elif species is not None:
        cmd += ["-s", species]
    cmd += list(extra_args)
    return cmd


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--query", required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument("--max-target-length", required=True)
    parser.add_argument("--species")
    parser.add_argument("--dist")
    parser.add_argument("rnahybrid_args", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)

    calibrated = args.dist is not None
    out_file = ensure_parent(out_path_for(os.environ["RNAHYBRID_OUT_DIR"], args.query, args.target, calibrated))

    extra = args.rnahybrid_args
    if extra and extra[0] == "--":            # argparse REMAINDER keeps the separator
        extra = extra[1:]

    cmd = build_rnahybrid_command(
        which_required("RNAhybrid", "rnahybrid", label="RNAhybrid"),
        args.query, args.target, args.max_target_length,
        species=args.species, distribution=args.dist, extra_args=extra,
    )
    with open(out_file, "w") as fh:
        subprocess.run(cmd, check=True, stdout=fh)


if __name__ == "__main__":
    main()
