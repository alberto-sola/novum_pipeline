import hashlib
import json
import math
import os
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from functools import partial
from pathlib import Path

from _common import which_required, iter_fasta_records, iter_fasta_headers, query_key, ensure_parent


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


#----- The one normalisation both the RNG seed and the provenance digest must agree on -----#
def _normalise_sequence(sequence):
    return sequence.strip().upper().replace("T", "U")


#----- The one hash of a sequence. sha256, never built-in hash() — that one is salted per
#      process by PYTHONHASHSEED, which would make the seed below irreproducible -----#
def _sequence_digest(sequence):
    return hashlib.sha256(_normalise_sequence(sequence).encode())


#----- Per-miRNA RNG seed: stable across processes and machines, so calibration is reproducible -----#
def derive_seed(base_seed, sequence):
    return (base_seed + int.from_bytes(_sequence_digest(sequence).digest()[:4], "big")) % 2**31


#----- Provenance digest of one sequence; shares derive_seed's normalisation -----#
def sequence_sha256(sequence):
    return _sequence_digest(sequence).hexdigest()


#----- Provenance digest of a whole input file -----#
def sha256_file(path):
    with open(path, "rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


#----- Record count for provenance; headers only, so no sequence is accumulated to be discarded -----#
def count_fasta_records(path):
    return sum(1 for _header in iter_fasta_headers(path))


#----- Identity of one input file: distinguishes "same input" from "same path" -----#
def _file_provenance(path):
    return {"path": str(path), "sha256": sha256_file(path), "n_records": count_fasta_records(path)}


#----- Provenance: distinguishes "same input" from "same path", which `command` alone cannot -----#
def build_inputs_block(query, target, k, max_target_length, length_arg, forced_helix,
                       max_internal_loop, max_bulge_loop, rng_seed):
    return {
        "query": _file_provenance(query),
        "target": _file_provenance(target),
        "params": {
            "k": k,
            "max_target_length": max_target_length,
            "length_arg": length_arg["value"],
            "forced_helix": forced_helix,
            "max_internal_loop": max_internal_loop,
            "max_bulge_loop": max_bulge_loop,
        },
        "rng_seed": rng_seed,
    }


#----- libfaketime pins RNAcalibrate's clock; the binary exposes no seed flag, so the seed
#      must be injected from outside the process. macOS differs in BOTH the variable name and the library filename. -----#
FAKETIME_LIBRARIES = {
    "linux": "libfaketime.so.1",
    "darwin": "libfaketime.1.dylib",
}


#----- Locates the interposition library beside the `faketime` wrapper on PATH (never a hardcoded CONDA_PREFIX) -----#
def resolve_faketime_library(platform_name):
    library_name = FAKETIME_LIBRARIES.get(platform_name)
    if library_name is None:
        raise RuntimeError(
            f"Unsupported platform for RNG seed pinning: {platform_name}. "
            "Set rnacalibrate.rng_seed to null to run unpinned."
        )

    try:
        wrapper = which_required("faketime")
    except RuntimeError as exc:
        raise RuntimeError(
            "rnacalibrate.rng_seed is set but libfaketime is not installed. "
            "Add `libfaketime` to Workflow/Envs/rnahybrid.yaml, or set rng_seed to null."
        ) from exc

    library = Path(wrapper).resolve().parent.parent / "lib" / "faketime" / library_name
    if not library.exists():
        raise RuntimeError(
            f"Found the faketime wrapper at {wrapper} but not {library_name} at {library}. "
            "Reinstall libfaketime in Workflow/Envs/rnahybrid.yaml."
        )
    return library


#----- Environment overlay that freezes the clock at `seed`; '@' means absolute-and-frozen, not an offset -----#
def build_faketime_env(seed, library_path, platform_name, base_env):
    if platform_name not in FAKETIME_LIBRARIES:
        raise RuntimeError(f"Unsupported platform for RNG seed pinning: {platform_name}")

    env = dict(base_env)
    env["FAKETIME_FMT"] = "%s"
    env["FAKETIME"] = f"@{seed}"
    if platform_name == "darwin":
        env["DYLD_INSERT_LIBRARIES"] = str(library_path)
        env["DYLD_FORCE_FLAT_NAMESPACE"] = "1"
    else:
        env["LD_PRELOAD"] = str(library_path)
    return env


#----- Cheap probe rungs. Below the fit's convergence threshold RNAcalibrate returns
#      `-nan -nan`, a CONSTANT indistinguishable from "the seed has no effect" — measured,
#      hsa-miR-2861 (19 nt, 89% GC) is NaN at every seed until k~1500 — so neither rung can
#      be the last one, and probe_ladder always ends at the run's own k -----#
PROBE_K = 5
PROBE_K_RETRY = 200


#----- Cheap rungs only while they are cheaper than the real run; k itself is always last, so
#      the verdict is never decided in a regime the run never enters -----#
def probe_ladder(k):
    return [rung for rung in (PROBE_K, PROBE_K_RETRY) if rung < k] + [k]


#----- Pinning only means anything when `-s` is on, because that is the sole RNG consumer -----#
def pinning_enabled(rng_seed, randomize_targets):
    return rng_seed is not None and bool(randomize_targets)


#----- Three-outcome verdict. The third run carries the weight: two identical outputs are equally consistent with "pinned" and "probe too small to tell". -----#
def classify_probe(first, repeat, other):
    if first != repeat:
        return "blocked"
    if first == other:
        return "blind"
    return "verified"


#----- Refuses to proceed unless the clock is demonstrably pinned AND controlling the output.
#      Always probes with randomize_targets=True: `-s` is the only consumer of the RNG, so a
#      probe without it would be deterministic at every seed and report a false "blind".
#      The caller must therefore only invoke this when randomize_targets is on.
#      COST: probes run over the whole query FASTA, so the k rung is ~one full calibration
#      each. Cheap in practice (k=5 settles it unless every miRNA is degenerate there), but
#      the worst case is ~3x a real run — deliberate, since a single-record probe is far
#      likelier to report a false "blind" on a miRNA that is NaN at low k. -----#
def verify_faketime(executable, query, target, k, max_target_length, length_arg, seed,
                    library_path, platform_name):
    #----- One probe run at a given seed and k, returning its raw stdout to compare -----#
    def probe(probe_seed, probe_k):
        command = build_command(
            executable=executable, query=query, target=target, k=probe_k,
            max_target_length=max_target_length, length_arg=length_arg,
            randomize_targets=True,
        )
        env = build_faketime_env(probe_seed, library_path, platform_name, os.environ)
        try:
            return subprocess.run(command, check=True, capture_output=True, text=True, env=env).stdout
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                f"RNG seed pinning probe failed to run on {platform_name}: "
                f"{' '.join(command)} exited with code {exc.returncode}. "
                f"Error output: {exc.stderr}. Set rnacalibrate.rng_seed to null to run unpinned, "
                "accepting non-reproducible calibration."
            ) from exc

    ladder = probe_ladder(k)
    for probe_k in ladder:
        first = probe(seed, probe_k)
        # An unpinned clock reads the wall clock at 1s resolution, so back-to-back probes
        # (~0.15s apart) are byte-identical whether or not libfaketime is doing anything.
        # Sleeping past a second boundary forces an UNPINNED clock to diverge; a PINNED one
        # is frozen, so first == repeat however long we wait — no false "blocked".
        time.sleep(1.1)
        repeat = probe(seed, probe_k)
        other = probe(seed + 1, probe_k)
        verdict = classify_probe(first, repeat, other)
        if verdict == "verified":
            return
        if verdict == "blocked":
            raise RuntimeError(
                f"rnacalibrate.rng_seed is set but the RNG clock is not pinned on "
                f"{platform_name}: two runs at seed {seed} differed. Library injection was "
                "likely blocked (on macOS, System Integrity Protection strips "
                "DYLD_INSERT_LIBRARIES). Set rnacalibrate.rng_seed to null to run unpinned, "
                "accepting non-reproducible calibration."
            )

    raise RuntimeError(
        f"Could not verify RNG seed pinning on {platform_name}: output was identical at seeds {seed} and "
        f"{seed + 1} even at k={ladder[-1]}, the k this run itself uses, so the probe cannot tell whether the seed "
        "controls the result. Refusing to claim reproducibility. Set rnacalibrate.rng_seed to null to run unpinned, "
        "accepting non-reproducible calibration."
    )


#----- Assembles the RNAcalibrate command line; optional flags (max_internal_loop, max_bulge_loop, seed, randomize) are appended only if set -----#
def build_command(executable, query, target, k, max_target_length, length_arg,
                  randomize_targets=False, max_internal_loop=None, max_bulge_loop=None,
                  seed=None):
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


#----- A degenerate (-nan) fit, as opposed to malformed or absent output. Its own class so the
#      caller can attach seed/k advice by TYPE — rewording a message must not delete it -----#
class DegenerateFitError(RuntimeError):
    pass


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
            raise DegenerateFitError(f"RNAcalibrate produced NaN parameters: {raw_line}")

        per_query.append({
            "query": fields[0],
            "sample_size": int(fields[1]),
            "xi": xi,
            "theta": theta,
        })

    if not per_query:
        raise RuntimeError("RNAcalibrate did not produce any calibration rows.")

    return {"per_query": per_query}


#----- Duplicate IDs would silently collapse in rnahybrid.py's load_per_query_distributions -----#
def iter_query_records(query_file):
    seen = set()

    for header, sequence in iter_fasta_records(query_file):
        name = query_key(header)
        if name in seen:
            raise RuntimeError(f"Duplicate query ID in FASTA: {name} ({query_file})")
        seen.add(name)
        yield header, sequence

    if not seen:
        raise RuntimeError(f"No FASTA records found in query file: {query_file}")


#----- A single-record invocation must yield exactly one fit; more means the split leaked -----#
def expect_single_row(parsed, mirna):
    rows = parsed["per_query"]
    if len(rows) != 1:
        raise RuntimeError(
            f"Expected exactly 1 calibration row for {mirna}, got {len(rows)}."
        )
    return rows[0]


#----- One miRNA's fit: its own single-record -q file, so it is always position 1 in
#      RNAcalibrate's RNG stream, and its own derived seed. Nothing here reads another
#      miRNA's state — that independence is what lets the driver run these concurrently -----#
def _calibrate_one(index, header, sequence, tmpdir, make_command,
                   rng_seed, library_path, platform_name, pin_clock):
    mirna = query_key(header)
    # Indexed, not named after the miRNA: two distinct IDs may carry the same sequence, and
    # concurrent jobs must never share a path.
    single_record = Path(tmpdir) / f"query_{index:05d}.fa"
    single_record.write_text(f">{header}\n{sequence}\n")
    command = make_command(query=str(single_record))

    derived_seed = None
    env = None
    if pin_clock:
        derived_seed = derive_seed(rng_seed, sequence)
        env = build_faketime_env(derived_seed, library_path, platform_name, os.environ)

    try:
        completed = subprocess.run(command, check=True, capture_output=True, text=True, env=env)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"RNAcalibrate failed for miRNA {mirna}: "
            f"{' '.join(command)} exited with code {exc.returncode}. "
            f"Error output: {exc.stderr}."
        ) from exc

    context = f"[miRNA {mirna}, derived_seed={derived_seed}]"
    try:
        parsed = parse_rnacalibrate_output(completed.stdout)
    except DegenerateFitError as exc:
        # Caught by TYPE, not by matching the message: a degenerate fit at a pinned seed is
        # reproducible, so re-running will not clear it.
        raise DegenerateFitError(
            f"{exc} {context}. A degenerate fit at a pinned seed is reproducible, so "
            "re-running will not clear it — change rnacalibrate.rng_seed or raise "
            "rnacalibrate.k."
        ) from exc
    except (RuntimeError, ValueError) as exc:
        # Malformed or absent output — the pinning/k advice above would misdirect here.
        raise RuntimeError(f"{exc} {context}.") from exc

    entry = expect_single_row(parsed, mirna)
    entry["sequence_sha256"] = sequence_sha256(sequence)
    entry["derived_seed"] = derived_seed
    return entry


#----- Top-level driver: one RNAcalibrate invocation PER miRNA. NOT an optimisation —
#      RNAcalibrate consumes ONE RNG stream across every query in a -q file, so a miRNA's
#      xi/theta depends on how many records precede it (measured, purely ordinal). Batching
#      these back into one call reintroduces that dependence and makes results
#      non-comparable across differing query sets. Costs ~0.9% (measured, k=2000).
#      The invocations are independent, so they run on `threads` cores at once; pool.map
#      preserves input order, keeping per_query in query-FASTA order -----#
def run_rnacalibrate(query, target, output_file, k, max_target_length, randomize_targets=False,
                     max_internal_loop=None, max_bulge_loop=None, seed=None, rng_seed=None,
                     threads=1):
    executable = which_required("RNAcalibrate", "rnacalibrate")
    stats = compute_target_length_stats(target)
    length_arg = build_length_arg(stats)
    platform_name = sys.platform

    pin_clock = pinning_enabled(rng_seed, randomize_targets)
    if rng_seed is not None and not randomize_targets:
        # Legal but inert (see pinning_enabled) — warn rather than raise; provenance
        # already records inputs.rng_seed as null here.
        print(
            f"rnacalibrate: rng_seed={rng_seed} is set but randomize_targets is false, so "
            "RNAcalibrate's -s flag (the RNG's only consumer) is never passed; the seed has "
            "no effect and this run is exactly as (non-)reproducible as rng_seed: null.",
            file=sys.stderr,
        )

    library_path = None
    if pin_clock:
        library_path = resolve_faketime_library(platform_name)
        verify_faketime(executable, query, target, k, max_target_length, length_arg,
                        rng_seed, library_path, platform_name)

    # One binding for every argv the run produces, so the per-miRNA invocations and the
    # provenance template below can never drift apart on a flag.
    make_command = partial(
        build_command,
        executable=executable, target=target, k=k, max_target_length=max_target_length,
        length_arg=length_arg, randomize_targets=randomize_targets,
        max_internal_loop=max_internal_loop, max_bulge_loop=max_bulge_loop, seed=seed,
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        # Materialised before dispatch so the duplicate-ID and empty-file guards in
        # iter_query_records still raise up front, not inside a worker.
        records = list(iter_query_records(query))
        calibrate = partial(
            _calibrate_one, tmpdir=tmpdir, make_command=make_command, rng_seed=rng_seed,
            library_path=library_path, platform_name=platform_name, pin_clock=pin_clock,
        )
        with ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
            per_query = list(pool.map(
                lambda job: calibrate(job[0], *job[1]), enumerate(records)
            ))

    ensure_parent(output_file).write_text(json.dumps(
        {
            # Informational only (nothing reads it): the argv invariant across the N
            # invocations, rather than any one miRNA's temp path.
            "command": make_command(query="<per-query>"),
            "target_length_stats": stats,
            "target_length_argument": length_arg,
            "inputs": build_inputs_block(
                query=query, target=target, k=k, max_target_length=max_target_length,
                length_arg=length_arg, forced_helix=seed,
                max_internal_loop=max_internal_loop, max_bulge_loop=max_bulge_loop,
                # The seed that was actually APPLIED, so provenance never claims a
                # pinned run that did not happen (e.g. randomize_targets: false).
                rng_seed=rng_seed if pin_clock else None,
            ),
            "calibration": {"per_query": per_query},
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
        rng_seed=snakemake.params.rng_seed,
        threads=snakemake.threads,
    )

if "snakemake" in globals():
    run_from_snakemake(snakemake)