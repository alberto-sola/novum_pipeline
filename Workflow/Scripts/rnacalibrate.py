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
from _rnacalibrate_config import assign_anchor, cell_edges, fitted_anchors, parse_anchors
from _gumbel_fit import (
    FitPoint,
    UnusableCurveError,
    classify_fit,
    fit_alpha_robust,
    predict_theta,
    residual_ratio,
    weighted_median,
)


#----- The lengths every downstream derivation reads: the `-l` stats and the cell occupancy
#      must describe the same population, so they share one pass and one `if seq` filter -----#
def record_lengths(target_file):
    lengths = [len(seq) for _header, seq in iter_fasta_records(target_file) if seq]
    if not lengths:
        raise RuntimeError(f"No FASTA records found in target file: {target_file}")
    return lengths


#----- Mean and population stdev of already-read record lengths (the input to RNAcalibrate's `-l`) -----#
def summarise_lengths(lengths):
    mean_length = sum(lengths) / len(lengths)
    if len(lengths) == 1:
        std_length = 0.0
    else:
        variance = sum((length - mean_length) ** 2 for length in lengths) / len(lengths)
        std_length = math.sqrt(variance)

    return {"count": len(lengths), "mean": mean_length, "std": std_length}


#----- The same, straight from a FASTA -----#
def compute_target_length_stats(target_file):
    return summarise_lengths(record_lengths(target_file))


#----- Rounds the (mean, std) pair to ints and formats them as RNAcalibrate's "<mean>,<std>" -l argument -----#
def build_length_arg(stats):
    mean_length = int(round(stats["mean"]))
    std_length = int(round(stats["std"]))
    return {
        "mean": mean_length,
        "std": std_length,
        "value": f"{mean_length},{std_length}",
    }


#----- The -l for one anchor. Fixing std/mean at 1/3 is the point: the null's WIDTH then
#      depends on the cell, not on how a genome's lengths fall inside it. divisor=2 is
#      Tier 2's widening rung -----#
def anchor_length_arg(anchor, divisor=3):
    std_length = int(round(anchor / divisor))
    return {"mean": anchor, "std": std_length, "value": f"{anchor},{std_length}"}


#----- Cell occupancy. An empty anchor is never fitted and gets no stratum — a genome with
#      no genes under 107 nt simply has no 76 cell -----#
def count_records_per_anchor(lengths, anchors):
    edges = cell_edges(anchors)
    counts = {anchor: 0 for anchor in anchors}
    for length in lengths:
        counts[assign_anchor(length, anchors, edges)] += 1
    return counts


#----- The one normalisation both the RNG seed and the provenance digest must agree on -----#
def _normalise_sequence(sequence):
    return sequence.strip().upper().replace("T", "U")


#----- The one hash of a sequence. sha256, never built-in hash() — that one is salted per
#      process by PYTHONHASHSEED, which would make the seed below irreproducible -----#
def _sequence_digest(sequence):
    return hashlib.sha256(_normalise_sequence(sequence).encode())


#----- Per-miRNA RNG seed, stable across processes and machines. With no anchor it returns
#      exactly what it always has, keeping the reference fit comparable with every run on
#      disk; an anchor folds in the cell and attempt, so strata never share a draw and a
#      re-draw is as reproducible as the first try -----#
def derive_seed(base_seed, sequence, anchor=None, attempt=0):
    digest = _sequence_digest(sequence)
    if anchor is not None:
        digest.update(f"|{anchor}|{attempt}".encode())
    return (base_seed + int.from_bytes(digest.digest()[:4], "big")) % 2**31


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
                       max_internal_loop, max_bulge_loop, rng_seed, length_anchors=None):
    return {
        "query": _file_provenance(query),
        "target": _file_provenance(target),
        "params": {
            "k": k,
            "max_target_length": max_target_length,
            "length_arg": length_arg["value"],
            "length_anchors": list(length_anchors) if length_anchors else None,
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

# Just past a wall-clock second. RNAcalibrate seeds from time() at 1s resolution, so this is
# what forces an UNPINNED clock to hand out a different draw; a pinned one is unaffected.
CLOCK_TICK_SLEEP_S = 1.1


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
        time.sleep(CLOCK_TICK_SLEEP_S)
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
def parse_rnacalibrate_output(stdout, reject_nan=True):
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
        # NaN means a degenerate sample; keep it out of the downstream -d. The anchor path
        # passes reject_nan=False: it classifies and retries rather than aborting.
        if reject_nan and (math.isnan(xi) or math.isnan(theta)):
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


#----- This miRNA's single-record -q file, so every invocation is position 1 in its own RNG
#      stream. Written ONCE in the main thread and only read after, so the anchor fits share
#      it without racing. Indexed, not named after the miRNA: two IDs may share a sequence -----#
def _write_query_file(tmpdir, index, header, sequence):
    path = Path(tmpdir) / f"query_{index:05d}.fa"
    path.write_text(f">{header}\n{sequence}\n")
    return path


#----- The derived seed and the faketime env for one fit; `seed_parts` separates the anchor
#      path's (anchor, attempt) draws. Both None when the clock is not pinned -----#
def _fit_env(pin_clock, rng_seed, sequence, library_path, platform_name, **seed_parts):
    if not pin_clock:
        return None, None
    derived_seed = derive_seed(rng_seed, sequence, **seed_parts)
    return derived_seed, build_faketime_env(derived_seed, library_path, platform_name, os.environ)


#----- The one RNAcalibrate invocation. A non-zero EXIT is the only failure it raises on —
#      a bad FIT is the caller's to classify. `subject` names the job in that error -----#
def _run_or_raise(command, env, subject):
    try:
        return subprocess.run(command, check=True, capture_output=True,
                              text=True, env=env).stdout
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"RNAcalibrate failed for {subject}: "
            f"{' '.join(command)} exited with code {exc.returncode}. "
            f"Error output: {exc.stderr}."
        ) from exc


#----- One miRNA's reference fit at the whole-file -l. Nothing here reads another miRNA's
#      state — that independence is what lets the driver run these concurrently -----#
def _calibrate_one(query_path, header, sequence, make_command,
                   rng_seed, library_path, platform_name, pin_clock):
    mirna = query_key(header)
    command = make_command(query=str(query_path))
    derived_seed, env = _fit_env(pin_clock, rng_seed, sequence, library_path, platform_name)
    stdout = _run_or_raise(command, env, f"miRNA {mirna}")

    context = f"[miRNA {mirna}, derived_seed={derived_seed}]"
    try:
        parsed = parse_rnacalibrate_output(stdout)
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


#----- Tier 2 rungs: (std divisor, k multiplier, resulting status). The INDEX is the attempt:
#      folded into derive_seed when pinned, or forced past a clock second when not (see
#      _calibrate_anchor), so every rung draws a fresh sample at an identical command — free,
#      because near the 2.0 cutoff the outcome is near a coin flip (three identical runs gave
#      nan, nan, converged). k buys no extra bins (B caps at 500), only a wider range -----#
REPAIR_LADDER = (
    (3, 1, "fitted"),               # attempt 0 — the configured setting
    (3, 1, "fitted"),               # attempt 1 — re-draw
    (3, 1, "fitted"),               # attempt 2 — re-draw
    (2, 1, "refitted_widened"),     # attempt 3 — widen to -l L,L/2
    (3, 3, "refitted_high_k"),      # attempt 4 — raise k 3x at the configured width
)

#----- The positive counterpart to REPAIR_LADDER's statuses: "measured", never "modelled".
#      resolve_strata filters on membership here, not on excluding "rejected" — so a future
#      modelled status can never sneak into the weighted least squares unnoticed -----#
FITTED_STATUSES = frozenset(status for _, _, status in REPAIR_LADDER)


#----- The one spelling of a stratum's keys. resolve_strata merges these with {**entry, ...},
#      so a key missing from one construction site would write a ragged rnacalibrate.json -----#
def _stratum(anchor, **fields):
    return {"anchor": anchor, "length_arg": None, "k": None, "sample_size": None,
            "xi": None, "theta": None, "derived_seed": None, "attempts": 0,
            "status": "rejected", "reason": None, "residual": None, **fields}


#----- A cell with no usable fit yet; carries WHY it ended up modelled rather than measured -----#
def _rejected_stratum(anchor, reason, attempts=0):
    return _stratum(anchor, reason=reason, attempts=attempts)


#----- One (miRNA x anchor) fit with the Tier 2 ladder. Never raises on a bad FIT — Tier 3
#      turns an exhausted ladder into an extrapolated cell. A non-zero EXIT does raise:
#      that is a broken invocation, not a degenerate sample -----#
def _calibrate_anchor(query_path, header, sequence, anchor, make_command, k,
                      rng_seed, library_path, platform_name, pin_clock):
    mirna = query_key(header)

    reason = "unusable_curve"
    for attempt, (divisor, k_multiplier, status) in enumerate(REPAIR_LADDER):
        length_arg = anchor_length_arg(anchor, divisor)
        attempt_k = k * k_multiplier
        command = make_command(query=str(query_path), length_arg=length_arg, k=attempt_k)

        derived_seed, env = _fit_env(pin_clock, rng_seed, sequence, library_path,
                                     platform_name, anchor=anchor, attempt=attempt)
        if not pin_clock and attempt > 0:
            # Unpinned, so RNAcalibrate seeds from time() at 1s resolution: an identical
            # command issued immediately would collapse into the same draw. Force the clock
            # past a second boundary instead, mirroring verify_faketime's same trick.
            time.sleep(CLOCK_TICK_SLEEP_S)

        stdout = _run_or_raise(command, env, f"miRNA {mirna} at anchor {anchor}")

        context = f"[miRNA {mirna}, anchor {anchor}, attempt {attempt}, derived_seed={derived_seed}]"
        try:
            parsed = parse_rnacalibrate_output(stdout, reject_nan=False)
        except (RuntimeError, ValueError) as exc:
            # Malformed or absent output from one of many concurrent jobs — without this
            # context the error names no miRNA, anchor, or seed.
            raise RuntimeError(f"{exc} {context}.") from exc
        row = expect_single_row(parsed, f"{mirna} @ anchor {anchor}")
        reason = classify_fit(row["sample_size"], row["xi"], row["theta"])
        if reason is None:
            return _stratum(
                anchor, length_arg=length_arg["value"], k=attempt_k,
                sample_size=row["sample_size"], xi=row["xi"], theta=row["theta"],
                derived_seed=derived_seed, attempts=attempt + 1, status=status,
            )

    return _rejected_stratum(anchor, reason, attempts=len(REPAIR_LADDER))


#----- Tier 3, per miRNA: fit its OWN alpha over the sound anchors, drop 3-sigma outliers,
#      fill the rest off the curve. The residual test also catches Tier 2's retry selection
#      bias. Returns (strata, alpha, spread); alpha is None under reference_fallback -----#
def resolve_strata(entries, query_length, reference):
    sound = [entry for entry in entries if entry["status"] in FITTED_STATUSES]
    points = [FitPoint(e["anchor"], query_length, e["theta"], e["sample_size"])
              for e in sound]

    try:
        fit = fit_alpha_robust(points)
    except UnusableCurveError as exc:
        # Never halt a panel over one low-GC miRNA — fall back to its whole-file fit for
        # every cell, which is exactly today's behaviour for it. "reason" stays in its
        # fixed, machine-read vocabulary; the free-text detail is only for a human, so it
        # goes to stderr instead.
        print(f"rnacalibrate: {reference['query']} — reference_fallback ({exc})",
              file=sys.stderr)
        return [
            {**entry, "length_arg": None, "k": None,
             "sample_size": reference["sample_size"], "xi": reference["xi"],
             "theta": reference["theta"], "derived_seed": reference["derived_seed"],
             "status": "reference_fallback", "reason": entry["reason"] or "unusable_curve",
             "residual": None}
            for entry in entries
        ], None, None

    # Read off the fit's own kept set rather than re-deriving it from `dropped`: which
    # anchors survived rejection is fit_alpha_robust's answer to give, not ours.
    survivor_anchors = {point.target_length for point in fit.kept}
    survivors = [e for e in sound if e["anchor"] in survivor_anchors]
    # xi drifts ~4% across a 20x length range but is noisy per fit, so take the survivors'
    # weighted median rather than their mean or any one anchor.
    xi_modelled = weighted_median([e["xi"] for e in survivors],
                                  [e["sample_size"] for e in survivors])

    strata = []
    for entry in entries:
        if entry["anchor"] in survivor_anchors:
            point = FitPoint(entry["anchor"], query_length, entry["theta"],
                             entry["sample_size"])
            strata.append({**entry,
                           "residual": residual_ratio(point, fit.alpha, fit.log_c)})
            continue
        strata.append({
            **entry, "length_arg": None, "k": None, "sample_size": None,
            "xi": xi_modelled,
            "theta": predict_theta(fit.alpha, fit.log_c, entry["anchor"], query_length),
            "derived_seed": None, "status": "extrapolated",
            "reason": entry["reason"] or "residual_outlier", "residual": None,
        })
    return strata, fit.alpha, fit.residual_spread


#----- Top-level driver: one RNAcalibrate invocation PER miRNA (plus, when length_anchors is
#      set, one per fitted anchor). The invocations are independent, so they run on `threads`
#      cores at once; pool.map preserves input order, keeping per_query in query-FASTA order -----#
def run_rnacalibrate(query, target, output_file, k, max_target_length, randomize_targets=False,
                     max_internal_loop=None, max_bulge_loop=None, seed=None, rng_seed=None,
                     threads=1, length_anchors=None):
    executable = which_required("RNAcalibrate", "rnacalibrate")
    lengths = record_lengths(target)
    stats = summarise_lengths(lengths)
    length_arg = build_length_arg(stats)
    platform_name = sys.platform
    anchors = parse_anchors(length_anchors)

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

    common = dict(
        executable=executable, target=target, max_target_length=max_target_length,
        randomize_targets=randomize_targets, max_internal_loop=max_internal_loop,
        max_bulge_loop=max_bulge_loop, seed=seed,
    )
    # The reference fit binds k and -l; the anchor path leaves both free for Tier 2 to vary.
    # One `common` for both, so they cannot drift apart on a flag.
    make_command = partial(build_command, k=k, length_arg=length_arg, **common)
    make_anchor_command = partial(build_command, **common)

    # One table per anchor; fit_targets and occupied are views of it, never re-derived from a second copy of the counts.
    anchor_block = None
    fit_targets = ()
    occupied = ()
    if anchors:
        counts = count_records_per_anchor(lengths, anchors)
        fittable = fitted_anchors(anchors)
        anchor_block = [
            {"anchor": anchor,
             "length_arg": (anchor_length_arg(anchor)["value"]
                            if counts[anchor] and anchor in fittable else None),
             "n_records": counts[anchor],
             "status": ("empty" if not counts[anchor]
                        else "attempted" if anchor in fittable
                        else "above_ceiling")}
            for anchor in anchors
        ]
        fit_targets = tuple(r["anchor"] for r in anchor_block if r["status"] == "attempted")
        occupied = tuple(r["anchor"] for r in anchor_block if r["n_records"])

    with tempfile.TemporaryDirectory() as tmpdir:
        # Materialised before dispatch so the duplicate-ID and empty-file guards in
        # iter_query_records still raise up front, not inside a worker.
        records = list(iter_query_records(query))
        # One -q file per miRNA, written here rather than in the workers: every fit of a
        # miRNA reads the same bytes, and a file only read after dispatch cannot race.
        query_files = [_write_query_file(tmpdir, index, header, sequence)
                       for index, (header, sequence) in enumerate(records)]

        calibrate = partial(
            _calibrate_one, make_command=make_command, rng_seed=rng_seed,
            library_path=library_path, platform_name=platform_name, pin_clock=pin_clock,
        )
        calibrate_anchor = partial(
            _calibrate_anchor, make_command=make_anchor_command, k=k,
            rng_seed=rng_seed, library_path=library_path, platform_name=platform_name,
            pin_clock=pin_clock,
        )

        # One flat job list, so all 72 x 7 = 504 fits run concurrently rather than 72.
        jobs = []
        for index, (header, sequence) in enumerate(records):
            jobs.append((None, index, header, sequence))
            for anchor in fit_targets:
                jobs.append((anchor, index, header, sequence))

        def run_job(job):
            anchor, index, header, sequence = job
            if anchor is None:
                return calibrate(query_files[index], header, sequence)
            return calibrate_anchor(query_files[index], header, sequence, anchor)

        with ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
            results = list(pool.map(run_job, jobs))

    # pool.map preserves order, so results line up with jobs and per_query stays in
    # query-FASTA order.
    per_query = []
    by_index = {}
    for (anchor, index, _header, _sequence), result in zip(jobs, results):
        if anchor is None:
            by_index[index] = (result, [])
        else:
            by_index[index][1].append(result)

    for index, (header, sequence) in enumerate(records):
        entry, fitted = by_index[index]
        if anchors:
            fitted_by_anchor = {row["anchor"]: row for row in fitted}
            entries = []
            for anchor in occupied:
                row = fitted_by_anchor.get(anchor)
                entries.append(row if row is not None else _rejected_stratum(anchor, "above_ceiling"))
            strata, alpha, spread = resolve_strata(entries, len(sequence), entry)
            entry["alpha"] = alpha
            entry["alpha_residual_spread"] = spread
            entry["strata"] = strata
            if any(s["status"] != "fitted" for s in strata):
                print(
                    f"rnacalibrate: {entry['query']} — "
                    + ", ".join(f"{s['anchor']}:{s['status']}"
                                for s in strata if s["status"] != "fitted"),
                    file=sys.stderr,
                )
        per_query.append(entry)

    # "anchors" leads when there is a ladder; rnahybrid keys the shape off its presence.
    calibration = {"per_query": per_query}
    if anchor_block is not None:
        calibration = {"anchors": anchor_block, **calibration}

    ensure_parent(output_file).write_text(json.dumps(
        {
            "command": make_command(query="<per-query>"),
            "target_length_stats": stats,
            "target_length_argument": length_arg,
            "inputs": build_inputs_block(
                query=query, target=target, k=k, max_target_length=max_target_length,
                length_arg=length_arg, forced_helix=seed,
                max_internal_loop=max_internal_loop, max_bulge_loop=max_bulge_loop,
                rng_seed=rng_seed if pin_clock else None,
                length_anchors=anchors,
            ),
            "calibration": calibration,
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
        length_anchors=snakemake.params.length_anchors,
        threads=snakemake.threads,
    )

if "snakemake" in globals():
    run_from_snakemake(snakemake)