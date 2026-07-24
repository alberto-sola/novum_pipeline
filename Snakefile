import os
import sys

#----- Handle to the config file -----#
configfile: "Config/config.yaml"

#----- Populate Snakefile variables with the config file -----#
targets             = config["targets"]
rnacalibrate_config = config.get("rnacalibrate", {}) or {}
rnahybrid_config    = config.get("rnahybrid", {}) or {}
intarna_config      = config.get("intarna", {}) or {}
results_dir         = config.get("results_dir", "Data/Results").rstrip("/")
max_target_length   = rnacalibrate_config.get("max_target_length", 50000)
shared_threads             = int(config.get("threads", 1))
shared_max_suboptimal_hits = config.get("max_suboptimal_hits")
shared_seed                = config.get("seed")
rnahybrid_max_hybrid_energy = rnahybrid_config.get("max_hybrid_energy")
intarna_max_hybrid_energy   = intarna_config.get("max_hybrid_energy")

SAMPLE_DIR         = f"{results_dir}" + "/{sample}"
SAMPLE_VARIANT_DIR = SAMPLE_DIR + "/{variant}"

#----- Per-sample query lookup: `queries:` mapping (keys must match `targets:`),
#      with single `query: <path>` accepted as a legacy broadcast across all targets -----#
def _resolve_queries(cfg, target_keys):
    queries_cfg = cfg.get("queries")
    if queries_cfg is not None:
        missing = set(target_keys) - set(queries_cfg)
        extra   = set(queries_cfg) - set(target_keys)
        if missing or extra:
            raise ValueError(
                "queries/targets keys must match. "
                f"Missing in queries: {sorted(missing)}; extra in queries: {sorted(extra)}."
            )
        return dict(queries_cfg)
    legacy = cfg.get("query")
    if legacy is None:
        raise ValueError("Config must define either `queries:` (mapping) or `query:` (single path).")
    return {sample: legacy for sample in target_keys}

queries = _resolve_queries(config, targets.keys())

#----- Optional plotting configuration -----#
plots_config  = config.get("plots", {}) or {}
plots_type    = plots_config.get("type")
plots_enabled = plots_type is not None
plots_slug    = (plots_type or "none").replace(",", "-").replace(" ", "")

#----- Variant axes: each arm maps a config "mode" to the output-tree literals it expands into -----#
W_CALIBRATION    = "w_calibration"
WO_CALIBRATION   = "wo_calibration"
W_ACCESSIBILITY  = "w_accessibility"
WO_ACCESSIBILITY = "wo_accessibility"
VARIANT_LABELS = {W_CALIBRATION: "calibrated", WO_CALIBRATION: "uncalibrated"}

#----- on→with-variant / off→without / both→both; one builder so the shape can't drift -----#
def _variant_map(with_literal, without_literal):
    return {"on": [with_literal], "off": [without_literal], "both": [with_literal, without_literal]}

_CALIBRATION_VARIANTS   = _variant_map(W_CALIBRATION,   WO_CALIBRATION)
_ACCESSIBILITY_VARIANTS = _variant_map(W_ACCESSIBILITY, WO_ACCESSIBILITY)

#----- Shared validate+lookup: turn a resolved mode into its variant list, or fail listing the allowed set -----#
def variants_for_mode(mode, mode_to_variants, label):
    if mode not in mode_to_variants:
        raise ValueError(
            f"{label} must be one of {list(mode_to_variants)}, got {mode!r}"
        )
    return mode_to_variants[mode]

#----- Shared variant resolver: read `key` from cfg (default "both"); YAML 1.1 parses bare on/off as bools, so coerce back -----#
def _variant_mode(cfg, key):
    mode = cfg.get(key, "both")
    if isinstance(mode, bool):
        mode = "on" if mode else "off"
    return mode

variants = variants_for_mode(
    _variant_mode(rnacalibrate_config, "calibration_variant"),
    _CALIBRATION_VARIANTS, "rnacalibrate.calibration_variant",
)
intarna_variants = variants_for_mode(
    _variant_mode(intarna_config, "accessibility_variant"),
    _ACCESSIBILITY_VARIANTS, "intarna.accessibility_variant",
)

#----- IntaRNA config validation at DAG-build time, so a bad config fails before the
#      hour-long RNAhybrid arm starts. _intarna_config is stdlib-only by contract: it is
#      imported here into the Snakemake DRIVER env, which ships none of the arms'
#      scientific deps. Never import a script that ends in run_from_snakemake() instead. -----#
sys.path.insert(0, os.path.join(workflow.basedir, "Workflow", "Scripts"))
from _intarna_config import (
    INTARNA_DEFAULT_ACCESSIBILITY_SEARCH_DEPTH,
    opt,
    reject_removed_keys,
    validate_intarna_config,
)

# Fatal before anything else: a removed key left in a config is inert but looks live.
reject_removed_keys(config)

# One validation call, one emit site: every fatal and every warning lives in the
# driver-importable module, where each is unit-tested.
for _warning in validate_intarna_config(
    intarna_config, shared_seed,
    rnahybrid_max_hybrid_energy=rnahybrid_max_hybrid_energy,
    max_suboptimal_hits=shared_max_suboptimal_hits,
    accessibility_on=W_ACCESSIBILITY in intarna_variants,
):
    logger.warning(_warning)

# opt(), not .get(default): an explicit `accessibility_search_depth:` null returns None
# from .get and would blow up inside the max() in intarna_outnumber.
intarna_search_depth = opt(intarna_config, "accessibility_search_depth",
                           INTARNA_DEFAULT_ACCESSIBILITY_SEARCH_DEPTH)


RNAHYBRID_VARIANT_RE = f"{W_CALIBRATION}|{WO_CALIBRATION}"
INTARNA_VARIANT_RE   = f"{W_ACCESSIBILITY}|{WO_ACCESSIBILITY}"


#----- Per-variant calibration routing: the {variant} wildcard is constrained to
#      the two calibration literals; calibration_input(wc) supplies the JSON only
#      for w_calibration so a single rnahybrid rule serves both branches -----#
wildcard_constraints:
    variant = "|".join([W_CALIBRATION, WO_CALIBRATION, W_ACCESSIBILITY, WO_ACCESSIBILITY])


def calibration_input(wc):
    if wc.variant == W_CALIBRATION:
        return f"{results_dir}/{wc.sample}/{W_CALIBRATION}/rnacalibrate.json"
    return []


#----- The one predicate all three IntaRNA resolvers below branch on. Under the per-rule
#      INTARNA_VARIANT_RE constraint the variant is one of exactly two literals, so
#      "accessibility on" and "!= wo_accessibility" coincide — spell it once -----#
def _accessibility_on(wc):
    return wc.variant == W_ACCESSIBILITY


#----- Per-variant accessibility routing: resolve the {variant} wildcard to IntaRNA's
#      --acc mode (N=none, C=constrained) here, mirroring calibration_input, so the
#      script consumes a resolved knob instead of re-declaring the variant literal -----#
def intarna_acc_mode(wc):
    return "C" if _accessibility_on(wc) else "N"


#----- --outMaxE filters TOTAL E, which equals E_hybrid only under acc=N. Emit it as a
#      tool-level pre-filter exactly where that holds, and defer to tidy_intarna's
#      E_hybrid gate under acc=C, where the same number would be a ~2x stricter bar -----#
def intarna_outmaxe(wc):
    return None if _accessibility_on(wc) else intarna_max_hybrid_energy


#----- acc=C ranks by total E while the gate is on E_hybrid, so the strongest duplex can
#      sit outside a top-1 list. Report deeper there and let tidy_intarna select. A floor,
#      never a ceiling: a null shared knob means "all hits" and must not narrow to the
#      depth — which is also why warning (e) tells the user only the cap can bound it -----#
def intarna_outnumber(wc):
    if not _accessibility_on(wc) or shared_max_suboptimal_hits is None:
        return shared_max_suboptimal_hits
    return max(intarna_search_depth, shared_max_suboptimal_hits)


#----- Consensus arm-variant selection: prefer the w_* variant of each arm,
#      fall back to wo_* when only that ran; resolved once at DAG-build time -----#
def _preferred_variant(available, preferred, fallback):
    return preferred if preferred in available else fallback

def _consensus_rnahybrid_annotated(wc):
    variant = _preferred_variant(variants, W_CALIBRATION, WO_CALIBRATION)
    return f"{results_dir}/{wc.sample}/{variant}/rnahybrid_annotated.csv"

def _consensus_intarna_annotated(wc):
    variant = _preferred_variant(intarna_variants, W_ACCESSIBILITY, WO_ACCESSIBILITY)
    return f"{results_dir}/{wc.sample}/{variant}/intarna_annotated.csv"


#----- output finale ‒ pipeline conclusion -----#
rule all:
    input:
        expand(
            f"{SAMPLE_VARIANT_DIR}/rnahybrid_enhanced.txt",
            sample=targets.keys(),
            variant=variants,
        ),
        expand(
            f"{SAMPLE_VARIANT_DIR}/intarna_enhanced.txt",
            sample=targets.keys(),
            variant=intarna_variants,
        ),
        expand(
            f"{SAMPLE_DIR}/consensus/consensus_enhanced.txt",
            sample=targets.keys(),
        ),
        *(expand(
            f"{SAMPLE_DIR}/plots_{plots_slug}.pdf",
            sample=targets.keys(),
        ) if plots_enabled else [])

#----- Dynamically calibrates the statistics based on the target sequence -----#
rule rnacalibrate:
    input:
        query=lambda wc: queries[wc.sample],
        target=lambda wc: targets[wc.sample]
    output:
        calibration=f"{SAMPLE_DIR}/" + W_CALIBRATION + "/rnacalibrate.json"
    conda:
        "Workflow/Envs/rnahybrid.yaml"
    params:
        k=rnacalibrate_config.get("k", 10000),
        max_target_length=max_target_length,
        randomize_targets=rnacalibrate_config.get("randomize_targets", False),
        max_internal_loop=rnahybrid_config.get("max_internal_loop"),
        max_bulge_loop=rnahybrid_config.get("max_bulge_loop"),
        seed=shared_seed
    script:
        "Workflow/Scripts/rnacalibrate.py"

#----- Run RNAhybrid on multiple cores -----#
rule rnahybrid:
    wildcard_constraints:
        variant = RNAHYBRID_VARIANT_RE
    input:
        query=lambda wc: queries[wc.sample],
        target=lambda wc: targets[wc.sample],
        calibration=calibration_input
    output:
        compact=f"{SAMPLE_VARIANT_DIR}/rnahybrid_output.tsv"
    conda:
        "Workflow/Envs/rnahybrid.yaml"
    params:
        species=rnahybrid_config.get("species", "3utr_human"),
        max_suboptimal_hits=shared_max_suboptimal_hits,
        max_internal_loop=rnahybrid_config.get("max_internal_loop"),
        max_bulge_loop=rnahybrid_config.get("max_bulge_loop"),
        max_hybrid_energy=rnahybrid_max_hybrid_energy,
        pvalue_threshold=rnahybrid_config.get("pvalue_threshold"),
        seed=shared_seed,
        distribution=rnahybrid_config.get("distribution"),
        max_target_length=max_target_length
    threads:
        shared_threads
    script:
        "Workflow/Scripts/rnahybrid.py"

#----- Filter, select, arrange the RNAhybrid's output -----#
rule tidy_rnahybrid:
    wildcard_constraints:
        variant = RNAHYBRID_VARIANT_RE
    input:
        compact=f"{SAMPLE_VARIANT_DIR}/rnahybrid_output.tsv"
    output:
        tidy=f"{SAMPLE_VARIANT_DIR}/tidy_output.csv"
    conda:
        "Workflow/Envs/postprocess.yaml"
    script:
        "Workflow/Scripts/tidy_rnahybrid.py"

#----- Parse metadata from the target file (genome) and merge them with the tidied output -----#
rule annotate_rnahybrid:
    wildcard_constraints:
        variant = RNAHYBRID_VARIANT_RE
    input:
        tidy=f"{SAMPLE_VARIANT_DIR}/tidy_output.csv",
        target=lambda wc: targets[wc.sample]
    output:
        annotated=f"{SAMPLE_VARIANT_DIR}/rnahybrid_annotated.csv"
    conda:
        "Workflow/Envs/postprocess.yaml"
    params:
        insert_after="P_value"
    script:
        "Workflow/Scripts/annotate.py"

#----- Visualize RNAhybrid alignments + all the data -----#
rule enhance_rnahybrid:
    wildcard_constraints:
        variant = RNAHYBRID_VARIANT_RE
    input:
        annotated=f"{SAMPLE_VARIANT_DIR}/rnahybrid_annotated.csv"
    output:
        enhanced=f"{SAMPLE_VARIANT_DIR}/rnahybrid_enhanced.txt"
    conda:
        "Workflow/Envs/postprocess.yaml"
    script:
        "Workflow/Scripts/enhance_rnahybrid.py"

#----- Run IntaRNA on multiple cores (native --threads, no external chunking) -----#
rule intarna:
    wildcard_constraints:
        variant = INTARNA_VARIANT_RE
    input:
        query=lambda wc: queries[wc.sample],
        target=lambda wc: targets[wc.sample]
    output:
        csv=f"{SAMPLE_VARIANT_DIR}/intarna_output.csv"
    conda:
        "Workflow/Envs/intarna.yaml"
    threads:
        shared_threads
    params:
        intarna             = intarna_config,
        acc                 = intarna_acc_mode,
        max_suboptimal_hits = intarna_outnumber,
        out_max_energy      = intarna_outmaxe,
        seed                = shared_seed
    script:
        "Workflow/Scripts/intarna.py"

#----- Filter, rename, and normalize the IntaRNA CSV output -----#
rule tidy_intarna:
    wildcard_constraints:
        variant = INTARNA_VARIANT_RE
    input:
        csv=f"{SAMPLE_VARIANT_DIR}/intarna_output.csv",
        target=lambda wc: targets[wc.sample]
    output:
        tidy=f"{SAMPLE_VARIANT_DIR}/intarna_tidy.csv"
    conda:
        "Workflow/Envs/intarna.yaml"
    params:
        # Authoritative gate on both variants: --outMaxE cannot express an E_hybrid bound
        # under acc=C. Idempotent on wo_accessibility, where the tool already applied it.
        max_hybrid_energy   = intarna_max_hybrid_energy,
        max_suboptimal_hits = shared_max_suboptimal_hits
    script:
        "Workflow/Scripts/tidy_intarna.py"

#----- Parse metadata from the target FASTA and merge with the IntaRNA tidy output -----#
rule annotate_intarna:
    wildcard_constraints:
        variant = INTARNA_VARIANT_RE
    input:
        tidy=f"{SAMPLE_VARIANT_DIR}/intarna_tidy.csv",
        target=lambda wc: targets[wc.sample]
    output:
        annotated=f"{SAMPLE_VARIANT_DIR}/intarna_annotated.csv"
    conda:
        "Workflow/Envs/postprocess.yaml"
    params:
        insert_after="E"
    script:
        "Workflow/Scripts/annotate.py"

#----- Visualize IntaRNA alignments and energy data as a per-record human-readable report -----#
rule enhance_intarna:
    wildcard_constraints:
        variant = INTARNA_VARIANT_RE
    input:
        annotated=f"{SAMPLE_VARIANT_DIR}/intarna_annotated.csv"
    output:
        enhanced=f"{SAMPLE_VARIANT_DIR}/intarna_enhanced.txt"
    conda:
        "Workflow/Envs/postprocess.yaml"
    script:
        "Workflow/Scripts/enhance_intarna.py"

#----- Intersect the two arms' annotated tables into the consensus (both-tools) set -----#
rule intersect:
    input:
        rnahybrid=_consensus_rnahybrid_annotated,
        intarna=_consensus_intarna_annotated
    output:
        annotated=f"{SAMPLE_DIR}/consensus/consensus_annotated.csv"
    conda:
        "Workflow/Envs/postprocess.yaml"
    script:
        "Workflow/Scripts/intersect.py"

#----- Render the consensus set as a per-record report: merged metadata + both duplexes -----#
rule enhance_consensus:
    input:
        annotated=f"{SAMPLE_DIR}/consensus/consensus_annotated.csv"
    output:
        enhanced=f"{SAMPLE_DIR}/consensus/consensus_enhanced.txt"
    conda:
        "Workflow/Envs/postprocess.yaml"
    script:
        "Workflow/Scripts/enhance_consensus.py"

#----- Build configurable ggplot2 plots from the annotated tables -----#
rule build_plots:
    input:
        annotated=expand(
            f"{results_dir}" + "/{{sample}}/{variant}/rnahybrid_annotated.csv",
            variant=variants,
        )
    output:
        pdf=f"{SAMPLE_DIR}/plots_{plots_slug}.pdf"
    conda:
        "Workflow/Envs/plots.yaml"
    params:
        type=plots_type,
        basesize=plots_config.get("basesize", 12),
        pvalue_threshold=plots_config.get("pvalue_threshold", []),
        locus=plots_config.get("locus", []),
        gene=plots_config.get("gene", []),
        protein=plots_config.get("protein", []),
        per_mirna_top_n=plots_config.get("per_mirna_top_n", 12),
        variant_labels=[VARIANT_LABELS[v] for v in variants]
    script:
        "Workflow/Scripts/build_plots.R"
