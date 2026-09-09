import os
import sys

#----- Handle to the config file -----#
configfile: "Config/config.yaml"

#----- IntaRNA config schema + validation, imported before anything reads the config so a bad
#      config fails at DAG-build time, before the hour-long RNAhybrid arm starts.
#      _intarna_config is stdlib-only by contract: it is imported here into the Snakemake
#      DRIVER env, which ships none of the arms' scientific deps. Never import a script that
#      ends in run_from_snakemake() instead -----#
sys.path.insert(0, os.path.join(workflow.basedir, "Workflow", "Scripts"))
from _intarna_config import (
    ACCESSIBILITY_PU_FLOOR_KEYS,
    INTARNA_DEFAULT_ACCESSIBILITY_SEARCH_DEPTH,
    opt,
    reject_removed_keys,
    tool_config,
    validate_intarna_config,
)

#----- Populate Snakefile variables with the config file. opt() wherever a default exists:
#      an explicit YAML null reaches .get as None and would flow on as the string "None".
#      Plain .get() only where null is itself meaningful — "all hits", --noSeed, ungated -----#
targets             = config["targets"]
rnacalibrate_config = config.get("rnacalibrate", {}) or {}
rnahybrid_config    = config.get("rnahybrid", {}) or {}
intarna_config      = config.get("intarna", {}) or {}
results_dir         = opt(config, "results_dir", "Data/Results").rstrip("/")
max_target_length   = opt(rnacalibrate_config, "max_target_length", 50000)
shared_threads             = int(opt(config, "threads", 1))
shared_max_suboptimal_hits = config.get("max_suboptimal_hits")
shared_seed                = config.get("seed")
rnahybrid_max_hybrid_energy = rnahybrid_config.get("max_hybrid_energy")
intarna_max_hybrid_energy   = intarna_config.get("max_hybrid_energy")

SAMPLE_DIR         = results_dir + "/{sample}"
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


#----- The miRNA FASTA for this sample, named once instead of a `lambda wc:` respelled at
#      every input slot below -----#
def sample_query(wc):
    return queries[wc.sample]


#----- The CDS FASTA for this sample -----#
def sample_target(wc):
    return targets[wc.sample]

#----- Optional plotting configuration -----#
plots_config  = config.get("plots", {}) or {}
plots_type    = plots_config.get("type")
plots_enabled = plots_type is not None
plots_slug    = (plots_type or "none").replace(",", "-").replace(" ", "")

#----- Variant axes: each arm has one boolean axis, encoded in the output tree as a
#      (with, without) pair of sibling directory names. Naming each pair once is what keeps
#      the wildcard regexes and the consensus fallback below from drifting apart -----#
W_CALIBRATION    = "w_calibration"
WO_CALIBRATION   = "wo_calibration"
W_ACCESSIBILITY  = "w_accessibility"
WO_ACCESSIBILITY = "wo_accessibility"
CALIBRATION_LITERALS   = (W_CALIBRATION,   WO_CALIBRATION)
ACCESSIBILITY_LITERALS = (W_ACCESSIBILITY, WO_ACCESSIBILITY)
VARIANT_LABELS = {W_CALIBRATION: "calibrated", WO_CALIBRATION: "uncalibrated"}

#----- The one resolver both arms go through: read `key` from cfg (default "both"), coerce
#      YAML 1.1's bare on/off back from bool, and expand on→with / off→without / both→both
#      into the arm's directory literals, failing with the allowed set on anything else -----#
def variants_for(cfg, key, label, literals):
    with_literal, without_literal = literals
    mode = cfg.get(key, "both")
    if isinstance(mode, bool):
        mode = "on" if mode else "off"
    modes = {"on": [with_literal], "off": [without_literal],
             "both": [with_literal, without_literal]}
    if mode not in modes:
        raise ValueError(f"{label} must be one of {list(modes)}, got {mode!r}")
    return modes[mode]

variants = variants_for(rnacalibrate_config, "calibration_variant",
                        "rnacalibrate.calibration_variant", CALIBRATION_LITERALS)
intarna_variants = variants_for(intarna_config, "accessibility_variant",
                                "intarna.accessibility_variant", ACCESSIBILITY_LITERALS)

#----- Fatal before anything else: a removed key left in a config is inert but looks live -----#
reject_removed_keys(config)

#----- One validation call, one emit site: every fatal and every warning lives in the
#      driver-importable module, where each is unit-tested -----#
for _warning in validate_intarna_config(
    intarna_config, shared_seed,
    rnahybrid_max_hybrid_energy=rnahybrid_max_hybrid_energy,
    max_suboptimal_hits=shared_max_suboptimal_hits,
    accessibility_on=W_ACCESSIBILITY in intarna_variants,
):
    logger.warning(_warning)

#----- A floor on --outNumber for the accessibility arm; opt() because an explicit null would
#      otherwise reach the max() in intarna_outnumber -----#
intarna_search_depth = opt(intarna_config, "accessibility_search_depth",
                           INTARNA_DEFAULT_ACCESSIBILITY_SEARCH_DEPTH)


RNAHYBRID_VARIANT_RE = "|".join(CALIBRATION_LITERALS)
INTARNA_VARIANT_RE   = "|".join(ACCESSIBILITY_LITERALS)

# Spelled once: rule rnacalibrate's output and calibration_input's return must agree, and a
# divergence would surface as a missing-input error rather than as the wiring bug it is.
RNACALIBRATE_JSON = SAMPLE_DIR + "/" + W_CALIBRATION + "/rnacalibrate.json"


#----- Every {variant} wildcard is one of the four arm literals. A safety net for rules added
#      later: each rule below already narrows to its own arm's two, so nothing currently
#      resolves through this union -----#
wildcard_constraints:
    variant = f"{RNAHYBRID_VARIANT_RE}|{INTARNA_VARIANT_RE}"


#----- Per-variant calibration routing: supplies the JSON only for w_calibration, so a
#      single rnahybrid rule serves both branches. Preserve the empty-list-vs-path
#      contract — rnahybrid.py coerces `[]` to None to pick the broadcast branch -----#
def calibration_input(wc):
    if wc.variant == W_CALIBRATION:
        return RNACALIBRATE_JSON.format(sample=wc.sample)
    return []


#----- The one predicate every IntaRNA resolver below branches on. Under the per-rule
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


#----- Per-side Pu floors: the mirror of intarna_outmaxe. That one withholds under acc=C;
#      these withhold under acc=N, which computes no Pu at all — `accessibility_variant:
#      both` must not carry a floor into the wo_ arm and empty it. One factory over the
#      shared key table, so the two sides cannot drift from each other or from warning (h) -----#
def _pu_floor_resolver(key):
    return lambda wc: intarna_config.get(key) if _accessibility_on(wc) else None

intarna_min_pu_target, intarna_min_pu_query = (
    _pu_floor_resolver(key) for key in ACCESSIBILITY_PU_FLOOR_KEYS
)


#----- Consensus arm-variant selection: prefer the w_* variant of each arm, fall back to
#      wo_* when only that ran. Constants, not input functions: the choice reads the config
#      and never the wildcards, so it genuinely resolves once at DAG-build time -----#
def _preferred_variant(available, literals):
    with_literal, without_literal = literals
    return with_literal if with_literal in available else without_literal

CONSENSUS_RNAHYBRID_VARIANT = _preferred_variant(variants, CALIBRATION_LITERALS)
CONSENSUS_INTARNA_VARIANT   = _preferred_variant(intarna_variants, ACCESSIBILITY_LITERALS)


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
        query=sample_query,
        target=sample_target
    output:
        calibration=RNACALIBRATE_JSON
    conda:
        "Workflow/Envs/rnahybrid.yaml"
    threads:
        # The per-miRNA invocations are independent (each is position 1 in its own RNG
        # stream), so the split runs concurrently — this rule was the pipeline's single
        # largest serial cost.
        shared_threads
    params:
        k=rnacalibrate_config.get("k", 10000),
        max_target_length=max_target_length,
        randomize_targets=rnacalibrate_config.get("randomize_targets", False),
        rng_seed=rnacalibrate_config.get("rng_seed"),
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
        query=sample_query,
        target=sample_target,
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
        target=sample_target
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
        query=sample_query,
        target=sample_target
    output:
        csv=f"{SAMPLE_VARIANT_DIR}/intarna_output.csv"
    conda:
        "Workflow/Envs/intarna.yaml"
    threads:
        shared_threads
    params:
        # tool_config(), not the whole block: params are a rerun trigger, so handing this rule
        # the downstream-only keys (the tidy gate, the per-side Pu floors) would make retuning
        # a cheap tidy_intarna knob invalidate every intarna_output.csv.
        intarna             = tool_config(intarna_config),
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
        target=sample_target
    output:
        tidy=f"{SAMPLE_VARIANT_DIR}/intarna_tidy.csv"
    resources:
        mem_mb = 1200
    conda:
        # postprocess, not intarna: this rule needs pandas and never invokes the binary,
        # and pinning it here is what let intarna.yaml drop pandas entirely.
        "Workflow/Envs/postprocess.yaml"
    params:
        # Authoritative gate on both variants: --outMaxE cannot express an E_hybrid bound
        # under acc=C. Idempotent on wo_accessibility, where the tool already applied it.
        max_hybrid_energy   = intarna_max_hybrid_energy,
        max_suboptimal_hits = shared_max_suboptimal_hits,
        # Withheld under acc=N by the resolvers above; applied before the cap, so they
        # select which site represents a pair rather than only rejecting the chosen one.
        min_target_unpaired_probability = intarna_min_pu_target,
        min_query_unpaired_probability  = intarna_min_pu_query
    script:
        "Workflow/Scripts/tidy_intarna.py"

#----- Parse metadata from the target FASTA and merge with the IntaRNA tidy output -----#
rule annotate_intarna:
    wildcard_constraints:
        variant = INTARNA_VARIANT_RE
    input:
        tidy=f"{SAMPLE_VARIANT_DIR}/intarna_tidy.csv",
        target=sample_target
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
        rnahybrid=f"{SAMPLE_DIR}/{CONSENSUS_RNAHYBRID_VARIANT}/rnahybrid_annotated.csv",
        intarna=f"{SAMPLE_DIR}/{CONSENSUS_INTARNA_VARIANT}/intarna_annotated.csv"
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
            f"{SAMPLE_VARIANT_DIR}/rnahybrid_annotated.csv",
            variant=variants,
            allow_missing=True,
        )
    output:
        pdf=f"{SAMPLE_DIR}/plots_{plots_slug}.pdf"
    conda:
        "Workflow/Envs/plots.yaml"
    params:
        # opt() because the shipped config writes `locus:`/`gene:`/`protein:` as bare nulls.
        type=plots_type,
        basesize=opt(plots_config, "basesize", 12),
        pvalue_threshold=opt(plots_config, "pvalue_threshold", []),
        locus=opt(plots_config, "locus", []),
        gene=opt(plots_config, "gene", []),
        protein=opt(plots_config, "protein", []),
        per_mirna_top_n=opt(plots_config, "per_mirna_top_n", 12),
        variant_labels=[VARIANT_LABELS[v] for v in variants]
    script:
        "Workflow/Scripts/build_plots.R"
