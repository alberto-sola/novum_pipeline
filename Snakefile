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
shared_max_total_energy    = config.get("max_total_energy")
shared_seed                = config.get("seed")

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


#----- Per-variant accessibility routing: resolve the {variant} wildcard to IntaRNA's
#      --acc mode (N=none, C=constrained) here, mirroring calibration_input, so the
#      script consumes a resolved knob instead of re-declaring the variant literal -----#
def intarna_acc_mode(wc):
    return "N" if wc.variant == WO_ACCESSIBILITY else "C"


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
        max_total_energy=shared_max_total_energy,
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
        max_suboptimal_hits = shared_max_suboptimal_hits,
        max_total_energy    = shared_max_total_energy,
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
