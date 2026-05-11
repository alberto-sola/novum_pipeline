#----- Handle to the config file -----#
configfile: "Config/config.yaml"

#----- Populate Snakefile variables with the config file -----#
targets             = config["targets"]
rnacalibrate_config = config.get("rnacalibrate", {}) or {}
rnahybrid_config    = config.get("rnahybrid", {}) or {}
intarna_config      = config.get("intarna", {}) or {}
results_dir         = config.get("results_dir", "Data/Results").rstrip("/")
max_target_length   = rnacalibrate_config.get("max_target_length", 50000)

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
build_plots_config  = config.get("build_plots", {}) or {}
build_plots_type    = build_plots_config.get("type")
build_plots_enabled = build_plots_type is not None
build_plots_slug    = (build_plots_type or "none").replace(",", "-").replace(" ", "")

#----- Resolve which variants to produce -----#
W_CALIBRATION  = "w_calibration"
WO_CALIBRATION = "wo_calibration"
VARIANT_LABELS = {W_CALIBRATION: "calibrated", WO_CALIBRATION: "uncalibrated"}

_MODE_TO_VARIANTS = {
    "calibrated":   [W_CALIBRATION],
    "uncalibrated": [WO_CALIBRATION],
    "both":         [W_CALIBRATION, WO_CALIBRATION],
}

def _resolve_mode(cfg):
    mode = cfg.get("mode")
    if mode is not None:
        if mode not in _MODE_TO_VARIANTS:
            raise ValueError(
                f"rnacalibrate.mode must be one of {list(_MODE_TO_VARIANTS)}, got {mode!r}"
            )
        return mode
    # Deprecated fallback: `enabled: true|false` maps to a single-variant run.
    enabled = cfg.get("enabled", False)
    return "calibrated" if enabled else "uncalibrated"

variants = _MODE_TO_VARIANTS[_resolve_mode(rnacalibrate_config)]

#----- Resolve which IntaRNA accessibility variants to produce -----#
W_ACCESSIBILITY  = "w_accessibility"
WO_ACCESSIBILITY = "wo_accessibility"

_ACC_MODE_TO_VARIANTS = {
    "on":   [W_ACCESSIBILITY],
    "off":  [WO_ACCESSIBILITY],
    "both": [W_ACCESSIBILITY, WO_ACCESSIBILITY],
}

def _resolve_acc_mode(cfg):
    mode = cfg.get("acc_mode", "both")
    # YAML 1.1 parses bare `on`/`off` as Python True/False — coerce back to strings.
    if isinstance(mode, bool):
        mode = "on" if mode else "off"
    if mode not in _ACC_MODE_TO_VARIANTS:
        raise ValueError(
            f"intarna.acc_mode must be one of {list(_ACC_MODE_TO_VARIANTS)}, got {mode!r}"
        )
    return mode

intarna_variants = _ACC_MODE_TO_VARIANTS[_resolve_acc_mode(intarna_config)]

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
        *(expand(
            f"{SAMPLE_DIR}/plots_" + build_plots_slug + ".pdf",
            sample=targets.keys(),
        ) if build_plots_enabled else [])

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
        u=rnahybrid_config.get("u"),
        v=rnahybrid_config.get("v"),
        seed=rnahybrid_config.get("seed")
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
        hits=rnahybrid_config.get("hits"),
        u=rnahybrid_config.get("u"),
        v=rnahybrid_config.get("v"),
        energy=rnahybrid_config.get("energy"),
        pvalue=rnahybrid_config.get("pvalue"),
        seed=rnahybrid_config.get("seed"),
        distribution=rnahybrid_config.get("distribution"),
        max_target_length=max_target_length
    threads:
        int(rnahybrid_config.get("threads", 1))
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
        int(intarna_config.get("threads", 1))
    params:
        intarna=intarna_config
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

#----- Build configurable ggplot2 plots from the annotated tables -----#
rule build_plots:
    input:
        annotated=expand(
            f"{results_dir}" + "/{{sample}}/{variant}/rnahybrid_annotated.csv",
            variant=variants,
        )
    output:
        pdf=f"{SAMPLE_DIR}/plots_" + build_plots_slug + ".pdf"
    conda:
        "Workflow/Envs/plots.yaml"
    params:
        type=build_plots_type,
        basesize=build_plots_config.get("basesize", 12),
        pvalue=build_plots_config.get("pvalue", []),
        locus=build_plots_config.get("locus", []),
        gene=build_plots_config.get("gene", []),
        protein=build_plots_config.get("protein", []),
        per_mirna_top_n=build_plots_config.get("per_mirna_top_n", 12),
        variant_labels=[VARIANT_LABELS[v] for v in variants]
    script:
        "Workflow/Scripts/build_plots.R"
