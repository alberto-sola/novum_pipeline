#----- Handle to the config file -----#
configfile: "Config/config.yaml"

#----- Populate Snakefile variables with the config file -----#
targets             = config["targets"]
rnacalibrate_config = config.get("rnacalibrate", {}) or {}
rnahybrid_config    = config.get("rnahybrid", {}) or {}
results_dir         = config.get("results_dir", "Data/Results").rstrip("/")
max_target_length   = rnacalibrate_config.get("max_target_length", 50000)

#----- Per-sample query lookup: `queries:` mapping (keys must match `targets:`),
#      with single `query: <path>` accepted as a legacy broadcast across all targets.
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


wildcard_constraints:
    variant = f"{W_CALIBRATION}|{WO_CALIBRATION}"


def calibration_input(wc):
    if wc.variant == W_CALIBRATION:
        return f"{results_dir}/{wc.sample}/{W_CALIBRATION}/rnacalibrate.json"
    return []


#----- output finale ‒ pipeline conclusion -----#
rule all:
    input:
        expand(
            f"{results_dir}" + "/{sample}/{variant}/rnahybrid_enhanced.txt",
            sample=targets.keys(),
            variant=variants,
        ),
        *(expand(
            f"{results_dir}" + "/{sample}/plots_" + build_plots_slug + ".pdf",
            sample=targets.keys(),
        ) if build_plots_enabled else [])

#----- Dynamically calibrates the statistics based on the target sequence -----#
rule rnacalibrate:
    input:
        query=lambda wc: queries[wc.sample],
        target=lambda wc: targets[wc.sample]
    output:
        calibration=f"{results_dir}" + "/{sample}/" + W_CALIBRATION + "/rnacalibrate.json"
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
    input:
        query=lambda wc: queries[wc.sample],
        target=lambda wc: targets[wc.sample],
        calibration=calibration_input
    output:
        compact=f"{results_dir}" + "/{sample}/{variant}/rnahybrid_output.tsv"
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
    input:
        compact=f"{results_dir}" + "/{sample}/{variant}/rnahybrid_output.tsv"
    output:
        tidy=f"{results_dir}" + "/{sample}/{variant}/tidy_output.csv"
    conda:
        "Workflow/Envs/postprocess.yaml"
    script:
        "Workflow/Scripts/tidy_rnahybrid.py"

#----- Parse metadata from the target file (genome) and merge them with the tidied output -----#
rule annotate_rnahybrid:
    input:
        tidy=f"{results_dir}" + "/{sample}/{variant}/tidy_output.csv",
        target=lambda wc: targets[wc.sample]
    output:
        annotated=f"{results_dir}" + "/{sample}/{variant}/rnahybrid_annotated.csv"
    conda:
        "Workflow/Envs/postprocess.yaml"
    script:
        "Workflow/Scripts/annotate_rnahybrid.py"

#----- Visualize RNAhybrid alignments + all the data -----#
rule enhance_rnahybrid:
    input:
        annotated=f"{results_dir}" + "/{sample}/{variant}/rnahybrid_annotated.csv"
    output:
        enhanced=f"{results_dir}" + "/{sample}/{variant}/rnahybrid_enhanced.txt"
    conda:
        "Workflow/Envs/postprocess.yaml"
    script:
        "Workflow/Scripts/enhance_rnahybrid.py"

#----- Build configurable ggplot2 plots from the annotated tables -----#
rule build_plots:
    input:
        annotated=expand(
            f"{results_dir}" + "/{{sample}}/{variant}/rnahybrid_annotated.csv",
            variant=variants,
        )
    output:
        pdf=f"{results_dir}" + "/{sample}/plots_" + build_plots_slug + ".pdf"
    conda:
        "Workflow/Envs/plots.yaml"
    params:
        type=build_plots_type,
        basesize=build_plots_config.get("basesize", 12),
        pvalue=build_plots_config.get("pvalue", 0.01),
        locus=build_plots_config.get("locus", []),
        gene=build_plots_config.get("gene", []),
        protein=build_plots_config.get("protein", []),
        variant_labels=[VARIANT_LABELS[v] for v in variants]
    script:
        "Workflow/Scripts/build_plots.R"