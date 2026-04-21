#----- Handle to the config file -----#
configfile: "Config/config.yaml"

#----- Populate Snakefile variables with the config file -----#
query = config["query"]
targets = config["targets"]
rnacalibrate_config = config["rnacalibrate"]
rnahybrid_config = config["rnahybrid"]
results_dir = config.get("results_dir", "Data/Results").rstrip("/")

#----- Resolve which variants to produce -----#
_MODE_TO_VARIANTS = {
    "calibrated":   ["w_calibration"],
    "uncalibrated": ["wo_calibration"],
    "both":         ["w_calibration", "wo_calibration"],
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
    variant = r"w_calibration|wo_calibration"


def calibration_input(wc):
    if wc.variant == "w_calibration":
        return f"{results_dir}/{wc.sample}/w_calibration/rnacalibrate.json"
    return []


#----- output finale ‒ pipeline conclusion -----#
rule all:
    input:
        expand(
            f"{results_dir}" + "/{sample}/{variant}/rnahybrid_enhanced.txt",
            sample=targets.keys(),
            variant=variants,
        )

#----- Dynamically calibrates the statistics based on the target sequence -----#
rule rnacalibrate:
    input:
        query=query,
        target=lambda wc: targets[wc.sample]
    output:
        calibration=f"{results_dir}" + "/{sample}/w_calibration/rnacalibrate.json"
    conda:
        "Workflow/Envs/rnahybrid.yaml"
    params:
        k=rnacalibrate_config.get("k", 5000),
        max_target_length=rnacalibrate_config.get("max_target_length", 50000),
        randomize_targets=rnacalibrate_config.get("randomize_targets", rnacalibrate_config.get("use_target_distribution", False)),
        u=rnahybrid_config.get("u"),
        v=rnahybrid_config.get("v"),
        seed=rnahybrid_config.get("seed")
    script:
        "Workflow/Scripts/rnacalibrate.py"

#----- Run RNAhybrid on multiple cores -----#
rule rnahybrid:
    input:
        query=query,
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
        distribution_file=lambda wc, input: input.calibration if input.calibration else None
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