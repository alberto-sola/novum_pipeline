#----- Handle to the config file -----#
configfile: "Config/config.yaml"

#----- Populate Snakefile variables with the config file -----#
query = config["query"]
targets = config["targets"]
run_rnacalibrate_config = config["calibration"]
run_rnahybrid_config = config["rnahybrid"]
results_dir = config.get("results_dir", "Data/Results")
calibration_enabled = bool(run_rnacalibrate_config.get("enabled", False))


def get_calibration_output(sample):
    return f"{results_dir}" + f"/{sample}/rnacalibrate.json"


#----- output finale ‒ pipeline conclusion -----#
rule all:
    input:
        expand(f"{results_dir}" + "/{sample}/rnahybrid_enhanced.txt", sample=targets.keys())

#----- Dynamically calibrate the statistics based on the target sequence -----#
rule run_rnacalibrate:
    input:
        query=query,
        target=lambda wc: targets[wc.sample]
    output:
        calibration=f"{results_dir}" + "/{sample}/rnacalibrate.json"
    conda:
        "Workflow/Envs/rnahybrid.yaml"
    params:
        k=run_rnacalibrate_config.get("k"),
        max_target_length=run_rnacalibrate_config.get("max_target_length", 50000),
        randomize_targets=run_rnacalibrate_config.get(
            "randomize_targets",
            run_rnacalibrate_config.get("use_target_distribution", False),
        ),
        u=run_rnahybrid_config.get("u"),
        v=run_rnahybrid_config.get("v"),
        seed=run_rnahybrid_config.get("seed")
    script:
        "Workflow/Scripts/rnacalibrate.py"

#----- Run RNAhybrid on multiple cores -----#
rule run_rnahybrid:
    input:
        query=query,
        target=lambda wc: targets[wc.sample],
        calibration=lambda wc: get_calibration_output(wc.sample) if calibration_enabled else []
    output:
        compact=f"{results_dir}" + "/{sample}/rnahybrid_output.tsv"
    conda:
        "Workflow/Envs/rnahybrid.yaml"
    params:
        species=run_rnahybrid_config.get("species"),
        hits=run_rnahybrid_config.get("hits"),
        u=run_rnahybrid_config.get("u"),
        v=run_rnahybrid_config.get("v"),
        energy=run_rnahybrid_config.get("energy"),
        pvalue=run_rnahybrid_config.get("pvalue"),
        seed=run_rnahybrid_config.get("seed"),
        distribution=run_rnahybrid_config.get("distribution"),
        distribution_file=lambda wc, input: input.calibration if input.calibration else None
    threads:
        int(run_rnahybrid_config.get("threads", 1))
    script:
        "Workflow/Scripts/rnahybrid.py"

#----- Filter, select, arrange the RNAhybrid's output -----#
rule tidy_rnahybrid:
    input:
        compact=f"{results_dir}" + "/{sample}/rnahybrid_output.tsv"
    output:
        tidy=f"{results_dir}" + "/{sample}/tidy_output.csv"
    params:
        script="Workflow/Scripts/tidy_rnahybrid.R"
    shell:
        """
        Rscript {params.script} --input {input.compact} --output {output.tidy}
        """

#----- Parse metadata from the target file (genome) and merge them with the tidied output -----#
rule annotate_rnahybrid:
    input:
        tidy=f"{results_dir}" + "/{sample}/tidy_output.csv",
        target=lambda wc: targets[wc.sample]
    output:
        annotated=f"{results_dir}" + "/{sample}/rnahybrid_annotated.csv"
    script:
        "Workflow/Scripts/annotate_rnahybrid.py"

#----- Visualize RNAhybrid alignments + all the data -----#
rule enhance_rnahybrid:
    input:
        annotated=f"{results_dir}" + "/{sample}/rnahybrid_annotated.csv"
    output:
        enhanced=f"{results_dir}" + "/{sample}/rnahybrid_enhanced.txt"
    script:
        "Workflow/Scripts/enhance_rnahybrid.py"
