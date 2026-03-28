configfile: "Config/config.yaml"

query = config["query"]
targets = config.get("targets")

if targets is None:
    targets = {"default": config["target"]}

run_rnahybrid_config = config.get("run_rnahybrid", {})
# tidy_rnahybrid_config = config.get("tidy_rnahybrid", {})
# annotate_rnahybrid = config.get("annotate_rnahybrid", {})

results_dir = config.get("results_dir", "Data/Results")


rule all:
    input:
        expand(f"{results_dir}" + "/{sample}/rnahybrid_enhanced.txt", sample=targets.keys())


rule run_rnahybrid:
    input:
        query=query,
        target=lambda wc: targets[wc.sample]
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
        distribution=run_rnahybrid_config.get("distribution")
    threads:
        int(run_rnahybrid_config.get("threads", 1))
    script:
        "Workflow/Scripts/rnahybrid.py"


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


rule annotate_rnahybrid:
    input:
        tidy=f"{results_dir}" + "/{sample}/tidy_output.csv",
        target=lambda wc: targets[wc.sample]
    output:
        annotated=f"{results_dir}" + "/{sample}/rnahybrid_annotated.csv"
    script:
        "Workflow/Scripts/annotate_rnahybrid.py"

rule enhance:
    input:
        annotated=f"{results_dir}" + "/{sample}/rnahybrid_annotated.csv"
    output:
        enhanced=f"{results_dir}" + "/{sample}/rnahybrid_enhanced.txt"
    script:
        "Workflow/Scripts/enhance_rnahybrid.py"