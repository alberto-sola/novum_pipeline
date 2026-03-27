configfile: "Config/config.yaml"

query = config["query"]
target = config["target"]

run_rnahybrid_config = config.get("run_rnahybrid", {})
# tidy_rnahybrid_config = config.get("tidy_rnahybrid", {})
# annotate_rnahybrid = config.get("annotate_rnahybrid", {})

outputs_config = config.get("outputs", {})
compact_output = outputs_config.get("compact", "Data/Results/rnahybrid_compact.tsv")
tidy_output = outputs_config.get("tidy", "Data/Results/rnahybrid_tidy.csv")
annotated_output = outputs_config.get("annotated", "Data/Results/rnahybrid_annotated.csv")
enhanced_output = outputs_config.get("enhanced", "Data/Results/rnahybrid_enhanced.txt")


rule all:
    input:
        enhanced_output


rule run_rnahybrid:
    input:
        query=query,
        target=target
    output:
        compact=compact_output
    conda:
        "Workflow/Envs/rnahybrid.yaml"
    params:
        species=run_rnahybrid_config.get("species", "3utr_human"),
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
        compact=compact_output
    output:
        tidy=tidy_output
    params:
        script="Workflow/Scripts/tidy_rnahybrid.R"
    shell:
        """
        Rscript {params.script} --input {input.compact} --output {output.tidy}
        """


rule annotate_rnahybrid:
    input:
        tidy=tidy_output,
        target=target
    output:
        annotated=annotated_output
    script:
        "Workflow/Scripts/annotate_rnahybrid.py"

rule enhance:
    input:
        annotated=annotated_output
    output:
        enhanced=enhanced_output
    script:
        "Workflow/Scripts/enhance_rnahybrid.py"
