configfile: "Config/config.yaml"

query = config["query"]
target = config["target"]
species = config.get("species", "3utr_human")

rnahybrid_config = config.get("rnahybrid", {})
outputs_config = config.get("outputs", {})

compact_output = outputs_config.get("compact", "Data/results/rnahybrid_compact.tsv")
tidy_output = outputs_config.get("tidy", "Data/results/rnahybrid_tidy.csv")


rule all:
    input:
        tidy_output


rule rnahybrid:
    input:
        query=query,
        target=target
    output:
        compact=compact_output
    conda:
        "Workflow/Envs/rnahybrid.yaml"
    params:
        species=species,
        hits=rnahybrid_config.get("hits"),
        u=rnahybrid_config.get("u"),
        v=rnahybrid_config.get("v"),
        energy=rnahybrid_config.get("energy"),
        pvalue=rnahybrid_config.get("pvalue"),
        seed=rnahybrid_config.get("seed")
    threads:
        int(rnahybrid_config.get("threads", 1))
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
        Rscript {params.script} \
          --input {input.compact} \
          --output {output.tidy}
        """
