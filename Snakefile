configfile: "Config/config.yaml"

def format_optional_args(**kwargs):
    parts = []
    for flag, value in kwargs.items():
        if value is not None:
            parts.extend([flag, str(value)])
    return " ".join(parts)

query = config["query"]
target = config["target"]
species = config["species"]

rnahybrid_config = config.get("rnahybrid", {})
tidy_rnahybrid_config = config.get("tidy_rnahybrid", {})
outputs_config = config.get("outputs", {})

compact_output = outputs_config.get("compact", "Data/results/rnahybrid_compact.tsv")
tidy_output = outputs_config.get("tidy", "Data/results/rnahybrid_tidy.tsv")


rule all:
    input:
        tidy_output


rule rnahybrid:
    input:
        query=query,
        target=target
    output:
        compact=compact_output
    params:
        species=species,
        extra_args=format_optional_args(
            **{
                "--hits": rnahybrid_config.get("hits"),
                "-u": rnahybrid_config.get("u"),
                "-v": rnahybrid_config.get("v"),
                "--energy": rnahybrid_config.get("energy"),
                "--pvalue": rnahybrid_config.get("pvalue"),
                "--seed": rnahybrid_config.get("seed"),
            }
        )
    threads:
        int(rnahybrid_config.get("threads", 1))
    shell:
        """
        python Workflow/Scripts/rnahybrid.py \
          --query {input.query} \
          --target {input.target} \
          --species {params.species} \
          --output {output.compact} \
          --threads {threads} \
          {params.extra_args}
        """


rule tidy_rnahybrid:
    input:
        compact=compact_output
    output:
        tidy=tidy_output
    params:
        script="Workflow/Scripts/tidy_rnahybrid.R",
        extra_args=format_optional_args(
            **{
                "--mirna": tidy_rnahybrid_config.get("mirna"),
                "--pvalue-cutoff": tidy_rnahybrid_config.get("pvalue_cutoff"),
            }
        )
    shell:
        """
        Rscript {params.script} \
          --input {input.compact} \
          --output {output.tidy} \
          {params.extra_args}
        """
