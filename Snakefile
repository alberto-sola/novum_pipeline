rule all:
    input:
        "output.tsv"

rule generate:
    output:
        "output.tsv"
    shell:
        "echo 'Hello Snakemake' > {output}"