# Novum Pipeline

Minimal Snakemake workflow for running RNAhybrid on miRNA queries against bacterial coding sequences, then tidying and annotating the results.

## Dependencies

Commands below assume a Linux bash shell on Debian/Ubuntu.

### GNU Parallel

```bash
sudo apt update
sudo apt install -y parallel
```

### Conda

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/miniconda3/bin/activate
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority strict
```

This produces the recommended channel order:

```text
bioconda
conda-forge
defaults
```

### RNAhybrid

```bash
conda install bioconda::rnahybrid
```

### Snakemake

```bash
conda create -n snakemake-modern python=3.11 snakemake-minimal
conda activate snakemake-modern
```

## Run

```bash
snakemake --use-conda --cores N
```

## Configuration

Edit [`Config/config.yaml`](/home/alber/Documents/novum_pipeline/Config/config.yaml) to change inputs, RNAhybrid parameters, and output file names.

Example:

```yaml
query: Data/Raw/validated_miRNAs.fa
target: Data/Raw/GCF_000005845.2_ASM584v2_cds_from_genomic_escherichia_coli.fna

run_rnahybrid:
  threads: 16
  species: 3utr_human
  hits: 3
  u: 3
  v: 3
  energy: -18
  pvalue: null
  seed: null

tidy_rnahybrid: {}

annotate_rnahybrid: {}

outputs:
  compact: Data/Results/rnahybrid_output.tsv
  tidy: Data/Results/tidy_output.csv
  annotated: Data/Results/rnahybrid_named.csv
```
