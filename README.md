# Novum Pipeline

Minimal Snakemake workflow for optionally calibrating RNAhybrid statistics with RNAcalibrate, then running RNAhybrid on miRNA queries against bacterial coding sequences, tidying the output, and annotating the hits.

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

Edit [`Config/config.yaml`](/home/alber/Documents/novum_pipeline/Config/config.yaml) to change inputs, optional RNAcalibrate settings, RNAhybrid parameters, and output file names.

Example:

```yaml
query: Data/Raw/validated_miRNAs.fa
targets:
  escherichia_coli: Data/Raw/GCF_000005845.2_ASM584v2_cds_from_genomic_escherichia_coli.fna

calibration:
  enabled: false
  k: 5000
  max_target_length: 50000
  randomize_targets: false

rnahybrid:
  threads: 16
  species: 3utr_human
  hits: 3
  u: 3
  v: 3
  energy: -18
  pvalue: null
  seed: null
  distribution: 2.769859,0.233703

tidy_rnahybrid: {}

annotate_rnahybrid: {}

results_dir: Data/Results/run_001
```

When `calibration.enabled: true`, the workflow writes one calibration artifact per target at `.../{sample}/rnacalibrate.json`. RNAhybrid will read the mean `xi/theta` distribution from that file and override the static `rnahybrid.distribution` value for that sample. `calibration.randomize_targets` maps to the RNAcalibrate `-s` flag and should usually stay `false` unless you have confirmed it produces valid fits for your inputs.
