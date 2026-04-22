# Novum Pipeline

Minimal Snakemake workflow for optionally calibrating RNAhybrid statistics with RNAcalibrate, then running RNAhybrid on miRNA queries against bacterial coding sequences, tidying the output, annotating the hits, enhancing them with alignment visualizations, and optionally rendering ggplot2 summary PDFs.

## Workflow

```text
         query.fa ─────────────────────┐
            │                          │
            ▼                          ▼
   ┌─────────────────┐        ┌─────────────────┐
   │  rnacalibrate   │ ─────▶ │    rnahybrid    │ 
   │  (optional)     │  xi/θ  │  (GNU Parallel) │
   └─────────────────┘        └────────┬────────┘
                                       │
                                       ▼
                              ┌─────────────────┐
                              │ tidy_rnahybrid  │
                              │    (Python)     │
                              └────────┬────────┘
                                       │
                                       ▼
                              ┌──────────────────┐
                              │annotate_rnahybrid│
                              │    (Python)      │
                              └────────┬─────────┘
                                       │
                            ┌──────────┴──────────┐
                            ▼                     ▼
                  ┌──────────────────┐  ┌──────────────────┐
                  │enhance_rnahybrid │  │   build_plots    │
                  │     (Python)     │  │  (R, optional)   │
                  └──────────────────┘  └──────────────────┘
```

1. **rnacalibrate** (optional) — estimates extreme-value distribution parameters (xi, theta) from the target sequence.
2. **rnahybrid** — runs RNAhybrid in parallel via GNU Parallel across all query/target pairs.
3. **tidy_rnahybrid** — filters, selects, and arranges raw RNAhybrid output into a tidy CSV.
4. **annotate_rnahybrid** — parses FASTA headers from the target genome and merges metadata with the tidy output.
5. **enhance_rnahybrid** — adds human-readable alignment visualizations to the annotated results.
6. **build_plots** (optional) — renders configurable ggplot2 PDFs (p-value distribution, position-vs-p-value, position-vs-energy) from the annotated tables. Skipped entirely when `build_plots.type` is unset.

## Dependencies

Commands below assume a Linux bash shell on Debian/Ubuntu.

### GNU Parallel

```bash
sudo apt update && sudo apt install -y parallel
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

Edit [`Config/config.yaml`](Config/config.yaml) to change inputs, optional RNAcalibrate settings, RNAhybrid parameters, and output file names.

Example:

```yaml
query: Data/Raw/validated_miRNAs.fa
targets:
  escherichia_coli: Data/Raw/GCF_000005845.2_ASM584v2_cds_from_genomic_escherichia_coli.fna

rnacalibrate:
  mode: both                # "calibrated" | "uncalibrated" | "both"
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

build_plots:                  # optional; omit the whole block to skip plotting
  type: all                   # pvalue_distribution | position_pvalue | position_energy | all | comma-separated subset
  basesize: 12
  pvalue: 0.01                # cutoff for position_pvalue and position_energy
  locus:                      # validated-hit locus_tags (highlighted in the plots)
    - LGG_RS02140
    - LGG_RS05490
  gene:                       # validated-hit gene_names (highlighted in the plots)
    - spaC
    - lexA
```

`rnacalibrate.mode` selects which variant(s) the pipeline produces in a single invocation:

- `calibrated` — runs `rnacalibrate` first and feeds its xi/theta estimate into `rnahybrid`. Outputs land under `{results_dir}/{sample}/w_calibration/`.
- `uncalibrated` — skips `rnacalibrate` entirely; `rnahybrid` falls back to `rnahybrid.distribution` (or its built-in default if unset). Outputs land under `{results_dir}/{sample}/wo_calibration/`.
- `both` — produces both trees in one run, so the two can be compared side by side.

When a calibrated variant runs, the calibration artifact is written at `{results_dir}/{sample}/w_calibration/rnacalibrate.json` and `rnahybrid` reads the estimated xi/theta distribution from it, overriding the static `rnahybrid.distribution` value for that sample. The `rnacalibrate.randomize_targets` option maps to the RNAcalibrate `-s` flag and should usually stay `false` unless you have confirmed it produces valid fits for your inputs.

> **Deprecated:** the older `rnacalibrate.enabled: true|false` flag is still honored when `mode` is absent (`true` → `calibrated`, `false` → `uncalibrated`), but new configs should use `mode`. Support for `enabled` will be removed in a future release.

The optional `build_plots` block controls a final ggplot2 rule that consumes the annotated CSVs and produces a single PDF per sample. The rule is opt-in: omit the block (or remove `type:`) and no plot job is scheduled, no R env is materialized.

- `type` selects which panels appear in the PDF. Valid values: `pvalue_distribution`, `position_pvalue`, `position_energy`, `all`, or any comma-separated subset (e.g. `"pvalue_distribution,position_pvalue"`). Panels appear in the order listed.
- `basesize` sets the ggplot2 base font size.
- `pvalue` is the cutoff applied to the `position_pvalue` and `position_energy` scatter points (the `pvalue_distribution` density plot ignores it).
- `locus` and `gene` are lists of `locus_tag` and `gene_name` values; matching rows are highlighted as "validated hits" with their own color/shape and listed in a shared legend.

Output filenames bake the type into the name (`plots_<type>.pdf`, with commas in `type` rendered as hyphens), so cycling through types in the same `results_dir` doesn't overwrite earlier PDFs. When `rnacalibrate.mode: both`, each plot is faceted by calibration variant; when only one variant is produced, the facet collapses to a single panel automatically.

## Citation

If you use this workflow in a publication, please cite this repository and the external tools it relies on:

- **RNAhybrid**: Kruger, J. and Rehmsmeier, M. (2006). RNAhybrid: microRNA target prediction easy, fast and flexible. *Nucleic Acids Research*, 34(Web Server issue), W451--W454. https://doi.org/10.1093/nar/gkl243
- **GNU Parallel**: Tange, O. (2023, November 22). GNU Parallel 20231122 ('Grindavik'). Zenodo. https://doi.org/10.5281/zenodo.10199085
