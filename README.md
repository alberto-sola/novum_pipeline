# Novum Pipeline

A Snakemake workflow that runs **RNAhybrid and IntaRNA in parallel** on miRNA queries against bacterial coding sequences, RNAhybrid optionally per-miRNA-calibrated via RNAcalibrate, IntaRNA optionally with accessibility correction, each arm producing its own tidied, annotated, and enhanced output tree. Cross-referencing the two tools' outputs per interaction pair is planned but not yet implemented.

## Workflow

```text
                          ┌── (RNAcalibrate) ─ RNAhybrid ─ tidy ─ annotate ─ enhance ─┬─ (plots)
queries + targets ────────┤                                                           │
                          └── IntaRNA ──────────────────── tidy ─ annotate ─ enhance ─┘  (no plots yet)
```

### RNAhybrid arm

1. **rnacalibrate** *(optional)* — fits an extreme-value distribution per query miRNA, emitting one (xi, theta) pair each.
2. **rnahybrid** — runs in parallel via GNU Parallel. In calibrated mode each (miRNA, target chunk) pair uses that miRNA's own (xi, theta).
3. **tidy_rnahybrid** — filters and arranges the raw output into a tidy CSV.
4. **annotate_rnahybrid** — parses FASTA headers from the target genome and merges metadata onto the tidy rows. Shares `Workflow/Scripts/annotate.py` with the IntaRNA arm (the `insert_after` column is parameterized per rule).
5. **enhance_rnahybrid** — appends human-readable alignment visualizations.
6. **build_plots** *(optional, RNAhybrid only)* — renders configurable ggplot2 PDFs. Skipped entirely when `build_plots.type` is unset. No IntaRNA / cross-arm plot stage exists yet.

### IntaRNA arm

7. **intarna** — runs IntaRNA with native `--threads`; no external chunking. The `wo_accessibility` variant locks the §10.1 RNAhybrid-emulation profile (`--noSeed --acc=N --intLoopMax=30 --mode=M`) for an apples-to-apples comparison.
8. **tidy_intarna** — filters and renames columns; computes `Position = Start1 / Gene_length` for positional comparability with the RNAhybrid arm.
9. **annotate_intarna** — parses FASTA headers and merges metadata (shares `annotate.py`, inserting after the `E` column).
10. **enhance_intarna** — per-record human-readable report with energy, accessibility, and seed metadata plus a duplex block rendered from `subseqDP`/`hybridDP`.

## Quick start

If you already have conda and Snakemake set up (otherwise see [Dependencies](#dependencies)):

1. Drop your input FASTAs into `Data/Raw/` and edit `Config/config.yaml` so each `queries:` and `targets:` key points at them. The two mappings must share keys.
2. Run:

   ```bash
   snakemake --use-conda --cores 8
   ```

3. Results land under `Data/Results/<sample>/`. Both arms run: RNAhybrid outputs appear under `{w,wo}_calibration/` and IntaRNA outputs under `{w,wo}_accessibility/`. See [Outputs](#outputs) for the full layout.

Need a graphical config editor? Skip ahead to the [UI](#config-editor-ui--optional) once you've installed dependencies.

## Dependencies

Commands assume a Linux bash shell on Debian/Ubuntu.

### Conda

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/miniconda3/bin/activate
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

Channel order matters: bioconda is built on top of conda-forge, so conda-forge must have higher priority to avoid ABI mismatches.

```text
conda-forge
bioconda
defaults
```

### Snakemake

```bash
conda create -n snakemake-modern python=3.11 snakemake-minimal
conda activate snakemake-modern
```

### Config editor (UI) — optional

A PyWebView desktop app at `UI/launcher.py` edits `Config/config.yaml` graphically and can kick off the pipeline directly. Install its extras into the same env that has `snakemake` on PATH, then launch:

```bash
conda activate snakemake-modern
pip install pywebview pyqt5 pyqtwebengine
python UI/launcher.py
```

The launcher resolves the repo root from its own location, so it works from any working directory. Snakemake's stdout/stderr are inherited by the terminal you launched it from — that's where to watch progress when you trigger runs from the UI.

## Configuration

Edit [`Config/config.yaml`](Config/config.yaml). It has five blocks: `queries`/`targets` (inputs), `rnacalibrate` (calibration mode + parameters), `rnahybrid` (RNAhybrid options), `intarna` (IntaRNA options, new), and `build_plots` (plotting, optional).

Example:

```yaml
queries:
  escherichia: Data/Raw/RNAs/validated_escherichia.fa
  bacteroides: Data/Raw/RNAs/validated_bacteroides.fa

targets:
  escherichia: Data/Raw/genomes/GCF_000005845.2_ASM584v2_cds_from_genomic_escherichia_coli.fna
  bacteroides: Data/Raw/genomes/GCF_014131755.1_ASM1413175v1_cds_from_genomic_bacteroides_thetaiotaomicron.fna

rnacalibrate:
  mode: both
  k: 10000
  max_target_length: 50000
  randomize_targets: true

rnahybrid:
  threads: 16
  species: 3utr_human
  hits: null
  u: null
  v: null
  energy: -18
  pvalue: null
  seed: null
  distribution: null

intarna:
  acc_mode: both
  threads: 16
  mode: H
  model: S
  int_len_max: 0
  int_loop_max: 10
  seed:
    bp: 7
    max_E: null
    max_E_hybrid: null
    min_Pu: null
    no_GU: false
    no_GU_end: false
    q_range: null
    t_range: null
    max_UP: null
    out_best_only: false
  accessibility:
    w: null
    l: null
    no_lp: false
    no_gu_end: false
  output:
    max_E: -18
    delta_E: null
    number: 1
    overlap: B
    min_Pu: null
    no_lp: false
    no_gu_end: false
    csv_cols: id1,id2,start1,end1,start2,end2,subseqDP,hybridDP,E,E_hybrid,ED1,ED2,Pu1,Pu2,seedStart1,seedEnd1,seedE,seedStart2,seedEnd2
  extra_args: []

build_plots:
  type: position_pvalue,pvalue_distribution
  basesize: 12
  pvalue: 0.01
  locus:
    - FE838_RS16060
    - FE838_RS16090
    - FE838_RS16070
    - FE838_RS16075
    - FE838_RS16085
    - FE838_RS16065
    - FE838_RS16080
  gene:
    - yegH,hsa-miR-1226-5p
    - rnpA,hsa-miR-4747-3p
  protein:

results_dir: Data/Results/
```

### Queries and targets

`queries:` and `targets:` are paired by key — each key names one (query, target) sample run, and the two mappings must use identical keys. The legacy single `query: <path>` form is still accepted: when `queries:` is absent, the same FASTA is broadcast against every entry in `targets:`.

### Calibration mode

`rnacalibrate.mode` selects which variant(s) the pipeline produces in a single invocation:

- **`calibrated`** — runs `rnacalibrate` first and feeds its per-miRNA (xi, theta) into `rnahybrid`. Outputs land under `{results_dir}/{sample}/w_calibration/`.
- **`uncalibrated`** — skips `rnacalibrate`; `rnahybrid` falls back to `rnahybrid.distribution`, or to its built-in `species` default if unset. Outputs land under `{results_dir}/{sample}/wo_calibration/`.
- **`both`** — produces both trees in one run for side-by-side comparison.

When calibration runs, the artifact at `{results_dir}/{sample}/w_calibration/rnacalibrate.json` contains one (xi, theta) row per miRNA in the query FASTA. RNAhybrid is then invoked once per (miRNA, target chunk) pair with that miRNA's own parameters, so every row's p-value reflects its own null model rather than a population mean. The static `rnahybrid.distribution` value is ignored whenever calibration runs.

`rnacalibrate.randomize_targets` maps to RNAcalibrate's `-s` flag and should usually stay `false` unless you have confirmed it produces valid fits for your inputs.

> **Deprecated:** `rnacalibrate.enabled: true|false` is still honored when `mode` is absent (`true` → `calibrated`, `false` → `uncalibrated`), but new configs should use `mode`. Support will be removed in a future release.

### IntaRNA & accessibility mode

`intarna.acc_mode` selects which accessibility variant(s) the IntaRNA arm produces — analogous to `rnacalibrate.mode` on the RNAhybrid side:

- **`on`** — runs IntaRNA with RNA accessibility correction enabled (`--acc=C`). Outputs land under `{results_dir}/{sample}/w_accessibility/`. The `mode`, `model`, `int_len_max`, `int_loop_max`, `seed.*`, and `accessibility.*` config keys all apply.
- **`off`** — runs the **§10.1 RNAhybrid-emulation profile**: IntaRNA is locked to `--noSeed --acc=N --intLoopMax=30 --mode=M`. The config keys `mode`, `model`, `int_len_max`, `int_loop_max`, `seed.*`, and `accessibility.*` are **silently ignored** for this variant. Only `output.*` and `extra_args` are honored. Outputs land under `{results_dir}/{sample}/wo_accessibility/`. This is the canonical apples-to-apples profile for comparison against RNAhybrid: no seed constraint, no accessibility correction, exact mode, energy threshold comparable to `rnahybrid.energy`.
- **`both`** — produces both `w_accessibility/` and `wo_accessibility/` trees in one run.

`intarna.threads` maps to IntaRNA's native `--threads` flag; the IntaRNA arm does not use GNU Parallel (unlike the RNAhybrid arm, which chunks targets externally).

`output.csv_cols` is the column whitelist passed to IntaRNA via `--outCsvCols`. The downstream `tidy_intarna` rule requires at minimum `id1, id2, start1, end1, start2, end2, E`; trimming `csv_cols` below that will break the IntaRNA arm. The `hybridDP` and `subseqDP` columns are also required if you want `enhance_intarna` to render alignment visualizations.

> **YAML 1.1 gotcha:** bare `on` and `off` are parsed by PyYAML as Python `True`/`False`. The Snakefile coerces them back to strings, so `acc_mode: on` works — but `acc_mode: "on"` (quoted) is more portable if you feed the config to other tools.

Cross-referencing the RNAhybrid and IntaRNA output tables (per-pair agreement, ranking deltas, etc.) is the next planned milestone on this branch and is not yet implemented.

### Plots

The `build_plots` block is opt-in: omit it (or remove `type:`) and no plot job is scheduled, no R env is materialized. Plots currently consume only the RNAhybrid `rnahybrid_annotated.csv` outputs; an IntaRNA / cross-referenced plot stage is future work.

To re-render plots after editing `build_plots` without re-running the upstream pipeline:

```bash
snakemake --use-conda --cores N --forcerun build_plots
```

- **`type`** — which panels to render: `pvalue_distribution`, `position_pvalue`, `position_energy`, `all`, or any comma-separated subset (e.g. `"pvalue_distribution,position_pvalue"`). Panels appear in the order listed.
- **`basesize`** — ggplot2 base font size.
- **`pvalue`** — cutoff applied to the scatter points in `position_pvalue` and `position_energy` (the density plot ignores it).
- **`locus`**, **`gene`**, **`protein`** — lists of validated hits to highlight. Each entry is either:
    - **bare** — `<target>` matches any miRNA hitting that target;
    - **paired** — `<target>,<miRNA>` matches only when the row's miRNA equals the named one.

  Bare and paired entries can mix freely in the same list. YAML strings with spaces or punctuation (typical for `protein_name`) must be quoted: `"30S ribosomal protein S2,hsa-miR-X"`. The `protein:` block is optional and requires a `protein_name` column in the annotated CSV.

The `pvalue_distribution` panel draws one dashed vertical line per validated entity (lowest p-value per `(miRNA, locus_tag)`); `position_pvalue` and `position_energy` plot every alignment of the matched pairs, all keyed to the same color and shape in a shared legend.

The output filename bakes the `type` value in (commas become hyphens), so cycling through types in the same `results_dir` doesn't overwrite earlier PDFs. When `rnacalibrate.mode: both`, panels are faceted by calibration variant; with one variant, the facet collapses automatically.

## Outputs

One tree per `targets:` key, under `{results_dir}/{sample}/`:

```text
{sample}/
├── w_calibration/                  # rnacalibrate.mode = "calibrated" or "both"
│   ├── rnacalibrate.json
│   ├── rnahybrid_output.tsv
│   ├── tidy_output.csv
│   ├── rnahybrid_annotated.csv
│   └── rnahybrid_enhanced.txt
├── wo_calibration/                 # rnacalibrate.mode = "uncalibrated" or "both"
│   ├── rnahybrid_output.tsv
│   ├── tidy_output.csv
│   ├── rnahybrid_annotated.csv
│   └── rnahybrid_enhanced.txt
├── w_accessibility/                # intarna.acc_mode = "on" or "both"
│   ├── intarna_output.csv
│   ├── intarna_tidy.csv
│   ├── intarna_annotated.csv
│   └── intarna_enhanced.txt
├── wo_accessibility/               # intarna.acc_mode = "off" or "both"
│   ├── intarna_output.csv          # locked to the §10.1 RNAhybrid-emulation profile
│   ├── intarna_tidy.csv
│   ├── intarna_annotated.csv
│   └── intarna_enhanced.txt
└── plots_<type>.pdf                # only when build_plots.type is set; RNAhybrid arm only
```

Plots currently consume only `rnahybrid_annotated.csv`; an IntaRNA / cross-referenced plot stage is future work.

## Citation

If you use this workflow in a publication, please cite this repository and the external tools it relies on:

- **RNAhybrid / RNAcalibrate**: Rehmsmeier, M., Steffen, P., Höchsmann, M., and Giegerich, R. (2004). Fast and effective prediction of microRNA/target duplexes. *RNA*, 10(10), 1507–1517. https://doi.org/10.1261/rna.5248604
- **RNAhybrid (web server)**: Krüger, J. and Rehmsmeier, M. (2006). RNAhybrid: microRNA target prediction easy, fast and flexible. *Nucleic Acids Research*, 34(Web Server issue), W451–W454. https://doi.org/10.1093/nar/gkl243
- **IntaRNA**: Mann, M., Wright, P. R., and Backofen, R. (2017). IntaRNA 2.0: enhanced and customizable prediction of RNA–RNA interactions. *Nucleic Acids Research*, 45(W1), W435–W439. https://doi.org/10.1093/nar/gkx279
- **Snakemake**: Mölder, F., Jablonski, K. P., Letcher, B., et al. (2021). Sustainable data analysis with Snakemake. *F1000Research*, 10, 33. https://doi.org/10.12688/f1000research.29032.2
- **GNU Parallel**: Tange, O. (2011). GNU Parallel: The command-line power tool. *;login: The USENIX Magazine*, 36(1), 42–47.
