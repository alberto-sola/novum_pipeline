# Novum Pipeline

Novum Pipeline predicts which bacterial genes a set of miRNAs are likely to bind. It runs two well-established RNA–RNA interaction tools — **RNAhybrid** and **IntaRNA** — side by side over the same miRNAs and genomes, then reports what each tool finds on its own together with the **consensus**: the miRNA–gene pairs both tools agree on.

The whole process is wired together with [Snakemake](https://snakemake.github.io/), so a single command takes you from raw FASTA sequences all the way to annotated, human-readable predictions. You point it at your miRNAs and your target genomes, and it takes care of the rest.

New here? Install the [dependencies](#dependencies), then jump to the [quick start](#quick-start). The [How it works](#how-it-works) section at the end explains the two arms and the consensus in more detail.

## Dependencies

Everything the pipeline needs — RNAhybrid, IntaRNA, the R plotting stack — is installed through conda automatically: Snakemake builds an isolated environment for each step on the first run. All you set up by hand is **conda** and **Snakemake**. Commands assume a bash shell.

### Conda

*(Skip this if you already have conda.)* Download and install Miniconda — pick the installer for your platform:

```bash
# Linux
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# macOS (Apple Silicon) — use the arm64 installer
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
bash Miniconda3-latest-MacOSX-arm64.sh
```

Then activate it and configure the channels (identical on every platform):

```bash
source ~/miniconda3/bin/activate
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

### Snakemake

Create and activate the environment you'll run the pipeline from:

```bash
conda create -n snakemake-modern python=3.11 snakemake-minimal
conda activate snakemake-modern
```

### Apple Silicon Macs — one extra step

> **This step is for M-series (arm64) Macs only.** On Linux and Intel Macs, skip straight to the [quick start](#quick-start).

RNAhybrid has no Apple-Silicon (`osx-arm64`) conda build, so on M-series Macs the pipeline resolves its whole conda toolchain as **`osx-64` (Intel)** and runs it under **Rosetta 2**. A helper script sets that up in one shot — run it once, with your env active:

```bash
conda activate snakemake-modern
./setup-macos.sh
```

`setup-macos.sh` installs Rosetta 2 if it is missing (and skips it if already present), then pins this env's conda subdir to `osx-64`. It is a safe no-op on Linux and Intel Macs.

- It's a **one-time step per env** — the `osx-64` setting persists in the env, so re-run it only if you recreate the env.
- The **first** pipeline run builds the Intel environments and is slower; later runs are normal. Steady-state is roughly 10–25% slower than native ARM, and only on the IntaRNA arm.
- Run it **before** your first pipeline run. If you already started one, delete `.snakemake/conda/` so the per-rule environments rebuild cleanly as `osx-64`.

### Config editor (optional)

The pipeline reads all its settings from `Config/config.yaml`, which you can edit by hand (see [Configuration](#configuration)) or through a small desktop app. The app is entirely optional — install it only if you want a graphical editor:

```bash
conda activate snakemake-modern
pip install pywebview pyqt5 pyqtwebengine
```

See [Config editor (UI)](#config-editor-ui) in the quick start for how to launch it.

## Quick start

With the dependencies installed and your env active (`conda activate snakemake-modern`):

1. **Point the pipeline at your inputs.** Open `Config/config.yaml` and fill in the `queries:` block (your miRNA FASTAs) and the `targets:` block (the genome CDS FASTAs to search), using one shared key per sample:

   ```yaml
   queries:
     my_sample: Data/Raw/RNAs/my_mirnas.fa
   targets:
     my_sample: Data/Raw/genomes/my_genome_cds.fna
   ```

   Drop your own FASTAs under `Data/Raw/`, or let the pipeline assemble them for you from `(taxon, miRNA)` pairs — see [Preparing inputs automatically](#preparing-inputs-automatically). Every other setting has a sensible default, so you can leave the rest of the file untouched for a first run.

2. **Run it:**

   ```bash
   snakemake --use-conda --cores 8
   ```

   The first run is slower because Snakemake builds one conda environment per tool; later runs reuse them.

3. **Read the results** under `Data/Results/<sample>/`. See [Outputs](#outputs) for the layout.

### Config editor (UI)

Prefer not to edit YAML by hand? Launch the optional desktop editor (install it first — see [Config editor (optional)](#config-editor-optional)):

```bash
python UI/launcher.py
```

It edits `Config/config.yaml` graphically and can start a run for you with a single click. It locates the repo from its own path, so you can launch it from any directory. When you trigger a run, Snakemake's output is written to `Data/Results/.pipeline.log` and shown as a live tail inside the app.

## Outputs

Each `(query, target)` sample gets its own results tree under `{results_dir}/{sample}/`. With every variant enabled, the full layout is:

```text
{sample}/
├── w_calibration/
│   ├── rnacalibrate.json
│   ├── rnahybrid_output.tsv
│   ├── tidy_output.csv
│   ├── rnahybrid_annotated.csv
│   └── rnahybrid_enhanced.txt
├── wo_calibration/
│   ├── rnahybrid_output.tsv
│   ├── tidy_output.csv
│   ├── rnahybrid_annotated.csv
│   └── rnahybrid_enhanced.txt
├── w_accessibility/
│   ├── intarna_output.csv
│   ├── intarna_tidy.csv
│   ├── intarna_annotated.csv
│   └── intarna_enhanced.txt
├── wo_accessibility/
│   ├── intarna_output.csv
│   ├── intarna_tidy.csv
│   ├── intarna_annotated.csv
│   └── intarna_enhanced.txt
├── consensus/
│   ├── consensus_annotated.csv
│   └── consensus_enhanced.txt
└── plots_<type>.pdf
```

A given run fills only the sub-folders for the variants it produces (see [Configuration](#configuration)):

- **`w_calibration/`**, **`wo_calibration/`** — the RNAhybrid arm's calibrated and uncalibrated variants. `rnacalibrate.json` appears only in the calibrated tree.
- **`w_accessibility/`**, **`wo_accessibility/`** — the IntaRNA arm's accessibility-on and accessibility-off variants.
- **`consensus/`** — always produced; the `(miRNA, target)` pairs both tools predict.
- **`plots_<type>.pdf`** — only when `plots.type` is set (RNAhybrid arm only).

Within each arm's folder the files build on one another: the raw tool output → a tidied table → the same table with gene metadata merged in (`*_annotated.csv`) → a human-readable report with alignment/duplex visualizations (`*_enhanced.txt`).

## Configuration

All behaviour is controlled by [`Config/config.yaml`](Config/config.yaml), organised into labelled blocks. Here is a map of what each one does — an overview, not the full file. The tunable knobs inside each block are detailed in [Options in detail](#options-in-detail).

```yaml
#----- miRNAs and genomes (one query FASTA per target, joined by key) -----#
queries:
targets:
```

Your inputs. `queries:` and `targets:` are matched by key — each key names one `(miRNA set, genome)` sample, and both blocks must use the same keys. **This is the only block you must fill in.**

```yaml
#----- Shared by both arms -----#
threads:
max_suboptimal_hits:
max_total_energy:
seed:
```

Settings both tools obey, kept in one place so the two arms can't drift apart: CPU threads, how many hits to keep per pair, the energy cutoff, and the seed constraint. Note that `max_total_energy` is not measured identically by the two tools — see [Shared settings](#shared-settings).

```yaml
#----- Parameters for RNAcalibrate -----#
rnacalibrate:
```

Turns per-miRNA statistical calibration on or off for the RNAhybrid arm and controls how its null model is fit. See [Calibration](#calibration-rnahybrid-arm).

```yaml
#----- Parameters for RNAhybrid -----#
rnahybrid:
```

RNAhybrid's own parameters — background species model, loop-size limits, p-value threshold, and a static energy distribution used when calibration is off.

```yaml
#----- Parameters for IntaRNA -----#
intarna:
```

IntaRNA's parameters, including whether to apply accessibility correction, the prediction mode and model, and detailed seed / accessibility / output sub-blocks. See [Accessibility](#accessibility-intarna-arm) and [Seed handling](#seed-handling).

```yaml
#----- Parameters for plots (optional) -----#
plots:
```

Opt-in plotting. Leave `type:` unset and no plots are produced (and no R environment is built). See [Plots](#plots).

```yaml
#----- Where outputs should be generated -----#
results_dir:
```

Where result trees are written (default `Data/Results/`).

### Providing inputs (queries and targets)

`queries:` and `targets:` are paired by key — each key names one `(query, target)` sample, and the two blocks must use identical keys:

```yaml
queries:
  escherichia: Data/Raw/RNAs/validated_escherichia.fa
targets:
  escherichia: Data/Raw/genomes/GCF_000005845.2_ASM584v2_cds_from_genomic_escherichia_coli.fna
```

The legacy single `query: <path>` form is still accepted: when `queries:` is absent, that one FASTA is broadcast against every entry in `targets:`.

### Preparing inputs automatically

Rather than hand-collect FASTAs, the `data-prep/` helper can assemble them from a list of `(taxon, miRNA)` pairs: it downloads bacterial CDS assemblies from NCBI and miRNA sequences from miRBase, then rewrites just the `queries:`/`targets:` blocks of `Config/config.yaml` (every other key and comment is preserved). It never runs Snakemake — you keep full control.

```bash
# one pair per line: "<taxon> <miRNA>", taxon first (tab, ':' or spaces separate)
python data-prep/prepare_inputs.py --pairs pairs.txt
python data-prep/prepare_inputs.py --pairs pairs.txt --dry-run          # preview only, no downloads/writes
python data-prep/prepare_inputs.py Veillonella_parvula:hsa-miR-200b-3p  # inline pair(s)
```

A taxon reaches the config only if its genome **and** at least one of its miRNAs resolved; unresolved entries are skipped and recorded in a report TSV rather than aborting the run. See [`data-prep/README.md`](data-prep/README.md) for the full option list and assembly-selection rules.

## Options in detail

The defaults are good for a first run. Reach for these when you want to tune the arms.

### Shared settings

`threads`, `max_suboptimal_hits`, `max_total_energy`, and `seed` sit at the top level so both tools read the same value. Two of them carry arm-specific meaning worth knowing:

- **`max_total_energy`** — In IntaRNA this filters *total* interaction energy (hybridization plus the accessibility penalties `ED1 + ED2`); in RNAhybrid it filters *pure hybridization* energy. The two are directly comparable only when IntaRNA's accessibility correction is off (`accessibility_variant: "off"`).
- **`max_suboptimal_hits`** — `null` means "keep every hit below `max_total_energy`, per pair" on both arms. With no energy cutoff set, that collapses RNAhybrid to best-hit-only per pair. Set it to `N` to cap both arms at `N` hits per pair. (IntaRNA's hard ceiling is 1000 hits/pair; RNAhybrid is uncapped, so the two can only differ for a pair with more than 1000 suboptimal hits below the cutoff.)

`seed` is covered under [Seed handling](#seed-handling).

### Calibration (RNAhybrid arm)

`rnacalibrate.calibration_variant` selects which variant(s) the pipeline produces in a single run:

- **`on`** — runs `rnacalibrate` first and feeds each miRNA's fitted `(xi, theta)` into RNAhybrid. Outputs land under `{sample}/w_calibration/`.
- **`off`** — skips calibration; RNAhybrid falls back to `rnahybrid.distribution`, or to its built-in `species` default if that is unset. Outputs land under `{sample}/wo_calibration/`.
- **`both`** — produces both trees in one run for side-by-side comparison.

When calibration runs, `{sample}/w_calibration/rnacalibrate.json` holds one `(xi, theta)` row per query miRNA, and RNAhybrid is invoked once per `(miRNA, target chunk)` pair with that miRNA's own parameters — so every p-value reflects its own null model rather than a population average. The static `rnahybrid.distribution` value is ignored whenever calibration runs.

`rnacalibrate.randomize_targets` maps to RNAcalibrate's `-s` flag and should usually stay `false` unless you have confirmed it produces valid fits for your inputs.

### Accessibility (IntaRNA arm)

`intarna.accessibility_variant` is the IntaRNA-side analogue of `calibration_variant`:

- **`"on"`** — accessibility correction enabled (`--acc=C`). Outputs land under `{sample}/w_accessibility/`.
- **`"off"`** — no accessibility correction (`--acc=N`). Outputs land under `{sample}/wo_accessibility/`.
- **`"both"`** — produces both trees in one run.

The two variants differ **only** in `--acc=C` vs `--acc=N`; every other IntaRNA key (`prediction_mode`, `model`, `max_interaction_length`, `max_loop_size`, `seed.*`, `accessibility.*`) applies to both. The `accessibility.*` flags are still passed under `--acc=N` (where they're a no-op) so the command line stays consistent.

The top-level `threads` maps to IntaRNA's native `--threads`. Unlike the RNAhybrid arm (which chunks targets across GNU Parallel), IntaRNA runs as a single multi-threaded process per sample.

`intarna.output.columns` is the column whitelist passed to IntaRNA via `--outCsvCols`. The downstream `tidy_intarna` rule needs at least `id1, id2, start1, end1, start2, end2, E`; trimming below that breaks the arm. Keep `hybridDP` and `subseqDP` too if you want `enhance_intarna` to draw duplex visualizations.

> **YAML 1.1 gotcha:** bare `on`/`off` parse as Python `True`/`False`. The Snakefile coerces them back for **both** `accessibility_variant` and `calibration_variant`, so `variant: on` works — but the quoted `"on"`/`"off"` form is more portable.

### Seed handling

Top-level `seed: "x,y"` (in query/miRNA coordinates) is the single declaration that drives both arms:

- **RNAhybrid** receives it verbatim as `-f x,y`.
- **IntaRNA** derives `--seedBP=y-x+1`, `--seedQRange="x-y"`, and `--seedMaxUP=0` from it. The derived `--seedBP` is validated to be in `[2, 20]`.

`enabled`, `length`, and `query_range` are **inherited from the top-level `seed`** — there is no separate enable toggle and no seedBP/range keys in the config. Setting `seed: null` (the default) emits `--noSeed` on the IntaRNA side and drops `-f` on the RNAhybrid side, so neither arm enforces a seed.

The remaining `intarna.seed.*` keys are genuine, independent overrides that apply only when a seed is enforced: `max_unpaired_bases` (`--seedMaxUP`, defaults to `0` under a seed for RNAhybrid parity), `target_range` (`--seedTRange`, target coordinates), `max_energy`, `max_hybrid_energy`, `min_unpaired_probability`, `forbid_gu`, `forbid_gu_at_ends`, and `report_best_only`.

#### Recipe: RNAhybrid-emulation profile

To make IntaRNA behave like RNAhybrid (`--noSeed --acc=N --intLoopMax=30 --mode=M`) for an apples-to-apples comparison, set explicitly:

```yaml
seed: null                      # top-level: no seed on either arm
intarna:
  accessibility_variant: "off"
  prediction_mode: M
  max_loop_size: 30
```

This yields no seed constraint, no accessibility correction, exact mode, and an energy threshold comparable to the shared `max_total_energy`.

### Plots

The `plots` block is opt-in: omit it (or clear `type:`) and no plot job runs and no R environment is built. Plots currently draw only from the RNAhybrid `rnahybrid_annotated.csv` outputs; plotting the IntaRNA or consensus tables is future work.

To re-render plots after editing the block, without re-running the upstream pipeline:

```bash
snakemake --use-conda --cores N --forcerun build_plots
```

- **`type`** — which panels to render. `all` renders every panel in registered order; otherwise pass any comma-separated subset (e.g. `"pvalue_distribution,position_pvalue"`). Available panels: `pvalue_distribution`, `position_pvalue`, `position_energy`, `volcano`, `pvalue_ecdf`, `calibration_delta`, `per_mirna`, `position_density`, `seed_class`. They appear in the order listed; `calibration_delta` is auto-skipped when only one calibration variant ran.
- **`basesize`** — ggplot2 base font size.
- **`pvalue_threshold`** — cutoff for the scatter points in `position_pvalue`, `position_energy`, and `volcano`; also switches `position_density` to a two-tier (significant / non-significant) stack. The other panels ignore it.
- **`per_mirna_top_n`** — how many miRNAs to keep in the `per_mirna` panel (top-N by hit count; defaults to 12).
- **`locus`**, **`gene`**, **`protein`** — lists of validated hits to highlight. Each entry is either **bare** (`<target>` — matches any miRNA hitting that target) or **paired** (`<target>,<miRNA>` — matches only when the row's miRNA equals the named one). Bare and paired entries can mix freely in one list. Quote entries with spaces or punctuation: `"30S ribosomal protein S2,hsa-miR-X"`. The `protein:` block is optional and needs a `protein_name` column in the annotated CSV.

The `pvalue_distribution` panel draws one dashed vertical line per validated entity (lowest p-value per `(miRNA, locus_tag)`); `position_pvalue` and `position_energy` plot every alignment of the matched pairs, keyed to a shared colour and shape legend. The output filename bakes in the `type` value (commas become hyphens), so cycling through types in the same `results_dir` never overwrites earlier PDFs. With `rnacalibrate.calibration_variant: both`, panels are faceted by calibration variant; with one variant the facet collapses automatically.

## How it works

Two prediction arms run in parallel from the same inputs, and a consensus arm combines them. The Snakefile is the source of truth for the exact rule graph; the scripts it drives live in `Workflow/Scripts/`.

```text
                    ┌ (RNAcalibrate) → RNAhybrid → tidy → annotate → enhance
                    │                                         │  └────────→ plots (optional)
queries + targets ──┤                                         ├→ intersect → enhance   (consensus)
                    │                                         │
                    └ IntaRNA ─────────────────────→ tidy → annotate → enhance
```

### RNAhybrid arm

1. **rnacalibrate** *(optional)* — fits an extreme-value distribution per query miRNA, emitting one `(xi, theta)` pair each.
2. **rnahybrid** — runs in parallel via GNU Parallel. In calibrated mode each `(miRNA, target chunk)` pair uses that miRNA's own `(xi, theta)`.
3. **tidy_rnahybrid** — filters and arranges the raw output into a tidy CSV.
4. **annotate_rnahybrid** — parses FASTA headers from the target genome and merges gene metadata onto the tidy rows. Shares `Workflow/Scripts/annotate.py` with the IntaRNA arm (the `insert_after` column is parameterized per rule).
5. **enhance_rnahybrid** — appends human-readable alignment visualizations.
6. **build_plots** *(optional, RNAhybrid only)* — renders configurable ggplot2 PDFs. Skipped entirely when `plots.type` is unset.

### IntaRNA arm

7. **intarna** — runs IntaRNA with native `--threads`; no external chunking. The `w_accessibility` and `wo_accessibility` variants differ only in `--acc=C` vs `--acc=N`.
8. **tidy_intarna** — filters and renames columns; computes `Position = Start1 / Gene_length` for positional comparability with the RNAhybrid arm.
9. **annotate_intarna** — parses FASTA headers and merges metadata (shares `annotate.py`, inserting after the `E` column).
10. **enhance_intarna** — per-record human-readable report with energy, accessibility, and seed metadata plus a duplex block rendered from `subseqDP`/`hybridDP`.

### Consensus arm

11. **intersect** — inner-joins the two arms' annotated tables on `(miRNA, target)`, keeping each arm's best hit per pair (lowest `P_value` for RNAhybrid, lowest `E` for IntaRNA). Emits `consensus/consensus_annotated.csv` — the pairs both tools predict. When an arm produced two variants, its `w_*` tree is used (calibrated / accessibility-on), falling back to `wo_*` when only that one ran.
12. **enhance_consensus** — renders the consensus set as a per-record report: merged metadata plus both duplexes (RNAhybrid ASCII alignment and IntaRNA dot-bracket).

## Citation

If you use this workflow in a publication, please cite this repository and the external tools it relies on:

- **RNAhybrid / RNAcalibrate**: Rehmsmeier, M., Steffen, P., Höchsmann, M., and Giegerich, R. (2004). Fast and effective prediction of microRNA/target duplexes. *RNA*, 10(10), 1507–1517. https://doi.org/10.1261/rna.5248604
- **RNAhybrid (web server)**: Krüger, J. and Rehmsmeier, M. (2006). RNAhybrid: microRNA target prediction easy, fast and flexible. *Nucleic Acids Research*, 34(Web Server issue), W451–W454. https://doi.org/10.1093/nar/gkl243
- **IntaRNA**: Mann, M., Wright, P. R., and Backofen, R. (2017). IntaRNA 2.0: enhanced and customizable prediction of RNA–RNA interactions. *Nucleic Acids Research*, 45(W1), W435–W439. https://doi.org/10.1093/nar/gkx279
- **Snakemake**: Mölder, F., Jablonski, K. P., Letcher, B., et al. (2021). Sustainable data analysis with Snakemake. *F1000Research*, 10, 33. https://doi.org/10.12688/f1000research.29032.2
- **GNU Parallel**: Tange, O. (2011). GNU Parallel: The command-line power tool. *;login: The USENIX Magazine*, 36(1), 42–47.
