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
pip install pywebview pyqt5 pyqtwebengine pyyaml
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

### Tests

Unit tests cover the data-prep parsers and the pipeline post-processing logic
(no Snakemake or conda env needed — just `pandas` and `pyyaml` in the active env):

```bash
pytest data-prep/tests Workflow/Scripts/tests UI/tests
```

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
seed:
```

Settings both tools obey, kept in one place so the two arms can't drift apart: CPU threads, how many hits to keep per pair, and the seed constraint. The energy cutoff is **not** shared — the two tools measure energy on different scales, so each arm carries its own (see [Energy cutoffs](#energy-cutoffs)).

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

IntaRNA's parameters, including whether to apply accessibility correction, the prediction mode, the interaction model, and the Nearest-Neighbor energy table (`energy_set`: `Turner99`, `Turner04`, or `Andronescu07`), plus detailed helix / seed / accessibility / output sub-blocks. The `helix.*` keys apply **only** under `model: B` and are ignored otherwise. See [Helix block](#helix-block-intarna-arm), [Accessibility](#accessibility-intarna-arm), and [Seed handling](#seed-handling).

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

`threads`, `max_suboptimal_hits`, and `seed` sit at the top level so both tools read the same value.

- **`max_suboptimal_hits`** — `null` means "keep every hit below the arm's energy cutoff, per pair" on both arms. With no cutoff set, that collapses RNAhybrid to best-hit-only per pair. Set it to `N` to cap both arms at `N` hits per pair. (IntaRNA's hard ceiling is 1000 hits/pair; RNAhybrid is uncapped, so the two can only differ for a pair with more than 1000 suboptimal hits below the cutoff.)

`seed` is covered under [Seed handling](#seed-handling).

### Energy cutoffs

Each arm gates on **hybridization** energy, but the two tools state it on different scales, so each carries its own key:

```yaml
rnahybrid:
  max_hybrid_energy: -18      # -e ; RNAhybrid's mfe
intarna:
  max_hybrid_energy: -12.9    # IntaRNA's E_hybrid column
```

IntaRNA's `E_hybrid` includes duplex initiation (~+4.1 kcal/mol under Turner04) plus terminal-AU and dangling-end terms; RNAhybrid's `mfe` includes none of them. Measured over 1,122 matched pairs the difference is a constant **+5.10 kcal/mol** (additive, not proportional: regression slope +0.0019). So `-18` on the literature's scale is `-12.9` on IntaRNA's — the same bar, stated twice.

**This is not a loosened threshold.** −18 is the value the source papers used; +5.10 is measured; −12.9 follows by arithmetic. The pipeline previously applied `-18` to IntaRNA's **total** energy (`E_hybrid + ED1 + ED2`), which demanded roughly −33 to −36 on RNAhybrid's scale — about twice the literature threshold — and silently emptied the consensus for 3 of 4 samples. Setting both keys to the same number reinstates that bug; the pipeline warns when they drift more than 1 kcal/mol from the expected relationship.

Accessibility (`ED1`, `ED2`, `Pu1`, `Pu2`) is **reported by default, and gated only through two optional per-side floors**: `intarna.min_target_unpaired_probability` (`Pu1`, gene side) and `intarna.min_query_unpaired_probability` (`Pu2`, miRNA side), both null unless set. They are kept separate because the two sides' `Pu` distributions sit orders of magnitude apart, and they are applied in `tidy_intarna` **before** the per-pair cap — so a floor selects *which* site represents a pair rather than merely asking whether the strongest site happened to be accessible. Filtering in the other order drops pairs that do have an accessible site; measured on four genomes, the two orders differ by up to 3.5x in surviving pairs. Under `accessibility_variant: off` the floors are withheld entirely (`--acc=N` computes no `Pu`) and the pipeline warns.

Treat the floors as a deliberate choice, not a default: published `Pu ≥ 0.001` thresholds were benchmarked on bacterial sRNA near start codons, and applying them to CDS interiors deletes wet-lab-confirmed interactions. The tool-side knobs `intarna.output.min_unpaired_probability` (`--outMinPu`) and `intarna.seed.min_unpaired_probability` also remain available and null by default. `--outMinPu` is **not** the same bar as the per-side floors — it requires *every interacting position* to clear the value, on *both* sides, during prediction, where a floor bounds the whole site's `Pu1`/`Pu2` per side afterwards; setting both compounds them, which the pipeline warns about. One caveat: because `--outMaxE` is withheld on `w_accessibility`, IntaRNA's own default of `0` applies there, imposing a weak implicit `ED1 + ED2 < |E_hybrid|` bound — it only bites when accessibility energy is large, but it means a site with e.g. `E_hybrid = -13.5` and `ED1 + ED2 = 14.0` survives on `wo_accessibility` yet is silently absent on `w_accessibility`.

`intarna.accessibility_search_depth` (default 20) raises `--outNumber` on the `w_accessibility` arm only. Under `--acc=C` IntaRNA ranks by total energy while the gate is on `E_hybrid`, so the strongest duplex can sit outside a top-1 list; reporting deeper makes it reachable. It is a floor, never a ceiling — with `max_suboptimal_hits: null` it is ignored.

### The consensus join

`consensus_annotated.csv` joins the two arms on `(miRNA, Gene)` only. Each arm contributes **one** row per pair, chosen independently by its own gated energy, and coordinates never constrain the merge — so a consensus row means *both tools found a qualifying site somewhere in this gene*, not *both tools agree on this site*. `Site_offset_nt` reports how far apart the two chosen sites are, in nucleotides, **signed** — positive when IntaRNA's site lies downstream of RNAhybrid's. `max_suboptimal_hits` therefore cannot change the consensus row count as long as it stays at or below `accessibility_search_depth`; it changes only which site each arm contributes. On `w_accessibility`, raising `max_suboptimal_hits` above `accessibility_search_depth` makes IntaRNA report deeper into its total-`E` ranking, and a pair whose only qualifying site sat below the old depth can newly appear — so above that threshold the row count can change too.

### Calibration (RNAhybrid arm)

`rnacalibrate.calibration_variant` selects which variant(s) the pipeline produces in a single run:

- **`on`** — runs `rnacalibrate` first and feeds each miRNA's fitted `(xi, theta)` into RNAhybrid. Outputs land under `{sample}/w_calibration/`.
- **`off`** — skips calibration; RNAhybrid falls back to `rnahybrid.distribution`, or to its built-in `species` default if that is unset. Outputs land under `{sample}/wo_calibration/`.
- **`both`** — produces both trees in one run for side-by-side comparison.

When calibration runs, `{sample}/w_calibration/rnacalibrate.json` holds one `(xi, theta)` row per query miRNA, and RNAhybrid is invoked once per `(miRNA, target chunk)` pair with that miRNA's own parameters — so every p-value reflects its own null model rather than a population average. The static `rnahybrid.distribution` value is ignored whenever calibration runs.

`rnacalibrate.randomize_targets` maps to RNAcalibrate's `-s` flag. With `-s` (`true`, the default), p-values are calibrated against random sequences generated from each target's dinucleotide distribution — a proper null model. With `false`, the real target sequences are used directly as the random database, which is usually **not** what you want. Keep it `true` unless you have a specific reason to calibrate against the targets themselves.

`rnacalibrate.rng_seed` pins the RNG that `-s` draws from, so two runs over identical inputs (same query, same target, same params) produce byte-identical `(xi, theta)` fits — RNAcalibrate itself exposes no seed flag, so the pin is injected from outside via `libfaketime`. `null` means "don't pin": RNAcalibrate seeds from the wall clock and the fit is not reproducible run to run. The shipped `Config/config.yaml` sets `rng_seed: 1`, so calibration is pinned and reproducible by default — set it to `null` yourself to run unpinned. The seed only does something when `randomize_targets` is `true`, since `-s` is the RNG's only consumer — set alongside `randomize_targets: false` it has no effect, and the pipeline warns on stderr rather than failing.

### Accessibility (IntaRNA arm)

`intarna.accessibility_variant` is the IntaRNA-side analogue of `calibration_variant`:

- **`"on"`** — accessibility correction enabled (`--acc=C`). Outputs land under `{sample}/w_accessibility/`.
- **`"off"`** — no accessibility correction (`--acc=N`). Outputs land under `{sample}/wo_accessibility/`.
- **`"both"`** — produces both trees in one run.

The two variants differ in `--acc=C` vs `--acc=N`, and — when set — in the per-side `Pu` floors, which are withheld on `wo_accessibility` because `--acc=N` computes no unpaired probabilities. Every other IntaRNA key (`prediction_mode`, `model`, `energy_set`, `max_interaction_length`, `max_loop_size`, `helix.*`, `seed.*`, `accessibility.*`) applies to both. The `accessibility.*` flags are still passed under `--acc=N` (where they're a no-op) so the command line stays consistent.

One helix key is variant-sensitive: `helix.min_unpaired_probability` (`--helixMinPu`) is **completely inert** under `--acc=N`, but bites hard under `--acc=C` — at `0.9` it eliminated every interaction in testing. The UI greys it out when the accessibility variant is `off`.

The top-level `threads` maps to IntaRNA's native `--threads`. Unlike the RNAhybrid arm (which chunks targets across GNU Parallel), IntaRNA runs as a single multi-threaded process per sample.

`intarna.output.columns` is the column whitelist passed to IntaRNA via `--outCsvCols`. `tidy_intarna` needs at least `id1, id2, start1, end1, start2, end2, E`. The **consensus** arm additionally uses the carried metadata (`E_hybrid, ED1, ED2, Pu1, Pu2, subseqDP, hybridDP, seedStart1..seedEnd2`); trimming those does not crash consensus (they are filled `NA` with a note), but the consensus table will lack them. Keep `hybridDP` and `subseqDP` if you want `enhance_intarna` duplex visualizations.

### Helix block (IntaRNA arm)

`intarna.model: B` selects IntaRNA's **helix-block** model, which decomposes an interaction into short stable helices — each at most `helix.max_bp` base pairs — joined by flexible interior loops, instead of extending a single seed. The six `intarna.helix.*` keys tune that decomposition, and they are **ignored under every other model** (`X`, `S`, `P`): a maximally restrictive helix block leaves the predicted energy unchanged under all three. Snakemake warns at start-up if you set them under the wrong model, and the UI greys the whole group out.

Model `B` is not a cosmetic switch. On a test duplex it predicted −30.3 kcal/mol where the default `X` gave −36.4, with no helix keys set at all.

**The seed must fit inside one helix.** IntaRNA hard-errors when the derived `--seedBP` exceeds `--helixMaxBP`, so the pipeline checks it itself and fails at DAG-build time — before the RNAhybrid arm burns an hour — with a message naming the config keys rather than the IntaRNA flags. Remember that the top-level `seed: x,y` derives `seedBP = y − x + 1` (so the shipped `seed: 2,7` is **6** base pairs, not 7). Both directions are reachable: widen the seed past `helix.max_bp`, or lower `helix.max_bp` below the seed's width. If you lower `helix.max_bp`, keep it ≥ the seed width.

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

This yields no seed constraint, no accessibility correction, exact mode, and — with `intarna.max_hybrid_energy` — an energy threshold directly comparable to `rnahybrid.max_hybrid_energy`, since `E ≡ E_hybrid` when accessibility is off.

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
3. **tidy_rnahybrid** — filters and arranges the raw output into a tidy CSV. Emits `Position` as *(1-based first paired base)* `/ Gene_length`, advancing past the unpaired 5' dangling target nucleotide RNAhybrid's alignment window opens on, so the fraction anchors exactly where IntaRNA's `start1` does.
4. **annotate_rnahybrid** — parses FASTA headers from the target genome and merges gene metadata onto the tidy rows. Shares `Workflow/Scripts/annotate.py` with the IntaRNA arm (the `insert_after` column is parameterized per rule).
5. **enhance_rnahybrid** — appends human-readable alignment visualizations.
6. **build_plots** *(optional, RNAhybrid only)* — renders configurable ggplot2 PDFs. Skipped entirely when `plots.type` is unset.

### IntaRNA arm

7. **intarna** — runs IntaRNA with native `--threads`; no external chunking. The `w_accessibility` and `wo_accessibility` variants differ in `--acc=C` vs `--acc=N`, and in whether the per-side `Pu` floors apply.
8. **tidy_intarna** — filters and renames columns; is the authoritative energy gate on both IntaRNA variants (gates on `E_hybrid`, then applies the optional per-side `Pu1`/`Pu2` floors, then ranks by `E_hybrid` and caps per pair — the floors precede the cap so they choose the representative site); computes `Position = Start1 / Gene_length` for positional comparability with the RNAhybrid arm.
9. **annotate_intarna** — parses FASTA headers and merges metadata (shares `annotate.py`, inserting after the `E` column).
10. **enhance_intarna** — per-record human-readable report with energy, accessibility, and seed metadata plus a duplex block rendered from `subseqDP`/`hybridDP`.

### Consensus arm

11. **intersect** — inner-joins the two arms' annotated tables on `(miRNA, target)`, keeping each arm's representative hit per pair by its best **gated energy** (lowest `Energy` for RNAhybrid, tie-broken by `P_value`; lowest `E_hybrid` for IntaRNA, tie-broken by `E`). Emits `consensus/consensus_annotated.csv` — the pairs both tools predict — and also emits `Site_offset_nt`, the **signed** distance in nucleotides from RNAhybrid's chosen site to IntaRNA's (positive = IntaRNA downstream). When an arm produced two variants, its `w_*` tree is used (calibrated / accessibility-on), falling back to `wo_*` when only that one ran.
12. **enhance_consensus** — renders the consensus set as a per-record report: merged metadata plus both duplexes (RNAhybrid ASCII alignment and IntaRNA dot-bracket).

## Citation

If you use this workflow in a publication, please cite this repository and the external tools it relies on:

- **RNAhybrid / RNAcalibrate**: Rehmsmeier, M., Steffen, P., Höchsmann, M., and Giegerich, R. (2004). Fast and effective prediction of microRNA/target duplexes. *RNA*, 10(10), 1507–1517. https://doi.org/10.1261/rna.5248604
- **RNAhybrid (web server)**: Krüger, J. and Rehmsmeier, M. (2006). RNAhybrid: microRNA target prediction easy, fast and flexible. *Nucleic Acids Research*, 34(Web Server issue), W451–W454. https://doi.org/10.1093/nar/gkl243
- **IntaRNA**: Mann, M., Wright, P. R., and Backofen, R. (2017). IntaRNA 2.0: enhanced and customizable prediction of RNA–RNA interactions. *Nucleic Acids Research*, 45(W1), W435–W439. https://doi.org/10.1093/nar/gkx279
- **Snakemake**: Mölder, F., Jablonski, K. P., Letcher, B., et al. (2021). Sustainable data analysis with Snakemake. *F1000Research*, 10, 33. https://doi.org/10.12688/f1000research.29032.2
- **GNU Parallel**: Tange, O. (2011). GNU Parallel: The command-line power tool. *;login: The USENIX Magazine*, 36(1), 42–47.
