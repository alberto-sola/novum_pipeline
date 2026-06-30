# data-prep — pipeline input bootstrap

Resolve + download bacterial CDS assemblies (NCBI) and build miRNA query FASTAs
(miRBase), then surgically write the `queries:`/`targets:` blocks of
`Config/config.yaml` (every other param and comment is preserved). It does **not**
run Snakemake — you keep manual control.

## Usage

    python data-prep/prepare_inputs.py --pairs pairs.txt              # full run
    python data-prep/prepare_inputs.py --pairs pairs.txt --dry-run    # preview, no writes/downloads of CDS
    python data-prep/prepare_inputs.py Veillonella_parvula:hsa-miR-200b-3p   # inline pair(s)

Input lines are **taxon-first**; the separator is a tab, a `:`, or spaces, and
case is flexible:

    Veillonella_parvula     hsa-miR-200b-3p
    Dysosmobacter_welbionis:hsa-mir-200c-3p
    fusobacterium_nucleatum hsa-mir-515-5p

The same taxon on multiple lines collects its miRNAs into one query FASTA.
A bare miRNA name (no species prefix) is tried as `hsa-`, then `mmu-`
(`--prefix-chain` to change). Unresolved miRNAs / taxa are skipped and recorded
in `prep_report.tsv` rather than aborting the run.

## Options

    --config PATH        config to edit          [Config/config.yaml]
    --rnas-dir DIR       query FASTA output       [Data/Raw/RNAs]
    --genomes-dir DIR    genome output            [Data/Raw/genomes]
    --mirnas-dir DIR     miRBase mature.fa cache  [Data/Raw/miRNAs]
    --prefix-chain CSV   species fallback order   [hsa,mmu]
    --report PATH        per-(taxon,miRNA) TSV    [data-prep/prep_report.tsv]
    --force              re-download even if outputs exist
    --dry-run            resolve + print the blocks; no CDS download, no config edit

## Notes

- A taxon reaches `config.yaml` only if its genome resolved **and** ≥1 miRNA
  resolved (keeps the written config runnable).
- Assembly selection prefers a RefSeq *reference* genome, then *representative*,
  then the latest/most-complete assembly — so long-tail taxa with no designated
  reference still resolve.
- `NCBI_API_KEY` / `NCBI_API_EMAIL` (optional) raise the NCBI rate limit.

See `docs/superpowers/specs/2026-06-30-data-prep-design.md` for the full design.
