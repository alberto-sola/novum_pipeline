"""prepare_inputs.py — bootstrap pipeline inputs from (taxon, miRNA) pairs:
resolve+download genomes, build miRNA query FASTAs, write Config/config.yaml."""
from __future__ import annotations

import argparse
from collections import Counter
import os
import sys
from pathlib import Path

import pairs as pairs_mod
from fetch_genomes import (make_session, resolve_assemblies, download_targets, target_dest,
                           _RATE_DELAY_NO_KEY, _RATE_DELAY_WITH_KEY)
from mirnas import (ensure_mature_fa, index_mature, resolve_mirna, write_query_fasta, query_fasta_path, MatureUnavailable)
from config_writer import update_path_blocks


#----- Emit a taxon only if its genome resolved AND >=1 miRNA resolved -----#
def select_emitted(taxa_order, genome_paths, mirna_resolved, rnas_dir):
    queries: dict[str, str] = {}
    targets: dict[str, str] = {}
    skipped: dict[str, str] = {}
    for taxon in taxa_order:
        gpath = genome_paths.get(taxon)
        mirnas = mirna_resolved.get(taxon) or []
        if gpath is None:
            skipped[taxon] = "no_assembly"
            continue
        if not mirnas:
            skipped[taxon] = "no_mirnas"
            continue
        queries[taxon] = str(query_fasta_path(rnas_dir, taxon))
        targets[taxon] = gpath
    return queries, targets, skipped


#----- Writes the per-(taxon,miRNA) TSV audit report; RNA_STATUS follows GENOME_PATH -----#
def _write_report(path, groups, genome_paths, hits_by_taxon, per_mirna, rna_status):
    with open(path, "w") as fh:
        fh.write("TAXON\tGENOME_STATUS\tACCESSION\tASSEMBLY_NAME\tGENOME_PATH\tRNA_STATUS\t"
                 "MIRNA_INPUT\tMIRNA_STATUS\tMIRNA_RESOLVED\n")
        for taxon in groups:
            hit = hits_by_taxon.get(taxon)
            gpath = genome_paths.get(taxon)
            if hit is None:
                gstatus, acc, asm = "no_assembly", "-", "-"
            elif gpath is None:
                gstatus, acc, asm = "download_failed", hit[0], hit[1]
            else:
                gstatus, acc, asm = "ok", hit[0], hit[1]
            for raw in groups[taxon]:
                status, resolved = per_mirna[(taxon, raw)]
                fh.write(f"{taxon}\t{gstatus}\t{acc}\t{asm}\t{gpath or '-'}\t"
                         f"{rna_status.get(taxon, '-')}\t{raw}\t{status}\t{resolved}\n")


#----- Bootstrap: resolve+download targets, build query FASTAs, write config and report -----#
def run(pairs_path=None, inline=None, config="Config/config.yaml",
        rnas_dir="Data/Raw/RNAs", genomes_dir="Data/Raw/genomes",
        mirnas_dir="Data/Raw/miRNAs", prefix_chain=("hsa", "mmu"),
        report="data-prep/prep_report.tsv", force=False, dry_run=False):
    text = ""
    if pairs_path:
        text += Path(pairs_path).read_text()
    if inline:
        text += ("\n" if text else "") + "\n".join(inline)
    groups, malformed = pairs_mod.parse_pairs(text)
    for lineno, raw in malformed:
        print(f"  ! malformed line {lineno}: {raw!r}", file=sys.stderr)
    if not groups:
        print("No valid pairs. Nothing to do.")
        return

    api_key = os.environ.get("NCBI_API_KEY")
    delay = _RATE_DELAY_WITH_KEY if api_key else _RATE_DELAY_NO_KEY
    session = make_session(api_key=api_key, contact_email=os.environ.get("NCBI_API_EMAIL"))

    # --- genomes ---
    name_to_key = {pairs_mod.taxon_to_query(k): k for k in groups}
    print("--- Resolving assemblies ---")
    hits_by_name = resolve_assemblies(list(name_to_key), session, delay)
    hits_by_taxon = {name_to_key[n]: h for n, h in hits_by_name.items()}
    if dry_run:
        # Optimistic: we can't know whether an assembly has an RNA file without
        # fetching it, so predict the merged name.
        genome_paths = {k: (str(target_dest(genomes_dir, h[0], h[1])) if h else None)
                        for k, h in hits_by_taxon.items()}
        rna_status = {}
    else:
        print("--- Downloading CDS + RNA ---")
        paths_by_name, rna_by_name = download_targets(hits_by_name, genomes_dir, session,
                                                      delay, force=force)
        genome_paths = {name_to_key[n]: p for n, p in paths_by_name.items()}
        rna_status = {name_to_key[n]: s for n, s in rna_by_name.items()}
        counts = Counter(s for s in rna_status.values() if s != "-")
        print("RNA: " + (", ".join(f"{n} {s}" for s, n in counts.most_common()) or "none"))

    # --- miRNAs ---
    print("--- Resolving miRNAs ---")
    # Spec: dry-run "ensures mature.fa (idempotent cache)" — download if absent
    # but never force a re-download during a preview.
    mature = ensure_mature_fa(mirnas_dir, session, force=force and not dry_run)
    index = index_mature(mature)
    mirna_resolved: dict[str, list[tuple[str, str]]] = {}
    per_mirna: dict[tuple[str, str], tuple[str, str]] = {}
    for taxon, raws in groups.items():
        resolved = []
        for raw in raws:
            hit = resolve_mirna(raw, index, prefix_chain=tuple(prefix_chain))
            if hit is None:
                per_mirna[(taxon, raw)] = ("not_found", "-")
            else:
                canonical, seq, used = hit
                per_mirna[(taxon, raw)] = ("resolved_via_fallback" if used else "ok", canonical)
                resolved.append((canonical, seq))
        mirna_resolved[taxon] = resolved

    # --- emit ---
    queries, targets, skipped = select_emitted(list(groups), genome_paths, mirna_resolved, rnas_dir)

    if dry_run:
        print("\n# queries:")
        for k, v in queries.items():
            print(f"  {k}: {v}")
        print("# targets:")
        for k, v in targets.items():
            print(f"  {k}: {v}")
        print("  (names assume rna_from_genomic is available; an assembly "
              "without it lands as *_cds_from_genomic.fna)")
    else:
        for taxon in queries:
            write_query_fasta(taxon, mirna_resolved[taxon], rnas_dir)
        update_path_blocks(config, queries, targets)
        print(f"\nWrote {len(queries)} taxa to {config}")

    if skipped:
        print(f"Skipped {len(skipped)}: " + ", ".join(f"{k} ({v})" for k, v in skipped.items()))
    else:
        print("Skipped 0")
    # Dry-run is a preview: don't write (and clobber) the authoritative report with
    # genome paths whose files were never fetched.
    if not dry_run:
        _write_report(report, groups, genome_paths, hits_by_taxon, per_mirna, rna_status)
        print(f"Report → {report}")


#----- CLI entry point: parse args, run, surface miRBase outages without a traceback -----#
def main() -> None:
    ap = argparse.ArgumentParser(description="Bootstrap pipeline inputs from (taxon, miRNA) pairs.")
    ap.add_argument("inline", nargs="*", help='inline "taxon:miRNA" pairs')
    ap.add_argument("--pairs", help="pairs file (taxon-first syntax)")
    ap.add_argument("--config", default="Config/config.yaml")
    ap.add_argument("--rnas-dir", default="Data/Raw/RNAs")
    ap.add_argument("--genomes-dir", default="Data/Raw/genomes")
    ap.add_argument("--mirnas-dir", default="Data/Raw/miRNAs")
    ap.add_argument("--prefix-chain", default="hsa,mmu")
    ap.add_argument("--report", default="data-prep/prep_report.tsv")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    if not args.pairs and not args.inline:
        ap.error("provide a --pairs FILE or inline pairs")
    try:
        run(pairs_path=args.pairs, inline=args.inline, config=args.config,
            rnas_dir=args.rnas_dir, genomes_dir=args.genomes_dir, mirnas_dir=args.mirnas_dir,
            prefix_chain=tuple(args.prefix_chain.split(",")), report=args.report,
            force=args.force, dry_run=args.dry_run)
    except MatureUnavailable as e:
        # External miRBase outage, not a bug — report cleanly, no traceback.
        print(f"\nERROR: {e}", file=sys.stderr)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
