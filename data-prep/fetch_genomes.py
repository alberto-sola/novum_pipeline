"""
fetch_genomes.py — resolve taxon names to NCBI RefSeq reference assemblies and
download a chosen file (CDS FASTA by default) from each genome data package.

Workflow:
  1. Read taxon names (one per line) from --names
  2. GET /datasets/v2/genome/taxon/{name}/dataset_report for each name (throttled)
  3. POST /datasets/v2/genome/download with the collected accessions
  4. Extract the requested per-accession file → <out_dir>/<Name>__<accession>.<ext>
  5. Write download_report.tsv with per-taxon status (TAXON, STATUS, ACCESSION, FILE_PATH)

Environment:
  NCBI_API_KEY     Optional. Raises rate limit from ~3 to ~10 req/s.
  NCBI_API_EMAIL   Optional. Appended to the User-Agent header (NCBI etiquette).

Usage:
  python fetch_genomes.py [--names bacterial_names.txt] [--out-dir cds]
                         [--report download_report.tsv]
                         [--include CDS_FASTA [CDS_FASTA ...]]
                         [--zip-filename cds_from_genomic.fna]
                         [--batch-size 200]
"""
from __future__ import annotations

import argparse
import os
import shutil
import tempfile
import time
import zipfile
from pathlib import Path
from urllib.parse import quote

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# NCBI Datasets v2 REST API. Switch this if NCBI moves the endpoint.
_API_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"

# NCBI etiquette: identify the client. Email comes from NCBI_API_EMAIL if set.
_DEFAULT_USER_AGENT = "ncbi-cds-downloader/1.0"

# Filename inside each ncbi_dataset/data/<accession>/ directory to extract.
# CDS FASTA by default; override via --zip-filename for protein.faa, rna.fna, genomic.gff, etc.
_DEFAULT_ZIP_FILENAME = "cds_from_genomic.fna"

# HTTP status codes worth retrying with backoff (rate limits + transient server errors).
_RETRIABLE_STATUS = (429, 500, 502, 503, 504)

# Per-request delay between taxon resolution calls. NCBI: 3 req/s without a key, 10 with one.
_RATE_DELAY_NO_KEY = 0.35
_RATE_DELAY_WITH_KEY = 0.12


#----- Builds a requests session with retry-on-429/5xx and NCBI-etiquette headers -----#
def make_session(api_key: str | None = None, contact_email: str | None = None) -> requests.Session:
    s = requests.Session()
    user_agent = _DEFAULT_USER_AGENT
    if contact_email:
        user_agent = f"{user_agent} ({contact_email})"
    s.headers["User-Agent"] = user_agent
    if api_key:
        s.headers["api-key"] = api_key
    # urllib3 handles exponential backoff and respects Retry-After headers, so we don't roll our own.
    retry = Retry(total=3, backoff_factor=2.0, status_forcelist=_RETRIABLE_STATUS)
    s.mount("https://", HTTPAdapter(max_retries=retry))
    return s


#----- Resolves each taxon to its top RefSeq reference accession (None if NCBI has no match) -----#
def resolve_accessions(names: list[str], session: requests.Session, delay: float) -> dict[str, str | None]:
    result: dict[str, str | None] = {}
    for i, name in enumerate(names, 1):
        print(f"  [{i:>3}/{len(names)}] {name!r} ... ", end="", flush=True)
        # NCBI requires URL-encoded taxon names in the path (handles spaces, parentheses).
        url = f"{_API_BASE}/genome/taxon/{quote(name, safe='')}/dataset_report"
        r = session.get(url, params={
            "filters.reference_only": "true",
            "filters.assembly_source": "refseq",
            "page_size": "1",
        }, timeout=30)
        if r.status_code != 200:
            print(f"HTTP {r.status_code}")
            result[name] = None
        else:
            reports = r.json().get("reports", [])
            if reports:
                # `accession` is canonical; `current_accession` is set when the assembly was renamed.
                acc = reports[0].get("accession") or reports[0].get("current_accession")
                print(acc)
                result[name] = acc
            else:
                print("no_reference")
                result[name] = None
        # Hand-rolled rate limit; the urllib3 retry above only triggers on actual 429 responses.
        time.sleep(delay)
    return result


#----- Bulk-downloads accessions in batches via POST /genome/download, streaming each ZIP to a temp file -----#
def download_batch(accessions: list[str], include_types: list[str], session: requests.Session, batch_size: int) -> list[Path]:
    zip_paths: list[Path] = []
    n_batches = (len(accessions) + batch_size - 1) // batch_size
    for b, start in enumerate(range(0, len(accessions), batch_size), 1):
        chunk = accessions[start : start + batch_size]
        print(f"  Batch {b}/{n_batches}: {len(chunk)} accessions ...", flush=True)
        r = session.post(
            f"{_API_BASE}/genome/download",
            json={"accessions": chunk, "include_annotation_type": include_types},
            stream=True,
            timeout=600,
        )
        r.raise_for_status()
        # Stream to disk: bacterial CDS sets are small, but multi-include or eukaryote payloads can be GB-scale.
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".zip")
        try:
            for data in r.iter_content(chunk_size=1 << 20):  # 1 MiB
                tmp.write(data)
        finally:
            tmp.close()
        zip_paths.append(Path(tmp.name))
    return zip_paths


#----- Walks each ZIP and copies the requested filename out of every accession dir to <Name>__<acc>.<ext> -----#
def extract_flat(zip_paths: list[Path], acc_to_name: dict[str, str], out_dir: Path, zip_filename: str = _DEFAULT_ZIP_FILENAME) -> dict[str, str | None]:
    out_dir.mkdir(parents=True, exist_ok=True)
    # Pre-fill with None so we can flag accessions present in the request but missing in the ZIP payload.
    result: dict[str, str | None] = {acc: None for acc in acc_to_name}
    out_suffix = Path(zip_filename).suffix or ".dat"

    for zip_path in zip_paths:
        try:
            with zipfile.ZipFile(zip_path) as zf:
                for member in zf.namelist():
                    parts = member.split("/")
                    # NCBI Datasets ZIP layout: ncbi_dataset/data/<accession>/<filename>
                    if (
                        len(parts) == 4
                        and parts[0] == "ncbi_dataset"
                        and parts[1] == "data"
                        and parts[3] == zip_filename
                    ):
                        acc = parts[2]
                        name = acc_to_name.get(acc)
                        if name is None:
                            continue
                        # Spaces in scientific names break shell scripts; underscore them in filenames.
                        safe_name = name.replace(" ", "_")
                        dest = out_dir / f"{safe_name}__{acc}{out_suffix}"
                        with zf.open(member) as src, dest.open("wb") as dst:
                            shutil.copyfileobj(src, dst)
                        result[acc] = str(dest)
                        print(f"    {dest.name}")
        finally:
            zip_path.unlink(missing_ok=True)

    return result


#----- Writes the per-taxon TSV audit report (one row per input name; uppercase column headers) -----#
def write_report(path: Path, names: list[str], resolved: dict[str, str | None], extraction: dict[str, str | None]) -> None:
    with path.open("w") as fh:
        fh.write("TAXON\tSTATUS\tACCESSION\tFILE_PATH\n")
        for name in names:
            acc = resolved.get(name)
            if acc is None:
                fh.write(f"{name}\tno_reference\t-\t-\n")
                continue
            dest = extraction.get(acc)
            if dest is None:
                fh.write(f"{name}\tmissing_in_zip\t{acc}\t-\n")
            else:
                fh.write(f"{name}\tok\t{acc}\t{dest}\n")


#----- CLI entry point: parse args, resolve, download, extract, report -----#
def main() -> None:
    ap = argparse.ArgumentParser(description="Download per-taxon files from NCBI Datasets by taxon name")
    ap.add_argument("--names", default="bacterial_names.txt",
                    help="Input: one taxon name per line (default: bacterial_names.txt)")
    ap.add_argument("--out-dir", default="genomes",
                    help="Output directory for extracted files (default: genomes/)")
    ap.add_argument("--report", default="download_report.tsv",
                    help="TSV audit report path (default: download_report.tsv)")
    ap.add_argument("--include", nargs="+", default=["CDS_FASTA"], metavar="TYPE",
                    help="Annotation types to download (default: CDS_FASTA). "
                         "Common values: GENOME_FASTA, RNA_FASTA, CDS_FASTA, PROT_FASTA, GENOME_GFF.")
    ap.add_argument("--zip-filename", default=_DEFAULT_ZIP_FILENAME,
                    help=f"File to extract from each accession dir (default: {_DEFAULT_ZIP_FILENAME}). "
                         "Examples: protein.faa, rna.fna, genomic.gff.")
    ap.add_argument("--batch-size", type=int, default=200,
                    help="Max accessions per POST request (default: 200)")
    args = ap.parse_args()

    base = Path(__file__).parent
    names_path = base / args.names
    out_dir = base / args.out_dir
    report_path = base / args.report

    api_key = os.environ.get("NCBI_API_KEY")
    contact_email = os.environ.get("NCBI_API_EMAIL")
    delay = _RATE_DELAY_WITH_KEY if api_key else _RATE_DELAY_NO_KEY

    names = [ln.strip() for ln in names_path.read_text().splitlines() if ln.strip()]
    print(f"Loaded {len(names)} names from {names_path}\n")

    session = make_session(api_key=api_key, contact_email=contact_email)

    print("--- Resolving reference accessions ---")
    resolved = resolve_accessions(names, session, delay)

    acc_to_name = {acc: name for name, acc in resolved.items() if acc is not None}
    no_ref = [n for n, a in resolved.items() if a is None]
    print(f"\nResolved: {len(acc_to_name)}  No reference: {len(no_ref)}")
    if no_ref:
        print("  Skipped:", ", ".join(no_ref))

    extraction: dict[str, str | None] = {}
    if acc_to_name:
        print("\n--- Downloading packages ---")
        zip_paths = download_batch(list(acc_to_name), args.include, session, args.batch_size)
        print(f"\n--- Extracting {args.zip_filename} ---")
        extraction = extract_flat(zip_paths, acc_to_name, out_dir, zip_filename=args.zip_filename)
        n_ok = sum(1 for v in extraction.values() if v is not None)
        print(f"\nExtracted {n_ok}/{len(acc_to_name)} files → {out_dir}/")

    write_report(report_path, names, resolved, extraction)
    print(f"Report  → {report_path}")


if __name__ == "__main__":
    main()
