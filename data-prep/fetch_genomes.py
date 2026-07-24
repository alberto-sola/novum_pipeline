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
import gzip
import os
import shutil
import tempfile
import time
import zipfile
import zlib
from pathlib import Path
from urllib.parse import quote

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# NCBI Datasets v2 REST API. Switch this if NCBI moves the endpoint.
_API_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"

# NCBI etiquette: identify the client. Email comes from NCBI_API_EMAIL if set.
_DEFAULT_USER_AGENT = "ncbi-cds-downloader/1.0"

# The CDS set NCBI ships inside each ncbi_dataset/data/<accession>/ directory.
_CDS_FILENAME = "cds_from_genomic.fna"

# Filename the standalone CLI extracts; override via --zip-filename for
# protein.faa, rna.fna, genomic.gff, etc.
_DEFAULT_ZIP_FILENAME = _CDS_FILENAME

# Datasets v2 serves no rna_from_genomic.fna for prokaryotes, so the RNA half comes
# from the FTP mirror, at a URL derived from (acc, asm) rather than queried.
_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
_RNA_SUFFIX = "rna_from_genomic.fna.gz"

# NCBI substitutes these characters in the assembly-name part of its paths.
_NAME_UNSAFE = str.maketrans({c: "_" for c in " #()/"})

# HTTP status codes worth retrying with backoff (rate limits + transient server errors).
_RETRIABLE_STATUS = (429, 500, 502, 503, 504)

# Per-request delay between taxon resolution calls. NCBI: 3 req/s without a key, 10 with one.
_RATE_DELAY_NO_KEY = 0.35
_RATE_DELAY_WITH_KEY = 0.12

# Assembly-selection ranking: lower tuple sorts first.
_CATEGORY_RANK = {"reference genome": 0, "representative genome": 1, "": 2}
_LEVEL_RANK = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}


#----- Picks the best assembly: category > completeness > newest release date -----#
def select_best(candidates: list[dict]) -> tuple[str, str] | None:
    if not candidates:
        return None

    def info(rep):
        return rep.get("assembly_info", {}) or {}

    # Newest-first pre-sort, then pick the best (category, level) with a stable
    # min — ties keep the newest because the pre-sort already ordered by date desc.
    ranked = sorted(candidates, key=lambda r: info(r).get("release_date", ""), reverse=True)
    best = min(
        ranked,
        key=lambda r: (
            _CATEGORY_RANK.get(info(r).get("refseq_category", "") or "", 3),
            _LEVEL_RANK.get(info(r).get("assembly_level", "") or "", 4),
        ),
    )
    acc = best.get("accession") or best.get("current_accession")
    if not acc:
        return None
    return acc, info(best).get("assembly_name", "")


#----- Resolves one taxon: reference first, else best of a candidate page -----#
def _resolve_one(name: str, session: requests.Session) -> tuple[str, str] | None:
    base = f"{_API_BASE}/genome/taxon/{quote(name, safe='')}/dataset_report"

    def reports_for(params):
        r = session.get(base, params=params, timeout=30)
        return r.json().get("reports", []) if r.status_code == 200 else None

    # Tier 1: a designated reference assembly (clean filter, one row).
    rep = reports_for({"filters.reference_only": "true",
                       "filters.assembly_source": "refseq", "page_size": "1"})
    if rep:
        return select_best(rep)
    # Tier 2: no reference — rank a page of RefSeq candidates client-side.
    rep = reports_for({"filters.assembly_source": "refseq", "page_size": "50"})
    return select_best(rep) if rep is not None else None


#----- Resolves each taxon to (accession, assembly_name) with reference->latest fallback -----#
def resolve_assemblies(names: list[str], session: requests.Session, delay: float) -> dict[str, tuple[str, str] | None]:
    result: dict[str, tuple[str, str] | None] = {}
    for i, name in enumerate(names, 1):
        print(f"  [{i:>3}/{len(names)}] {name!r} ... ", end="", flush=True)
        try:
            hit = _resolve_one(name, session)
        except requests.RequestException as exc:
            print(f"error ({exc.__class__.__name__})")
            hit = None
        else:
            print(f"{hit[0]} ({hit[1]})" if hit else "no_assembly")
        result[name] = hit
        time.sleep(delay)
    return result


#----- "<accession>_<assembly name>", the stem NCBI uses in both paths and filenames -----#
def assembly_stem(acc: str, asm: str) -> str:
    return f"{acc}_{asm.translate(_NAME_UNSAFE)}"


#----- FTP URL of an assembly's rna_from_genomic.fna.gz, derived from (acc, asm) -----#
def rna_ftp_url(acc: str, asm: str) -> str:
    # GCF_000005845.2 -> prefix "GCF", digits "000005845" -> .../GCF/000/005/845/
    prefix, _, rest = acc.partition("_")
    digits = rest.split(".", 1)[0]
    stem = assembly_stem(acc, asm)
    return (f"{_FTP_BASE}/{prefix}/{digits[0:3]}/{digits[3:6]}/{digits[6:9]}"
            f"/{stem}/{stem}_{_RNA_SUFFIX}")


#----- On-disk name for a target FASTA; `with_rna` picks the merged vs CDS-only name -----#
def target_dest(out_dir, acc: str, asm: str, with_rna: bool = True) -> Path:
    kind = "cds_rna" if with_rna else "cds"
    return Path(out_dir) / f"{assembly_stem(acc, asm)}_{kind}_from_genomic.fna"


#----- Appends a gzipped FASTA onto `part_path`, healing a missing newline; returns records added -----#
def append_gz_fasta(part_path, gz_path) -> int:
    part_path = Path(part_path)
    records = 0
    with part_path.open("rb+") as dst:
        dst.seek(0, os.SEEK_END)
        rollback_to = dst.tell()  # pristine CDS-only length; restore here on any failure
        # A CDS chunk missing its trailing newline would glue onto the first RNA header.
        if rollback_to:
            dst.seek(-1, os.SEEK_END)
            if dst.read(1) != b"\n":
                dst.write(b"\n")
        try:
            with gzip.open(gz_path, "rb") as src:
                for line in src:
                    if line.startswith(b">"):
                        records += 1
                    dst.write(line)
        except Exception:
            # A gzip failing mid-stream has already written partial RNA; roll back so
            # the caller's *_cds_from_genomic.fna name stays honest.
            dst.seek(rollback_to)
            dst.truncate()
            raise
    return records


#----- Yields (accession, open member) for each matching file in the ZIPs; deletes each ZIP when done -----#
def iter_package_members(zip_paths, filename):
    for zip_path in zip_paths:
        try:
            with zipfile.ZipFile(zip_path) as zf:
                for member in zf.namelist():
                    parts = member.split("/")
                    # NCBI Datasets ZIP layout: ncbi_dataset/data/<accession>/<filename>
                    if (len(parts) == 4 and parts[0] == "ncbi_dataset"
                            and parts[1] == "data" and parts[3] == filename):
                        with zf.open(member) as src:
                            yield parts[2], src
        finally:
            zip_path.unlink(missing_ok=True)


#----- Extracts each accession's cds_from_genomic.fna from the batched ZIPs into a <stem>.part file -----#
def _download_cds_parts(accessions, asm_by_acc, out_dir, session, batch_size):
    parts: dict[str, Path] = {}
    zip_paths = download_batch(accessions, ["CDS_FASTA"], session, batch_size)
    for acc, src in iter_package_members(zip_paths, _CDS_FILENAME):
        asm = asm_by_acc.get(acc)
        if asm is None:
            continue
        # Staged as .part so an interrupted run can't leave a half-merged
        # file the next run mistakes for complete.
        part = out_dir / f"{assembly_stem(acc, asm)}.part"
        with part.open("wb") as dst:
            shutil.copyfileobj(src, dst)
        parts[acc] = part
    return parts


#----- Streams a response body to a temp file the caller owns; self-cleans a failed transfer -----#
def _stream_to_tempfile(response, suffix: str) -> Path:
    with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as fh:
        path = Path(fh.name)
        try:
            for chunk in response.iter_content(chunk_size=1 << 20):  # 1 MiB
                fh.write(chunk)
        except Exception:
            path.unlink(missing_ok=True)
            raise
    return path


#----- Fetches the assembly's RNA FASTA from FTP and appends it to `part`; returns merged/missing/error -----#
def _append_rna(part, acc, asm, session):
    url = rna_ftp_url(acc, asm)
    tmp = None
    try:
        # The api-key header is scoped to the API host; strip it for the FTP mirror.
        r = session.get(url, headers={"api-key": None}, stream=True, timeout=300)
        if r.status_code == 404:
            print(f"    no rna_from_genomic for {acc}")
            return "missing"
        r.raise_for_status()
        tmp = _stream_to_tempfile(r, ".gz")
        append_gz_fasta(part, tmp)
        return "merged"
    except (requests.RequestException, OSError, EOFError, zlib.error) as exc:
        # A corrupt gzip transfers cleanly and fails only in the reader, raising
        # EOFError/zlib.error — neither subclasses OSError, so both need naming.
        # The URL is printed so a derivation bug can't look like "NCBI has no RNA".
        print(f"    RNA fetch failed for {acc}: {exc.__class__.__name__} — {url}")
        return "error"
    finally:
        if tmp is not None:
            tmp.unlink(missing_ok=True)


#----- Merges RNA onto a staged .part, renames it to its final name, drops a stale sibling -----#
def _finalize_target(part, acc, asm, out_dir, session, delay):
    status = _append_rna(part, acc, asm, session)
    merged = status == "merged"
    dest = target_dest(out_dir, acc, asm, with_rna=merged)
    part.replace(dest)
    if merged:
        # Drop any pre-migration CDS-only file so one per accession survives.
        # One-directional by design: a degraded result never deletes a merged
        # sibling, since stale merged data beats it and self-heals via the cache.
        target_dest(out_dir, acc, asm, with_rna=False).unlink(missing_ok=True)
        print(f"    {dest.name}")
    time.sleep(delay)
    return status, dest


#----- Downloads CDS (API) + rna_from_genomic (FTP) per accession, merged into one target FASTA -----#
def download_targets(hits: dict[str, tuple[str, str] | None], out_dir, session,
                     delay: float, force: bool = False,
                     batch_size: int = 200) -> tuple[dict[str, str | None], dict[str, str]]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Dedup by accession: several taxa can resolve to the same assembly, so we
    # download each accession once and map the file back to every taxon below.
    asm_by_acc: dict[str, str] = {}
    for hit in hits.values():
        if hit is not None:
            asm_by_acc.setdefault(hit[0], hit[1])

    # Only the merged name counts as a cache hit. A CDS-only file is a degraded
    # result, so we retry it — which also migrates pre-merge trees automatically.
    final_by_acc: dict[str, Path] = {}
    status_by_acc: dict[str, str] = {}
    todo: list[str] = []
    for acc, asm in asm_by_acc.items():
        merged = target_dest(out_dir, acc, asm, with_rna=True)
        if not force and merged.exists():
            final_by_acc[acc], status_by_acc[acc] = merged, "cached"
        else:
            todo.append(acc)

    if todo:
        parts = _download_cds_parts(todo, asm_by_acc, out_dir, session, batch_size)
        for acc in todo:
            part = parts.get(acc)
            if part is None:
                # CDS never landed; there is nothing to merge onto and RNA is
                # never the sole content of a target FASTA.
                status_by_acc[acc] = "-"
                continue
            status_by_acc[acc], final_by_acc[acc] = _finalize_target(
                part, acc, asm_by_acc[acc], out_dir, session, delay)

    # Map every resolved taxon to its accession's file (None if it never landed),
    # so taxa sharing one accession all receive the same path. final_by_acc holds
    # only confirmed successes — a cache hit or a completed rename.
    paths: dict[str, str | None] = {}
    rna_status: dict[str, str] = {}
    for taxon, hit in hits.items():
        if hit is None:
            continue
        dest = final_by_acc.get(hit[0])
        paths[taxon] = str(dest) if dest is not None else None
        rna_status[taxon] = status_by_acc.get(hit[0], "-")
    return paths, rna_status


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
        zip_paths.append(_stream_to_tempfile(r, ".zip"))
    return zip_paths


#----- Walks each ZIP and copies the requested filename out of every accession dir to <Name>__<acc>.<ext> -----#
def extract_flat(zip_paths: list[Path], acc_to_name: dict[str, str], out_dir: Path, zip_filename: str = _DEFAULT_ZIP_FILENAME) -> dict[str, str | None]:
    out_dir.mkdir(parents=True, exist_ok=True)
    # Pre-fill with None so we can flag accessions present in the request but missing in the ZIP payload.
    result: dict[str, str | None] = {acc: None for acc in acc_to_name}
    out_suffix = Path(zip_filename).suffix or ".dat"

    for acc, src in iter_package_members(zip_paths, zip_filename):
        name = acc_to_name.get(acc)
        if name is None:
            continue
        # Spaces in scientific names break shell scripts; underscore them in filenames.
        safe_name = name.replace(" ", "_")
        dest = out_dir / f"{safe_name}__{acc}{out_suffix}"
        with dest.open("wb") as dst:
            shutil.copyfileobj(src, dst)
        result[acc] = str(dest)
        print(f"    {dest.name}")

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
