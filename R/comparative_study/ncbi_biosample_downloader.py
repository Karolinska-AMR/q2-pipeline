#!/usr/bin/env python3
"""
ncbi_biosample_downloader.py
============================
Downloads BioSample metadata and associated sequencing run files from NCBI/ENA.

Given a CSV of BioSample accessions (e.g. SAMN*), this script:
  1. Fetches rich metadata for each BioSample via NCBI E-utilities.
  2. Finds linked SRA experiments/runs via NCBI elink.
  3. Enriches run metadata (platform, layout, library strategy) from NCBI SRA.
  4. Fetches ENA FTP/HTTPS FASTQ download links via the ENA Portal API.
  5. Downloads FASTQ files (with resume, retry, and skip logic).
  6. Writes a consolidated metadata CSV and a download summary.

Usage:
    python ncbi_biosample_downloader.py \
        --input_csv biosamples.csv \
        --id_column biosample_id \
        --outdir ./output \
        --email you@example.com \
        [--api_key XXXXXX] \
        [--metadata_only] \
        [--threads 4]

Requirements:
    pip install requests biopython tqdm pandas

Author: Generated script — adapt freely.
"""

import argparse
import csv
import json
import logging
import os
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests
import pandas as pd
from tqdm import tqdm

# ── Optional Biopython (used for Entrez helpers) ──────────────────────────────
try:
    from Bio import Entrez
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False

# ── Constants ─────────────────────────────────────────────────────────────────
NCBI_BASE   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENA_PORTAL  = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_SEARCH  = "https://www.ebi.ac.uk/ena/portal/api/search"
RETRY_LIMIT = 5           # max HTTP retries per request
RETRY_DELAY = 2           # seconds between retries (doubles on each failure)
RATE_LIMIT  = 0.34        # seconds between NCBI calls (≤3/s without API key)
RATE_LIMIT_KEY = 0.1      # seconds between NCBI calls with API key (≤10/s)
CHUNK_SIZE  = 1024 * 1024 # 1 MiB for streaming downloads

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════════════
# 1.  HTTP HELPERS
# ═══════════════════════════════════════════════════════════════════════════════

def http_get(url: str, params: dict = None, timeout: int = 60) -> requests.Response:
    """GET with exponential-backoff retry logic."""
    delay = RETRY_DELAY
    for attempt in range(1, RETRY_LIMIT + 1):
        try:
            resp = requests.get(url, params=params, timeout=timeout)
            if resp.status_code == 429:          # rate-limited
                wait = int(resp.headers.get("Retry-After", delay))
                log.warning("Rate-limited; sleeping %ds", wait)
                time.sleep(wait)
                continue
            resp.raise_for_status()
            return resp
        except requests.RequestException as exc:
            if attempt == RETRY_LIMIT:
                raise
            log.warning("Attempt %d/%d failed (%s); retrying in %ds",
                        attempt, RETRY_LIMIT, exc, delay)
            time.sleep(delay)
            delay *= 2
    raise RuntimeError(f"All {RETRY_LIMIT} attempts failed for {url}")


# ═══════════════════════════════════════════════════════════════════════════════
# 2.  NCBI E-UTILITIES WRAPPERS
# ═══════════════════════════════════════════════════════════════════════════════

class NCBIClient:
    """Thin wrapper around NCBI E-utilities REST endpoints."""

    def __init__(self, email: str, api_key: Optional[str] = None):
        self.email   = email
        self.api_key = api_key
        self.delay   = RATE_LIMIT_KEY if api_key else RATE_LIMIT
        self._last   = 0.0

        # Also configure Biopython Entrez if available
        if HAS_BIOPYTHON:
            Entrez.email   = email
            Entrez.api_key = api_key or ""
            Entrez.tool    = "ncbi_biosample_downloader"

    def _base_params(self) -> dict:
        p: dict = {"email": self.email, "tool": "ncbi_biosample_downloader"}
        if self.api_key:
            p["api_key"] = self.api_key
        return p

    def _throttle(self):
        """Politely sleep to stay within NCBI rate limits."""
        elapsed = time.time() - self._last
        if elapsed < self.delay:
            time.sleep(self.delay - elapsed)
        self._last = time.time()

    # ── esearch ───────────────────────────────────────────────────────────────
    def esearch(self, db: str, term: str) -> list[str]:
        """Return a list of UIDs matching *term* in *db*."""
        self._throttle()
        params = {**self._base_params(),
                  "db": db, "term": term,
                  "retmax": 500, "retmode": "json"}
        try:
            data = http_get(f"{NCBI_BASE}/esearch.fcgi", params=params).json()
            return data.get("esearchresult", {}).get("idlist", [])
        except Exception as exc:
            log.warning("esearch failed for '%s' in %s (%s); returning empty list",
                        term, db, exc)
            return []

    # ── efetch ────────────────────────────────────────────────────────────────
    def efetch(self, db: str, ids: list[str], rettype: str = "xml",
               retmode: str = "xml") -> str:
        """Fetch records for a list of UIDs; returns raw text."""
        self._throttle()
        params = {**self._base_params(),
                  "db": db, "id": ",".join(ids),
                  "rettype": rettype, "retmode": retmode}
        return http_get(f"{NCBI_BASE}/efetch.fcgi", params=params).text

    # ── elink ─────────────────────────────────────────────────────────────────
    def elink(self, dbfrom: str, db: str, ids: list[str]) -> list[str]:
        """Return linked UIDs from *dbfrom* → *db*."""
        self._throttle()
        params = {**self._base_params(),
                  "dbfrom": dbfrom, "db": db,
                  "id": ",".join(ids), "retmode": "json"}
        try:
            resp = http_get(f"{NCBI_BASE}/elink.fcgi", params=params)
            data = resp.json()
        except Exception as exc:
            log.warning("elink JSON parse failed (%s); returning empty list", exc)
            return []
        
        uids: list[str] = []
        for linkset in data.get("linksets", []):
            for lsd in linkset.get("linksetdbs", []):
                uids.extend(lsd.get("links", []))
        return list(set(uids))

    # ── esummary ──────────────────────────────────────────────────────────────
    def esummary(self, db: str, ids: list[str]) -> dict:
        """Return esummary JSON for a list of UIDs."""
        self._throttle()
        params = {**self._base_params(),
                  "db": db, "id": ",".join(ids), "retmode": "json"}
        try:
            resp = http_get(f"{NCBI_BASE}/esummary.fcgi", params=params)
            data = resp.json()
            return data.get("result", {})
        except Exception as exc:
            log.warning("esummary JSON parse failed for %d UIDs (%s); returning empty dict",
                        len(ids), exc)
            return {}


# ═══════════════════════════════════════════════════════════════════════════════
# 3.  BIOSAMPLE METADATA PARSING
# ═══════════════════════════════════════════════════════════════════════════════

def _xml_text(xml: str, tag: str, default: str = "") -> str:
    """Naive regex-based XML text extraction (no lxml dependency required)."""
    m = re.search(rf"<{tag}[^>]*>(.*?)</{tag}>", xml, re.DOTALL | re.IGNORECASE)
    return m.group(1).strip() if m else default


def _attr_value(xml: str, attr_name: str) -> str:
    """
    Extract attribute value from BioSample XML attribute blocks like:
      <Attribute attribute_name="strain" ...>value</Attribute>
    Tries multiple case variations to handle inconsistencies.
    """
    for variant in [attr_name, attr_name.replace("_", " "), attr_name.replace(" ", "_")]:
        pattern = rf'<Attribute\s[^>]*attribute_name="{re.escape(variant)}"[^>]*>(.*?)</Attribute>'
        m = re.search(pattern, xml, re.DOTALL | re.IGNORECASE)
        if m:
            return m.group(1).strip()
    return ""


def _id_value(xml: str, db_label: str) -> str:
    """
    Extract value from <Id db_label="...">value</Id> blocks.
    Used for extracting 'Sample name' and other ID fields.
    """
    pattern = rf'<Id\s+[^>]*db_label="{re.escape(db_label)}"[^>]*>(.*?)</Id>'
    m = re.search(pattern, xml, re.IGNORECASE)
    return m.group(1).strip() if m else ""


def _attr_first(*values: str) -> str:
    """Return the first non-empty value from a list of candidates."""
    return next((v for v in values if v), "")


def parse_biosample_xml(xml_block: str, accession: str) -> dict:
    """
    Parse a single <BioSample>…</BioSample> XML block into a flat dict.
    Handles messy/missing fields gracefully.
    """
    # ── Accession & sample name ───────────────────────────────────────────────
    acc_match = re.search(r'accession="(SAMN\d+|SAME\d+|SAMD\d+)"', xml_block)
    biosample_acc = acc_match.group(1) if acc_match else accession

    # Sample name: try multiple sources
    sample_name = _id_value(xml_block, "Sample name")
    if not sample_name:
        sample_name = _xml_text(xml_block, "SampleName")
    if not sample_name:
        m = re.search(r'<Id\s+db="SRA"\s*>(.*?)</Id>', xml_block)
        sample_name = m.group(1).strip() if m else ""

    # ── Organism ──────────────────────────────────────────────────────────────
    organism = _xml_text(xml_block, "OrganismName")
    if not organism:
        organism = _xml_text(xml_block, "ScientificName")

    # ── Strain / isolate ─────────────────────────────────────────────────────
    strain = _attr_first(
        _attr_value(xml_block, "strain"),
        _attr_value(xml_block, "isolate"),
        _attr_value(xml_block, "isolation_source"),
    )

    # ── Collection date ───────────────────────────────────────────────────────
    collection_date = _attr_first(
        _attr_value(xml_block, "collection_date"),
        _attr_value(xml_block, "collection date"),
    )

    # ── Geographic location ───────────────────────────────────────────────────
    geo_loc = _attr_first(
        _attr_value(xml_block, "geo_loc_name"),
        _attr_value(xml_block, "geographic location"),
        _attr_value(xml_block, "country"),
    )

    # ── Host / source ─────────────────────────────────────────────────────────
    host = _attr_first(
        _attr_value(xml_block, "host"),
        _attr_value(xml_block, "isolation_source"),
        _attr_value(xml_block, "source"),
    )

    # ── Source material identifiers ───────────────────────────────────────────
    source_material_id = _attr_first(
        _attr_value(xml_block, "source_material_id"),
        _attr_value(xml_block, "source material identifiers"),
    )
    
  # ── Comments / Treatment info ─────────────────────────────────────────────────
    comment = _xml_text(
            _xml_text(_xml_text(xml_block, "Description"), "Comment"),
            "Paragraph"
        )
    
      


    # ── BioProject ───────────────────────────────────────────────────────────
    bp_match = re.search(r'<Id\s+db="BioProject"\s*>(.*?)</Id>', xml_block)
    bioproject = bp_match.group(1).strip() if bp_match else ""

    return {
        "biosample_accession": biosample_acc,
        "sample_name":         sample_name,
        "organism":            organism,
        "strain_isolate":      strain,
        "collection_date":     collection_date,
        "geo_loc_name":        geo_loc,
        "host_source":         host,
        "source_material_id":  source_material_id,
        "bioproject_accession": bioproject,
        "comment":comment,
        # SRA fields filled in later
        "sra_experiment":      "",
        "sra_run":             "",
        "platform":            "",
        "library_strategy":    "",
        "library_layout":      "",
        "paired_end":          "",
        "fastq_ftp":           "",
        "fastq_aspera":        "",
        "fastq_bytes":         "",
        "fastq_md5":           "",
    }


def fetch_biosample_records(ncbi: NCBIClient, accessions: list[str]) -> list[dict]:
    """
    Fetch BioSample XML for a list of accessions and return parsed dicts.
    Processes in batches of 20 to stay well within URL-length limits.
    """
    records: list[dict] = []
    BATCH = 20

    for i in range(0, len(accessions), BATCH):
        batch = accessions[i : i + BATCH]
        log.info("Fetching BioSample metadata for batch %d–%d of %d",
                 i + 1, min(i + BATCH, len(accessions)), len(accessions))

        # Step 1: esearch to get internal NCBI UIDs
        uid_map: dict[str, str] = {}     # accession → uid
        for acc in batch:
            uids = ncbi.esearch("biosample", f"{acc}[accn]")
            if uids:
                uid_map[acc] = uids[0]
            else:
                log.warning("  No UID found for %s — skipping", acc)

        if not uid_map:
            continue

        # Step 2: efetch XML in one call
        xml_text = ncbi.efetch("biosample", list(uid_map.values()),
                               rettype="xml", retmode="xml")
  
        # Step 3: split into individual BioSample blocks and parse
        blocks = re.findall(r"<BioSample[\s>].*?</BioSample>", xml_text,
                            re.DOTALL | re.IGNORECASE)

        parsed_accs: set[str] = set()
        for block in blocks:
            rec = parse_biosample_xml(block, "")
            if rec["biosample_accession"]:
                records.append(rec)
                parsed_accs.add(rec["biosample_accession"])

        # Emit placeholder rows for accessions that returned no XML block
        for acc in batch:
            if acc not in parsed_accs:
                log.warning("  No XML block parsed for %s", acc)
                records.append({
                    "biosample_accession": acc,
                    "sample_name": "", "organism": "", "strain_isolate": "",
                    "collection_date": "", "geo_loc_name": "", "host_source": "",
                    "bioproject_accession": "", "sra_experiment": "",
                    "sra_run": "", "platform": "", "library_strategy": "",
                    "library_layout": "", "paired_end": "",
                    "fastq_ftp": "", "fastq_aspera": "",
                    "fastq_bytes": "", "fastq_md5": "",
                })
    return records


# ═══════════════════════════════════════════════════════════════════════════════
# 4.  SRA RUN METADATA (via NCBI esummary)
# ═══════════════════════════════════════════════════════════════════════════════

def fetch_sra_runs_for_biosample(ncbi: NCBIClient, biosample_acc: str) -> list[dict]:
    """
    Use NCBI elink (biosample → sra) + esummary to get run-level metadata.
    Returns a list of run dicts per BioSample.
    """
    # Search for the BioSample UID
    bs_uids = ncbi.esearch("biosample", f"{biosample_acc}[accn]")
    if not bs_uids:
        return []

    # Link biosample → sra
    sra_uids = ncbi.elink("biosample", "sra", bs_uids)
    if not sra_uids:
        return []

    # Fetch SRA summaries
    BATCH = 100
    run_records: list[dict] = []
    for i in range(0, len(sra_uids), BATCH):
        batch_uids = sra_uids[i : i + BATCH]
        summary = ncbi.esummary("sra", batch_uids)

        for uid, item in summary.items():
            if uid == "uids":
                continue

            # ExpXml contains nested XML with run/experiment details
            exp_xml   = item.get("expxml", "")
            runs_xml  = item.get("runs",   "")

            experiment_acc = _xml_text(exp_xml, "Experiment")
            if not experiment_acc:
                m = re.search(r'acc="(SRX\d+|ERX\d+|DRX\d+)"', exp_xml)
                experiment_acc = m.group(1) if m else ""

            platform = _xml_text(exp_xml, "Platform")
            if not platform:
                m = re.search(r'<Platform\s+instrument_model="([^"]+)"', exp_xml)
                platform = m.group(1) if m else item.get("platform", "")

            lib_strategy = _xml_text(exp_xml, "Strategy")
            lib_layout   = _xml_text(exp_xml, "Layout")
            if not lib_layout:
                m = re.search(r"<PAIRED|<SINGLE", exp_xml, re.IGNORECASE)
                lib_layout = "PAIRED" if (m and "PAIRED" in m.group()) else (
                             "SINGLE" if m else "")

            paired_end = "yes" if lib_layout.upper() == "PAIRED" else (
                         "no"  if lib_layout.upper() == "SINGLE" else "")

            # Parse individual run accessions from <Runs> block
            run_accs = re.findall(r'acc="(SRR\d+|ERR\d+|DRR\d+)"', runs_xml)
            if not run_accs:
                run_accs = [""]

            for run_acc in run_accs:
                run_records.append({
                    "sra_experiment":   experiment_acc,
                    "sra_run":          run_acc,
                    "platform":         platform,
                    "library_strategy": lib_strategy,
                    "library_layout":   lib_layout,
                    "paired_end":       paired_end,
                })

    return run_records


# ═══════════════════════════════════════════════════════════════════════════════
# 5.  ENA PORTAL API — FASTQ FILES (not SRA archives)
# ═══════════════════════════════════════════════════════════════════════════════

ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"

def fetch_ena_fastq_links(biosample_acc: str) -> list[dict]:
    """
    Query ENA Portal API for FASTQ files linked to a BioSample.
    Returns list of dicts with actual FASTQ file URLs (not SRA archives).
    """
    params = {
        "accession": biosample_acc,
        "result": "read_run",
        "fields": "run_accession,experiment_accession,sample_accession,"
                  "study_accession,scientific_name,library_strategy,"
                  "library_layout,instrument_platform,instrument_model,"
                  "fastq_ftp,fastq_aspera,fastq_bytes,fastq_md5",
        "format": "json",
        "limit": 500,
    }
    try:
        resp = http_get(ENA_PORTAL, params=params, timeout=30)
        runs = resp.json()
        if isinstance(runs, list):
            # Filter to only runs that have actual FASTQ files
            return [r for r in runs if r.get("fastq_ftp")]
        return []
    except Exception as exc:
        log.debug("ENA Portal FASTQ lookup failed for %s: %s", biosample_acc, exc)
        return []


def construct_sra_backup_urls(run_acc: str) -> dict:
    """
    Construct SRA archive URLs as FALLBACK (for runs without FASTQ on ENA).
    These are for conversion with fasterq-dump, not direct use.
    """
    # AWS S3 SRA archive
    aws_url = f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{run_acc}/{run_acc}"
    
    # NCBI SRA archive
    ncbi_url = f"https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos4/sra-pub-zq-1/{run_acc[:6]}/{run_acc[6:11]}/{run_acc}/{run_acc}.lite.1"
    
    return {
        "run_accession": run_acc,
        "aws_url": aws_url,
        "ncbi_url": ncbi_url,
        "sra_archive": True,  # Flag indicating this is an archive, not FASTQ
    }


# ═══════════════════════════════════════════════════════════════════════════════
# 6.  MERGE METADATA RECORDS
# ═══════════════════════════════════════════════════════════════════════════════

def merge_records(bs_record: dict,
                  sra_runs: list[dict],
                  ena_fastq_runs: list[dict]) -> list[dict]:
    """
    Produce one row per run by merging BioSample and SRA/ENA data.
    Priority: ENA FASTQ > SRA metadata > BioSample
    """
    merged: list[dict] = []

    # Build lookup from run_accession → ENA FASTQ data
    ena_by_run: dict[str, dict] = {
        r.get("run_accession", ""): r 
        for r in ena_fastq_runs
        if r.get("run_accession")
    }

    # Gather all unique run accessions
    all_run_accs: list[str] = list(set(
        r.get("sra_run", "") for r in sra_runs if r.get("sra_run")
    ))

    # Fallback: if no run at all, still write one metadata row
    if not all_run_accs:
        return [{**bs_record}]

    # SRA metadata lookup by run accession
    sra_by_run: dict[str, dict] = {r.get("sra_run", ""): r for r in sra_runs}

    for run_acc in all_run_accs:
        sra = sra_by_run.get(run_acc, {})
        ena = ena_by_run.get(run_acc, {})

        row = dict(bs_record)   # copy base BioSample fields
        row["sra_experiment"]   = ena.get("experiment_accession") or sra.get("sra_experiment") or ""
        row["sra_run"]          = run_acc
        row["platform"]         = ena.get("instrument_platform") or sra.get("platform") or ""
        row["library_strategy"] = ena.get("library_strategy") or sra.get("library_strategy") or ""
        row["library_layout"]   = ena.get("library_layout") or sra.get("library_layout") or ""
        layout = row["library_layout"].upper()
        row["paired_end"]       = "yes" if layout == "PAIRED" else ("no" if layout == "SINGLE" else "")

        # ENA FASTQ files (primary)
        row["fastq_ftp"]    = ena.get("fastq_ftp", "") or ""
        row["fastq_aspera"] = ena.get("fastq_aspera", "") or ""
        row["fastq_bytes"]  = ena.get("fastq_bytes", "") or ""
        row["fastq_md5"]    = ena.get("fastq_md5", "") or ""

        # SRA archive URLs (fallback if no FASTQ)
        row["aws_url"]  = ""
        row["ncbi_url"] = ""
        row["ena_ftp"]  = ""

        merged.append(row)

    return merged


# ═══════════════════════════════════════════════════════════════════════════════
# 7.  FILE DOWNLOAD
# ═══════════════════════════════════════════════════════════════════════════════

def _file_exists_and_complete(path: Path, expected_bytes: Optional[int]) -> bool:
    """Return True if the file already exists and has the expected size (and is not zero-byte)."""
    if not path.exists():
        return False
    
    size = path.stat().st_size
    
    # Never skip zero-byte files
    if size == 0:
        return False
    
    # If we have expected size, check it matches
    if expected_bytes and size != expected_bytes:
        return False
    
    return True


def download_file(url: str, dest: Path, expected_bytes: Optional[int] = None) -> bool:
    """
    Download *url* to *dest* using wget or curl (supports HTTP/HTTPS/FTP/S3).
    Returns True on success (file exists and is non-zero).
    Prefers wget, falls back to curl.
    """
    if _file_exists_and_complete(dest, expected_bytes):
        log.info("    Skipping (already complete): %s", dest.name)
        return True

    dest.parent.mkdir(parents=True, exist_ok=True)

    # Delete zero-byte files before retrying
    if dest.exists() and dest.stat().st_size == 0:
        log.debug("    Removing zero-byte file: %s", dest.name)
        dest.unlink()

    # Try wget first (better resume support)
    cmd_wget = [
        "wget",
        "--continue",      # Resume partial downloads
        "--quiet",         # Minimal output
        "-O", str(dest),   # Output file
        url
    ]

    try:
        result = subprocess.run(cmd_wget, capture_output=True, text=True, timeout=300)
        if result.returncode == 0:
            # Validate file is not empty
            if dest.exists() and dest.stat().st_size > 0:
                log.debug("    wget succeeded, file size: %d bytes", dest.stat().st_size)
                return True
            elif dest.exists():
                log.warning("    wget created zero-byte file; deleting and retrying")
                dest.unlink()
                # Fall through to curl
            else:
                log.warning("    wget claimed success but file not created")
    except FileNotFoundError:
        pass  # wget not available, try curl below
    except subprocess.TimeoutExpired:
        log.warning("    wget timeout for %s", url)
        return False
    except Exception as exc:
        log.debug("    wget failed: %s", exc)

    # Fallback: try curl
    cmd_curl = [
        "curl",
        "-C", "-",        # Resume partial downloads
        "-L",             # Follow redirects
        "-f",             # Fail on HTTP errors
        "-s",             # Silent mode
        "-o", str(dest),  # Output file
        url
    ]

    try:
        result = subprocess.run(cmd_curl, capture_output=True, text=True, timeout=300)
        if result.returncode == 0:
            # Validate file is not empty
            if dest.exists() and dest.stat().st_size > 0:
                log.debug("    curl succeeded, file size: %d bytes", dest.stat().st_size)
                return True
            elif dest.exists():
                log.warning("    curl created zero-byte file; deleting")
                dest.unlink()
                log.error("    Download produced empty file from %s", url)
                return False
            else:
                log.error("    curl claimed success but file not created")
                return False
        else:
            log.debug("    curl failed with code %d for %s", result.returncode, url)
            if dest.exists() and dest.stat().st_size == 0:
                dest.unlink()
    except FileNotFoundError:
        log.error("    Neither wget nor curl found. Install one for downloads.")
        return False
    except subprocess.TimeoutExpired:
        log.warning("    curl timeout for %s", url)
        return False
    except Exception as exc:
        log.error("    Download failed for %s: %s", url, exc)
        return False

    return False


def download_via_sratoolkit(run_acc: str, dest_dir: Path) -> bool:
    """
    Fallback: use fasterq-dump (SRA Toolkit) to download a run.
    Requires 'fasterq-dump' to be on PATH.
    """
    dest_dir.mkdir(parents=True, exist_ok=True)
    existing = list(dest_dir.glob(f"{run_acc}*.fastq*"))
    if existing:
        log.info("    SRA files already present for %s — skipping", run_acc)
        return True

    cmd = [
        "fasterq-dump",
        "--outdir", str(dest_dir),
        "--split-files",
        "--threads", "4",
        "--progress",
        run_acc,
    ]
    log.info("    Running fasterq-dump for %s", run_acc)
    try:
        result = subprocess.run(cmd, capture_output=False, check=True)
        return result.returncode == 0
    except FileNotFoundError:
        log.error("    fasterq-dump not found. Install SRA Toolkit or use --metadata_only.")
        return False
    except subprocess.CalledProcessError as exc:
        log.error("    fasterq-dump failed for %s (exit %d)", run_acc, exc.returncode)
        return False


def download_run(row: dict, runs_dir: Path,
                 failed_log: list[dict]) -> int:
    """
    Download FASTQ files for a single run row.
    If no FASTQ available, optionally fall back to SRA archive + fasterq-dump.
    Returns the number of files successfully downloaded or already present.
    """
    bs_acc  = row.get("biosample_accession", "UNKNOWN")
    run_acc = row.get("sra_run", "")
    if not run_acc:
        log.warning("  Row for %s has no run accession — skipping download", bs_acc)
        return 0

    run_dir = runs_dir / bs_acc / run_acc
    run_dir.mkdir(parents=True, exist_ok=True)

    # ── Try ENA FASTQ links first ─────────────────────────────────────────────
    fastq_ftp = row.get("fastq_ftp", "") or ""
    
    if fastq_ftp:
        # ENA separates paired-end files with semicolons
        ftp_links = [u.strip() for u in fastq_ftp.split(";") if u.strip()]
        
        success_count = 0
        for ftp_url in ftp_links:
            # Extract filename from URL
            filename = Path(ftp_url.split("/")[-1]).name
            dest = run_dir / filename
            
            # Check if file already exists and is complete
            if dest.exists() and dest.stat().st_size > 0:
                log.info("  ✓ Already present: %s (%s)", filename, run_acc)
                success_count += 1
                continue
            
            log.info("  Downloading FASTQ: %s", filename)
            log.debug("    URL (FTP): %s", ftp_url)
            
            # Convert FTP to HTTPS for ENA
            https_url = ftp_url.replace("ftp://", "https://") \
                               .replace("ftp.sra.ebi.ac.uk", "era.ebi.ac.uk", 1)
            
            # Try HTTPS first
            log.debug("    Trying HTTPS: %s", https_url)
            ok = download_file(https_url, dest)
            
            # Try FTP if HTTPS failed
            if not ok:
                log.info("    HTTPS failed; trying FTP: %s", ftp_url)
                ok = download_file(ftp_url, dest)
            
            if ok:
                log.info("  ✓ Downloaded: %s (%d bytes)", filename, dest.stat().st_size)
                success_count += 1
            else:
                log.error("  ✗ Failed: %s", filename)
                failed_log.append({
                    "biosample": bs_acc, "run": run_acc,
                    "url": ftp_url, "reason": "FASTQ download failed (check URL validity)"
                })
        
        if success_count > 0:
            return success_count

    # ── Fallback: SRA Toolkit (fasterq-dump) ──────────────────────────────────
    log.info("  No FASTQ for %s/%s — trying fasterq-dump", bs_acc, run_acc)
    ok = download_via_sratoolkit(run_acc, run_dir)
    if not ok:
        failed_log.append({
            "biosample": bs_acc, "run": run_acc,
            "url": f"sra:{run_acc}", "reason": "fasterq-dump failed"
        })
        return 0
    # Count downloaded files
    return len(list(run_dir.glob("*.fastq*")))


# ═══════════════════════════════════════════════════════════════════════════════
# 8.  MAIN ORCHESTRATION
# ═══════════════════════════════════════════════════════════════════════════════

def load_accessions(csv_path: str, id_column: str) -> list[str]:
    """Read BioSample accessions from a CSV (or one-per-line plain text)."""
    df = pd.read_csv(csv_path, header=None if id_column == "0" else "infer",
                     dtype=str, skip_blank_lines=True)

    # Accept the first column if id_column not found
    if id_column not in df.columns:
        col = df.columns[0]
        log.warning("Column '%s' not found; using first column '%s'", id_column, col)
    else:
        col = id_column

    accs = df[col].dropna().str.strip().tolist()
    # Keep only recognisable BioSample accession patterns
    valid = [a for a in accs if re.match(r"^SAM[NED]\d+$", a, re.IGNORECASE)]
    invalid = set(accs) - set(valid)
    if invalid:
        log.warning("Ignoring %d unrecognised accession(s): %s",
                    len(invalid), ", ".join(sorted(invalid)[:5]))
    return valid


def write_metadata_csv(records: list[dict], out_path: Path):
    """Write all merged records to a CSV file."""
    columns = [
        "biosample_accession", "sample_name", "comment","organism", "strain_isolate",
        "collection_date", "geo_loc_name", "host_source", "source_material_id",
        "bioproject_accession", "sra_experiment", "sra_run", "platform",
        "library_strategy", "library_layout", "paired_end", "fastq_ftp",
        "fastq_aspera", "fastq_bytes", "fastq_md5",
    ]
    df = pd.DataFrame(records, columns=columns)
    df.to_csv(out_path, index=False)
    log.info("Metadata written → %s  (%d rows)", out_path, len(df))


def write_failed_log(failed: list[dict], out_path: Path):
    if not failed:
        return
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["biosample", "run", "url", "reason"])
        w.writeheader()
        w.writerows(failed)
    log.warning("Failed downloads logged → %s  (%d entries)", out_path, len(failed))


# ═══════════════════════════════════════════════════════════════════════════════
# 9.  CLI
# ═══════════════════════════════════════════════════════════════════════════════

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="ncbi_biosample_downloader",
        description="Download BioSample metadata and SRA/ENA FASTQ runs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Metadata only:
  python ncbi_biosample_downloader.py \\
      --input_csv biosamples.csv --id_column biosample_id \\
      --outdir ./output --email you@example.com --metadata_only

  # Full download with NCBI API key:
  python ncbi_biosample_downloader.py \\
      --input_csv biosamples.csv --id_column biosample_id \\
      --outdir ./output --email you@example.com --api_key XXXXXX
        """,
    )
    p.add_argument("--input_csv",  required=True,
                   help="CSV file with BioSample accessions.")
    p.add_argument("--id_column",  default="biosample_id",
                   help="Column name containing BioSample accessions (default: biosample_id).")
    p.add_argument("--outdir",     required=True,
                   help="Output directory (will be created if absent).")
    p.add_argument("--email",      required=True,
                   help="Your e-mail address (required by NCBI E-utilities).")
    p.add_argument("--api_key",    default="",
                   help="NCBI API key (optional; raises rate limit to 10 req/s).")
    p.add_argument("--metadata_only", action="store_true",
                   help="Fetch metadata only; skip all file downloads.")
    p.add_argument("--threads",    type=int, default=4,
                   help="Parallel download threads (default: 4).")
    p.add_argument("--verbose",    action="store_true",
                   help="Enable DEBUG logging.")
    return p


def main():
    args   = build_parser().parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # ── 0. Validate input ─────────────────────────────────────────────────────
    log.info("═" * 60)
    log.info("NCBI BioSample Downloader")
    log.info("Input CSV : %s", args.input_csv)
    log.info("Output dir: %s", outdir)
    log.info("═" * 60)

    accessions = load_accessions(args.input_csv, args.id_column)
    if not accessions:
        sys.exit("No valid BioSample accessions found in input CSV. Aborting.")
    log.info("Loaded %d BioSample accession(s)", len(accessions))

    # ── 1. Initialise NCBI client ─────────────────────────────────────────────
    ncbi = NCBIClient(email=args.email, api_key=args.api_key)

    # ── 2. Fetch BioSample metadata ───────────────────────────────────────────
    log.info("\n[Step 1/4] Fetching BioSample metadata from NCBI …")
    bs_records = fetch_biosample_records(ncbi, accessions)
    bs_by_acc  = {r["biosample_accession"]: r for r in bs_records}

    # ── 3. Fetch SRA run metadata per BioSample ───────────────────────────────
    log.info("\n[Step 2/4] Fetching SRA run metadata from NCBI …")
    sra_by_acc: dict[str, list[dict]] = {}
    for acc in tqdm(accessions, desc="SRA runs", unit="sample"):
        sra_by_acc[acc] = fetch_sra_runs_for_biosample(ncbi, acc)

    # ── 4. Fetch ENA FASTQ download links ──────────────────────────────────────
    log.info("\n[Step 3/4] Fetching ENA FASTQ links …")
    ena_fastq_by_acc: dict[str, list[dict]] = {}
    for acc in tqdm(accessions, desc="ENA FASTQ", unit="sample"):
        ena_fastq_by_acc[acc] = fetch_ena_fastq_links(acc)
        time.sleep(0.1)   # polite ENA rate limit

    # ── 5. Merge all records into a flat table ────────────────────────────────
    log.info("\n[Step 4/4] Merging records …")
    all_rows: list[dict] = []
    samples_with_runs: set[str] = set()

    for acc in accessions:
        bs_rec   = bs_by_acc.get(acc, {"biosample_accession": acc})
        sra_runs = sra_by_acc.get(acc, [])
        ena_fastq = ena_fastq_by_acc.get(acc, [])
        rows     = merge_records(bs_rec, sra_runs, ena_fastq)
        all_rows.extend(rows)
        if any(r.get("sra_run") for r in rows):
            samples_with_runs.add(acc)

    # ── 6. Write metadata CSV ─────────────────────────────────────────────────
    meta_path = outdir / "biosample_metadata.csv"
    write_metadata_csv(all_rows, meta_path)

    # ── 7. Download runs ───────────────────────────────────────────────────────
    runs_dir      = outdir / "runs"
    failed_log: list[dict] = []
    total_downloaded = 0
    total_failed     = 0

    if args.metadata_only:
        log.info("\n--metadata_only flag set; skipping downloads.")
    else:
        log.info("\n[Downloading] Fetching FASTQ files → %s", runs_dir)
        runs_dir.mkdir(parents=True, exist_ok=True)

        # Filter to rows that actually have a run accession
        download_rows = [r for r in all_rows if r.get("sra_run")]

        with ThreadPoolExecutor(max_workers=args.threads) as pool:
            futures = {
                pool.submit(download_run, row, runs_dir, failed_log): row
                for row in download_rows
            }
            for future in tqdm(as_completed(futures),
                               total=len(futures),
                               desc="Runs", unit="run"):
                n = future.result()
                total_downloaded += n
                if n == 0:
                    total_failed += 1

        # Write failed-download log
        failed_path = outdir / "failed_downloads.csv"
        write_failed_log(failed_log, failed_path)

    # ── 8. Summary ────────────────────────────────────────────────────────────
    log.info("\n%s", "═" * 60)
    log.info("SUMMARY")
    log.info("  BioSamples processed          : %d", len(accessions))
    log.info("  BioSamples with runs found    : %d", len(samples_with_runs))
    total_runs = len([r for r in all_rows if r.get("sra_run")])
    log.info("  Total run rows in metadata    : %d", total_runs)
    if not args.metadata_only:
        log.info("  Files successfully downloaded : %d", total_downloaded)
        log.info("  Runs with download failures   : %d", total_failed)
    log.info("  Metadata CSV                  : %s", meta_path)
    if not args.metadata_only:
        log.info("  Run files directory           : %s", runs_dir)
        if failed_log:
            log.info("  Failed downloads log          : %s",
                     outdir / "failed_downloads.csv")
    log.info("═" * 60)


if __name__ == "__main__":
    main()
