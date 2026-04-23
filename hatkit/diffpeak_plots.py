#!/usr/bin/env python3
"""
volcano_refactored.py — HATkit module
======================================
Generates volcano plots, MA plots, and replicate-correlation regression from
DESeq2 results (diffpeak or diffgene output).  Optionally highlights genomic
regions of interest and checks for obfuscated HBG1/HBG2 annotations.

Requires: bedtools (must be on PATH)
"""

from __future__ import annotations

import argparse
import gzip
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from adjustText import adjust_text
from scipy.stats import linregress, spearmanr
from tqdm import tqdm

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Configuration dataclass — replaces all module-level globals
# ---------------------------------------------------------------------------

@dataclass
class VolcanoConfig:
    """
    Carries every user-facing option through the analysis pipeline.
    Eliminates the need for global variables.
    """
    file: str
    regions: list[str] = field(default_factory=list)
    reverse_fc: bool = False
    use_pvalue: bool = False
    lfc_thresh: float = 2.0
    fdr_thresh: float = 0.01
    title: Optional[str] = None
    label_volcano: bool = False
    num_labels: int = 10
    label_specific: list[str] = field(default_factory=list)
    label_promoter_only: bool = False
    promoter_window: int = 2000
    override_hbg: bool = False
    dge: bool = False
    genome: str = "hg19"
    debug: bool = False

    # Derived at runtime
    outname: str = ""
    region_names: list[str] = field(default_factory=list)

    def __post_init__(self):
        self.outname = os.path.splitext(os.path.basename(self.file))[0]
        self.region_names = [
            os.path.splitext(os.path.basename(r))[0].split(".")[0]
            for r in self.regions
        ]
        if self.promoter_window not in (1000, 2000, 5000):
            raise ValueError("--promoterWindow must be 1000, 2000, or 5000")

    @property
    def y_col(self) -> str:
        return "pVal" if self.use_pvalue else "adj_pVal"

    @property
    def y_label(self) -> str:
        return "-log(p-val)" if self.use_pvalue else "-logFDR"

    # Highlight colour palette (extend as needed)
    HIGHLIGHT_COLORS = ["orange", "green", "brown"]


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def read_tsv(path: str, **kwargs) -> pd.DataFrame:
    """Read a tab-separated file, returning a DataFrame."""
    return pd.read_csv(path, sep="\t", **kwargs)


def write_tsv(df: pd.DataFrame, path: str, **kwargs) -> None:
    """Write a DataFrame to a tab-separated file."""
    df.to_csv(path, sep="\t", index=False, **kwargs)


def safe_copy(src: str, dst: str) -> None:
    """Copy *src* to *dst* using Python (no subprocess)."""
    shutil.copy2(src, dst)


def convert_scientific_notation(input_path: str, output_path: str) -> None:
    """
    Rewrite *input_path* to *output_path*, expanding any scientific-notation
    numeric fields to standard decimal notation (10 d.p.).

    Parameters
    ----------
    input_path  : source file path
    output_path : destination file path
    """
    pattern = re.compile(r"^[+-]?\d+(\.\d*)?[eE][+-]?\d+$")
    with open(input_path) as fin, open(output_path, "w") as fout:
        for line in fin:
            fields = line.rstrip("\n").split("\t")
            for i, f in enumerate(fields):
                if pattern.match(f):
                    try:
                        fields[i] = f"{float(f):.10f}"
                    except ValueError:
                        pass
            fout.write("\t".join(fields) + "\n")


# ---------------------------------------------------------------------------
# Genome / annotation helpers
# ---------------------------------------------------------------------------

UCSC_BASE = "https://hgdownload.soe.ucsc.edu/goldenPath"

def _stream_download(url: str, dest: Path) -> None:
    """Stream-download *url* to *dest*, showing a progress bar."""
    with requests.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        with open(dest, "wb") as fh, tqdm(
            total=total, unit="B", unit_scale=True, desc=dest.name
        ) as bar:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    fh.write(chunk)
                    bar.update(len(chunk))


def ensure_gtf(genome: str = "hg19", cache_dir: Optional[str] = None) -> Path:
    """
    Return the path to the refGene GTF (gzipped), downloading it if absent.

    Parameters
    ----------
    genome    : UCSC genome assembly (e.g. "hg19", "hg38")
    cache_dir : where to cache the file; defaults to ~/.hatkit/cache/
    """
    cache = Path(cache_dir or Path.home() / ".hatkit" / "cache")
    cache.mkdir(parents=True, exist_ok=True)
    gtf = cache / f"{genome}.refGene.gtf.gz"
    if not gtf.exists():
        url = f"{UCSC_BASE}/{genome}/bigZips/genes/{genome}.refGene.gtf.gz"
        logger.info("Downloading GTF: %s", url)
        _stream_download(url, gtf)
    return gtf


def build_transcript_map(gtf_path: Path) -> dict[str, tuple[str, str]]:
    """
    Parse a refGene GTF and return a transcript-ID → (gene_name, strand) map.

    Parameters
    ----------
    gtf_path : path to gzipped GTF file

    Returns
    -------
    dict mapping transcript_id → (gene_name, strand)
    """
    tx_map: dict[str, tuple[str, str]] = {}
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9 or parts[2] != "transcript":
                continue
            strand = parts[6]
            attrs = {}
            for attr in parts[8].split(";"):
                kv = attr.strip().split(" ", 1)
                if len(kv) == 2:
                    attrs[kv[0]] = kv[1].strip('"')
            tx_id = attrs.get("transcript_id")
            if tx_id:
                gene = attrs.get("gene_name") or attrs.get("gene_id") or "Unknown_gene"
                tx_map[tx_id] = (gene, strand)
    return tx_map


def fetch_gene_info(
    header_row: str,
    tx_map: dict[str, tuple[str, str]],
) -> tuple[str, str]:
    """
    Extract (gene_name, strand) for a FASTA header like '>NM_006625_up_1000_chr…'.

    Parameters
    ----------
    header_row : FASTA header string
    tx_map     : transcript map from :func:`build_transcript_map`

    Returns
    -------
    (gene_name, strand) or ("Unknown_gene", ".")
    """
    m = re.match(r"^>([A-Z]{2}_[0-9]+)", header_row)
    if not m:
        return "Unknown_gene", "."
    return tx_map.get(m.group(1), ("Unknown_gene", "."))


def ensure_promoter_bed(
    genome: str = "hg19",
    window: int = 2000,
    cache_dir: Optional[str] = None,
) -> Path:
    """
    Return the path to a promoter BED file, building it from UCSC upstream
    FASTA data if it does not already exist.

    Parameters
    ----------
    genome    : assembly (e.g. "hg19")
    window    : upstream window in bp — must be 1000, 2000, or 5000
    cache_dir : cache directory (default ~/.hatkit/cache/)

    Returns
    -------
    Path to upstream{window}.bed
    """
    cache = Path(cache_dir or Path.home() / ".hatkit" / "cache")
    cache.mkdir(parents=True, exist_ok=True)
    bed_path = cache / f"{genome}_upstream{window}.bed"

    if bed_path.exists() and bed_path.stat().st_size > 0:
        return bed_path

    fa_path = cache / f"{genome}_upstream{window}.fa"
    if not fa_path.exists() or fa_path.stat().st_size == 0:
        url = f"{UCSC_BASE}/{genome}/bigZips/upstream{window}.fa.gz"
        gz_path = cache / f"{genome}_upstream{window}.fa.gz"
        logger.info("Downloading upstream FASTA: %s", url)
        _stream_download(url, gz_path)
        subprocess.run(["gunzip", "-f", str(gz_path)], check=True)

    gtf = ensure_gtf(genome, cache_dir)
    tx_map = build_transcript_map(gtf)

    logger.info("Writing promoter BED (window=%d). This runs once per genome/window.", window)
    total = sum(1 for _ in fa_path.open())
    with fa_path.open() as fa, bed_path.open("w") as bed:
        for line in tqdm(fa, total=total, desc="Building promoter BED"):
            if not line.startswith(">"):
                continue
            m = re.search(r"(\w+):(\d+)-(\d+)", line)
            if m:
                chrom, start, end = m.group(1), m.group(2), m.group(3)
                gene, strand = fetch_gene_info(line.strip(), tx_map)
                bed.write(f"{chrom}\t{start}\t{end}\t{gene}\t{strand}\n")

    return bed_path


def intersects_promoter(
    gene_row: pd.Series,
    genome: str = "hg19",
    window: int = 2000,
    cache_dir: Optional[str] = None,
) -> bool:
    """
    Return True if *gene_row* overlaps the promoter BED for *genome*/*window*.

    Parameters
    ----------
    gene_row  : row from the DESeq2 DataFrame (needs 'chrom','start','end')
    genome    : assembly
    window    : promoter window in bp
    cache_dir : cache directory
    """
    bed_path = ensure_promoter_bed(genome, window, cache_dir)
    chrom = gene_row["chrom"]
    g_start = int(gene_row["start"])
    g_end = int(gene_row["end"])
    with bed_path.open() as fh:
        for line in fh:
            b_chrom, b_start, b_end, *_ = line.strip().split("\t")
            if b_chrom == chrom and g_start <= int(b_end) and g_end >= int(b_start):
                return True
    return False


# ---------------------------------------------------------------------------
# HBG override (bedtools-based, no temp shell scripts)
# ---------------------------------------------------------------------------

HBG_PROMOTERS = [
    ("chr11", 5241401, 5243401, "HBG1"),
    ("chr11", 5245100, 5247100, "HBG2"),
]


def check_hbg_label(
    input_df: pd.DataFrame,
) -> tuple[Optional[pd.DataFrame], bool, bool]:
    """
    Check whether HBG1/HBG2 promoter regions are present but obfuscated in
    *input_df* (e.g. merged into an adjacent peak during HOMER annotation).

    Runs ``bedtools intersect`` entirely in memory via named temporary files;
    no shell scripts are written to disk.

    Parameters
    ----------
    input_df : formatted DESeq2 DataFrame (must have chrom/start/end columns)

    Returns
    -------
    (intersect_df_or_None, hbg1_found, hbg2_found)
    """
    _check_tool("bedtools")

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        a_bed = tmp / "query.bed"
        b_bed = tmp / "hbg_promoters.bed"
        out_bed = tmp / "intersect.bed"

        # Write query BED (chrom/start/end + remaining columns)
        query_cols = ["chrom", "start", "end"] + [
            c for c in input_df.columns if c not in ("chrom", "start", "end")
        ]
        # Filter out unplaced contigs
        filtered = input_df[
            ~input_df["chrom"].astype(str).str.contains(r"chrUn_|random|-", regex=True)
        ]
        filtered[query_cols].to_csv(a_bed, sep="\t", header=False, index=False)

        # Write HBG promoter BED
        with b_bed.open("w") as fh:
            for chrom, start, end, name in HBG_PROMOTERS:
                fh.write(f"{chrom}\t{start}\t{end}\t{name}\n")

        result = subprocess.run(
            ["bedtools", "intersect", "-a", str(a_bed), "-b", str(b_bed), "-wao"],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            logger.warning("bedtools intersect failed: %s", result.stderr)
            return None, False, False

        lines = [l for l in result.stdout.strip().splitlines() if l]
        if not lines:
            return None, False, False

        n_query_cols = len(query_cols)
        extra = ["hbg_chr", "hbg_start", "hbg_end", "hbg_name", "overlap"]
        all_cols = query_cols + extra

        rows = [line.split("\t") for line in lines]
        df = pd.DataFrame(rows, columns=all_cols)
        df["overlap"] = pd.to_numeric(df["overlap"], errors="coerce").fillna(0)
        df = df[df["overlap"] > 0]

        if df.empty:
            return None, False, False

        df = df.loc[df.groupby("hbg_name")["overlap"].idxmax()].set_index("hbg_name")
        return df, "HBG1" in df.index, "HBG2" in df.index


# ---------------------------------------------------------------------------
# DESeq2 input parsing
# ---------------------------------------------------------------------------

# Standard column sets
_PEAK_COLS = [
    "title", "chrom", "start", "end", "strand", "peakScore", "regionSize",
    "annotation", "detailed_annot", "dist2TSS", "promID", "entrezID",
    "nearestUnigene", "nearestRefseq", "nearestEnsemble", "geneName",
    "geneAlias", "geneDesc", "geneType",
    "backgroundrep1TagDensity", "backgroundrep2TagDensity",
    "treatrep1TagDensity", "treatrep2TagDensity",
    "logFC", "pVal", "adj_pVal",
]
_PEAK_KEEP = [
    "chrom", "start", "end", "strand", "geneName",
    "backgroundrep1TagDensity", "backgroundrep2TagDensity",
    "treatrep1TagDensity", "treatrep2TagDensity",
    "logFC", "pVal", "adj_pVal",
]
_DGE_COLS = [
    "targetID", "pVal", "adj_pVal", "geneName", "num_agg", "X",
    "treatrep1TagDensity", "treatrep2TagDensity",
    "backgroundrep1TagDensity", "backgroundrep2TagDensity",
    "treat_avgTD", "background_avgTD", "logFC",
]
_DGE_KEEP = [
    "pVal", "adj_pVal", "geneName",
    "treatrep1TagDensity", "treatrep2TagDensity",
    "backgroundrep1TagDensity", "backgroundrep2TagDensity",
    "logFC", "treat_avgTD", "background_avgTD",
]
_NUMERIC_COLS = [
    "pVal", "adj_pVal", "logFC",
    "treatrep1TagDensity", "treatrep2TagDensity",
    "backgroundrep1TagDensity", "backgroundrep2TagDensity",
    "treat_avgTD", "background_avgTD",
]


def _log10_transform_column(series: pd.Series) -> pd.Series:
    """Apply log10 to values >= 1; set values < 1 to 0."""
    return series.apply(lambda v: np.log10(v) if v >= 1 else 0.0)


def _extract_sample_label(col_name: str) -> str:
    """
    Extract a human-readable sample label from a HOMER tag-density column name.

    HOMER column names follow the pattern:
        <condition>_<replicate>_tag/<path>  or  <condition>_<replicate>_tag/...

    We strip anything from the first '/' onward (path noise), then join the
    first two underscore-delimited fields so that e.g.
        '14_treated_2_tag/...'  →  '14_treated'
        'untreated_rep1_tag'    →  'untreated_rep1'

    Parameters
    ----------
    col_name : raw column name from the TSV header

    Returns
    -------
    Label string composed of the first two '_'-delimited tokens.
    """
    clean = col_name.split("/")[0]       # drop path suffix
    parts = clean.split("_")
    # Return first two tokens if available, otherwise whatever exists
    return "_".join(parts[:2]) if len(parts) >= 2 else parts[0]


def _validate_replicate_labels(
    label_a: str, label_b: str, group: str
) -> None:
    """
    Warn if the two replicate columns for a sample group do not share the
    same derived label, which would indicate an unexpected column ordering
    or mixed samples.

    Parameters
    ----------
    label_a : label derived from the first replicate column
    label_b : label derived from the second replicate column
    group   : human-readable group name ('treat' or 'background') for the
              warning message
    """
    if label_a != label_b:
        logger.warning(
            "Sample label mismatch in %s columns: '%s' vs '%s'. "
            "Check column ordering or use --reverseFC if groups are swapped.",
            group, label_a, label_b,
        )


def parse_deseq2(path: str, dge: bool = False) -> tuple[pd.DataFrame, dict]:
    """
    Read a DESeq2 result file and return a standardised DataFrame plus a
    metadata dict describing the sample columns found.

    For peak files the last seven data columns before logFC/pVal/adj_pVal are:
        [-7] background rep1 tag density
        [-6] background rep2 tag density
        [-5] treat rep1 tag density
        [-4] treat rep2 tag density
        [-3] logFC
        [-2] pVal
        [-1] adj_pVal

    Sample labels are derived by splitting each column name on '_' and joining
    the first two tokens (e.g. '14_treated_2_tag/...' → '14_treated').

    Parameters
    ----------
    path : path to the DESeq2 result TSV (or CSV for DGE)
    dge  : True for differential gene expression output

    Returns
    -------
    (df, meta)
    meta keys:
        't1_label'  – derived label for the first column group (cols 19–20)
        't2_label'  – derived label for the second column group (cols 21–22)
        't1_cols'   – raw column names for group 1
        't2_cols'   – raw column names for group 2
    """
    if dge:
        raw = pd.read_csv(path, sep=",")
        # _DGE_COLS maps: col[6]=treatrep1, col[7]=treatrep2,
        #                  col[8]=backgroundrep1, col[9]=backgroundrep2
        # t1 = background group (cols 8,9); t2 = treat/test group (cols 6,7)
        t1_cols = [raw.columns[8], raw.columns[9]]   # background rep1, rep2
        t2_cols = [raw.columns[6], raw.columns[7]]   # treat rep1, rep2
        t1_label = _extract_sample_label(t1_cols[0])
        t2_label = _extract_sample_label(t2_cols[0])
        meta = {
            "t1_label": t1_label,
            "t2_label": t2_label,
            "t1_cols": t1_cols,
            "t2_cols": t2_cols,
        }
        raw.columns = _DGE_COLS
        df = raw[_DGE_KEEP].copy()
        for col in _NUMERIC_COLS:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")

    else:
        raw = pd.read_csv(path, sep="\t")
        # Peak output: cols at fixed positions 19–22 are the four tag-density
        # columns.  Positions may vary if HOMER adds extra annotation fields,
        # so we also support reading from the tail of the header.
        # The last 7 columns (before logFC/pVal/adjpVal) are:
        #   …, bg_rep1, bg_rep2, treat_rep1, treat_rep2, logFC, pVal, adj_pVal
        # which in zero-based absolute indices are cols 19–22 for standard
        # HOMER output with 26 total columns.
        cols = list(raw.columns)
        # Use tail-relative positions so the code is robust to extra columns
        bg_rep1_col   = cols[-7]   # background rep 1
        bg_rep2_col   = cols[-6]   # background rep 2
        treat_rep1_col = cols[-5]  # treat rep 1
        treat_rep2_col = cols[-4]  # treat rep 2

        t1_label_a = _extract_sample_label(bg_rep1_col)
        t1_label_b = _extract_sample_label(bg_rep2_col)
        t2_label_a = _extract_sample_label(treat_rep1_col)
        t2_label_b = _extract_sample_label(treat_rep2_col)

        _validate_replicate_labels(t1_label_a, t1_label_b, "background")
        _validate_replicate_labels(t2_label_a, t2_label_b, "treat")

        meta = {
            "t1_label": t1_label_a,   # background group label
            "t2_label": t2_label_a,   # treat group label
            "t1_cols": [bg_rep1_col, bg_rep2_col],
            "t2_cols": [treat_rep1_col, treat_rep2_col],
        }
        raw.columns = _PEAK_COLS
        df = raw[_PEAK_KEEP].copy()
        df["background_avgTD"] = (
            df["backgroundrep1TagDensity"] + df["backgroundrep2TagDensity"]
        ) / 2
        df["treat_avgTD"] = (
            df["treatrep1TagDensity"] + df["treatrep2TagDensity"]
        ) / 2

    # Fill object-column NAs
    for col in df.select_dtypes(include="object").columns:
        df[col] = df[col].fillna("N/A")

    return df, meta


def resolve_background_treat(
    df: pd.DataFrame, meta: dict, reverse_fc: bool = False
) -> tuple[str, str, list[str], list[str]]:
    """
    Determine which sample group is 'background' and which is 'treat' based
    on the sign of the first-row LFC and the relative average tag densities,
    then return clean human-readable labels derived directly from the column
    names (not single-token prefixes).

    The logic mirrors the original script:
        - positive LFC + avg1 > avg2  →  t2 is treat  (higher = treat)
        - positive LFC + avg1 ≤ avg2  →  t1 is treat
        - negative LFC + avg1 > avg2  →  t1 is treat  (higher = background)
        - negative LFC + avg1 ≤ avg2  →  t2 is treat

    Parameters
    ----------
    df         : DataFrame from :func:`parse_deseq2`
    meta       : metadata dict from :func:`parse_deseq2`
    reverse_fc : if True, swap background and treat assignment

    Returns
    -------
    (back_label, treat_label, back_cols, treat_cols)
    """
    lfc_positive = pd.to_numeric(df["logFC"], errors="coerce").iloc[0] > 0
    avg1 = pd.to_numeric(df["background_avgTD"], errors="coerce").iloc[0]
    avg2 = pd.to_numeric(df["treat_avgTD"], errors="coerce").iloc[0]

    t1_label = meta["t1_label"]   # background group (cols -7/-6)
    t2_label = meta["t2_label"]   # treat group (cols -5/-4)
    t1_cols  = meta["t1_cols"]
    t2_cols  = meta["t2_cols"]

    if lfc_positive:
        if avg1 > avg2:
            back, treat, bcols, tcols = t2_label, t1_label, t2_cols, t1_cols
        else:
            back, treat, bcols, tcols = t1_label, t2_label, t1_cols, t2_cols
    else:
        if avg1 > avg2:
            back, treat, bcols, tcols = t1_label, t2_label, t2_cols, t1_cols
        else:
            back, treat, bcols, tcols = t2_label, t1_label, t1_cols, t2_cols

    if reverse_fc:
        back, treat, bcols, tcols = treat, back, tcols, bcols

    # Strip path suffixes from raw column names (kept for replicate regression)
    bcols = [c.split("/")[0] for c in bcols]
    tcols = [c.split("/")[0] for c in tcols]

    return back, treat, bcols, tcols


# ---------------------------------------------------------------------------
# Bedtools intersection
# ---------------------------------------------------------------------------

def _check_tool(name: str) -> None:
    """Raise RuntimeError if *name* is not on PATH."""
    if shutil.which(name) is None:
        raise RuntimeError(
            f"'{name}' not found on PATH. Please install it and try again."
        )


def bed_intersect(
    query_df: pd.DataFrame,
    region_path: str,
    query_cols: Optional[list[str]] = None,
) -> pd.DataFrame:
    """
    Run ``bedtools intersect -wo`` between *query_df* and a BED region file.

    Parameters
    ----------
    query_df    : DataFrame with at least chrom/start/end columns
    region_path : path to a BED file of regions of interest
    query_cols  : column names for *query_df*; defaults to df.columns

    Returns
    -------
    Intersect DataFrame with overlap column appended, or empty DataFrame.
    """
    _check_tool("bedtools")
    query_cols = list(query_df.columns) if query_cols is None else query_cols

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        a = tmp / "query.bed"
        b = tmp / "region.bed"
        out = tmp / "intersect.bed"

        # Write only chrom/start/end for the region BED
        region_df = pd.read_csv(region_path, sep="\t", header=None)
        region_df.iloc[:, :3].to_csv(b, sep="\t", header=False, index=False)

        # Write query without header
        query_df[query_cols[:len(query_df.columns)]].to_csv(
            a, sep="\t", header=False, index=False
        )

        result = subprocess.run(
            ["bedtools", "intersect", "-wo", "-a", str(a), "-b", str(b)],
            capture_output=True, text=True,
        )
        if result.returncode != 0 or not result.stdout.strip():
            return pd.DataFrame()

        rows = [line.split("\t") for line in result.stdout.strip().splitlines()]
        extra = ["chrom2", "start2", "end2", "overlap"]
        all_cols = query_cols + extra
        hdf = pd.DataFrame(rows, columns=all_cols)
        # Cast numeric columns — bedtools stdout is always plain text
        for col in _NUMERIC_COLS:
            if col in hdf.columns:
                hdf[col] = pd.to_numeric(hdf[col], errors="coerce")
        return hdf


def collect_region_highlights(
    df: pd.DataFrame,
    region_files: list[str],
) -> list[pd.DataFrame]:
    """
    Intersect *df* with each BED file in *region_files* and return a list of
    intersection DataFrames (one per region file).

    Parameters
    ----------
    df           : standardised DESeq2 DataFrame
    region_files : list of BED file paths

    Returns
    -------
    list of DataFrames (may be empty DataFrames if no intersection)
    """
    highlights = []
    for rfile in region_files:
        hit = bed_intersect(df, rfile)
        highlights.append(hit)
    return highlights


# ---------------------------------------------------------------------------
# Labelling helpers (shared between plotVolcano and plotHighlight)
# ---------------------------------------------------------------------------

def _add_top_labels(
    ax: plt.Axes,
    df: pd.DataFrame,
    y_col: str,
    n: int,
) -> None:
    """
    Add auto-adjusted text labels for the top-*n* up- and down-regulated
    entries in *df*.

    Parameters
    ----------
    ax    : matplotlib Axes
    df    : DataFrame indexed by geneName with logFC / y_col columns
    y_col : column to use for the y-axis
    n     : how many labels per direction
    """
    for direction in ("up", "down"):
        subset = df[df["logFC"] > 0] if direction == "up" else df[df["logFC"] < 0]
        subset = subset.sort_values(y_col, ascending=True).head(n)
        texts = []
        for gene, row in subset.iterrows():
            texts.append(
                ax.text(row["logFC"], -np.log10(row[y_col]), gene, fontsize=5)
            )
        if texts:
            adjust_text(
                texts,
                force_text=(0.5, 0.2),
                arrowprops=dict(arrowstyle="-", color="black", lw=0.5),
            )


def _add_specific_labels(
    ax: plt.Axes,
    df: pd.DataFrame,
    genes: list[str],
    y_col: str,
    cfg: VolcanoConfig,
    hbg_override: Optional[pd.DataFrame] = None,
    hbg1_found: bool = False,
    hbg2_found: bool = False,
) -> None:
    """
    Plot and label specific genes by name, handling HBG override logic and
    optional promoter-region filtering.

    Parameters
    ----------
    ax           : matplotlib Axes
    df           : DataFrame indexed by geneName
    genes        : list of gene names to highlight
    y_col        : column for y-axis values
    cfg          : VolcanoConfig
    hbg_override : HBG intersection DataFrame (from check_hbg_label)
    hbg1/2_found : whether HBG1/2 were found in the override check
    """
    texts: list = []

    for gene in genes:
        try:
            info = df.loc[gene]
        except KeyError:
            logger.info("%s not found in DataFrame — skipping", gene)
            continue

        # --- HBG override path ---
        if gene in ("HBG1", "HBG2") and cfg.override_hbg:
            hbg_found = hbg1_found if gene == "HBG1" else hbg2_found
            gene_in_df = df.index.isin([gene]).any()
            if not gene_in_df:
                if hbg_override is None or not hbg_found:
                    continue
                hbg_row = hbg_override.loc[gene]
                xs = hbg_row["logFC"] if isinstance(hbg_row["logFC"], pd.Series) else [hbg_row["logFC"]]
                ys = hbg_row[y_col] if isinstance(hbg_row[y_col], pd.Series) else [hbg_row[y_col]]
                for xi, yi in zip(xs, ys):
                    texts.append(ax.text(float(xi), -np.log10(float(yi)), gene, fontsize=5, fontweight="bold"))
                continue

        # --- Promoter-only path ---
        if cfg.label_promoter_only and not cfg.dge:
            rows_to_label = []
            if isinstance(info, pd.DataFrame):
                for _, row in info.iterrows():
                    if intersects_promoter(row, cfg.genome, cfg.promoter_window):
                        rows_to_label.append(row)
            else:
                if intersects_promoter(info, cfg.genome, cfg.promoter_window):
                    rows_to_label.append(info)
            if not rows_to_label:
                logger.debug("%s does not intersect a promoter — skipping", gene)
                continue
            for row in rows_to_label:
                ax.scatter(row["logFC"], -np.log10(row[y_col]), s=4, color="black")
                texts.append(ax.text(row["logFC"], -np.log10(row[y_col]), gene, fontsize=5, fontweight="bold"))
            continue

        # --- Standard path ---
        xs = info["logFC"] if isinstance(info, pd.Series) else info["logFC"]
        ys = info[y_col] if isinstance(info, pd.Series) else info[y_col]

        if isinstance(xs, pd.Series):
            for xi, yi in zip(xs, ys):
                texts.append(ax.text(float(xi), -np.log10(float(yi)), gene, fontsize=5, fontweight="bold"))
        else:
            ax.scatter(float(xs), -np.log10(float(ys)), s=4, color="black")
            texts.append(ax.text(float(xs), -np.log10(float(ys)), gene, fontsize=5, fontweight="bold"))

    if texts:
        adjust_text(texts, arrowprops=dict(arrowstyle="-", color="black", lw=1))


# ---------------------------------------------------------------------------
# Threshold helpers
# ---------------------------------------------------------------------------

def classify_regulation(
    df: pd.DataFrame, y_col: str, lfc_thresh: float, fdr_thresh: float
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Split *df* into up- and down-regulated subsets.

    Parameters
    ----------
    df         : DataFrame with logFC and y_col columns
    y_col      : 'pVal' or 'adj_pVal'
    lfc_thresh : log-fold-change threshold
    fdr_thresh : significance threshold

    Returns
    -------
    (up_df, down_df)
    """
    up = df[(df["logFC"] >= lfc_thresh) & (df[y_col] <= fdr_thresh)]
    down = df[(df["logFC"] <= -lfc_thresh) & (df[y_col] <= fdr_thresh)]
    return up, down


def _scatter_sig(ax: plt.Axes, df: pd.DataFrame, y_col: str, **kw) -> None:
    """Scatter-plot points using -log10(y_col) as y-axis."""
    y = pd.to_numeric(df[y_col], errors="coerce")
    ax.scatter(df["logFC"], y.apply(lambda v: -np.log10(v)), **kw)


def _draw_thresholds(ax: plt.Axes, lfc: float, fdr: float) -> None:
    """Draw dashed LFC and FDR threshold lines."""
    ax.axvline(-lfc, color="black", linestyle="--", alpha=0.6)
    ax.axvline(lfc, color="black", linestyle="--", alpha=0.6)
    ax.axhline(-np.log10(fdr), color="black", linestyle="--", alpha=0.6)


# ---------------------------------------------------------------------------
# Volcano plot
# ---------------------------------------------------------------------------

def plot_volcano(
    df: pd.DataFrame,
    cfg: VolcanoConfig,
    highlights: Optional[list[pd.DataFrame]] = None,
    output_dir: str = ".",
) -> tuple[float, float, float, float]:
    """
    Produce the main volcano plot (base + full with highlights/labels).

    Parameters
    ----------
    df         : standardised DESeq2 DataFrame
    cfg        : VolcanoConfig
    highlights : list of region-highlight DataFrames
    output_dir : where to save output files

    Returns
    -------
    (xl, xr, yl, yu) — axis limits, passed to plot_highlight()
    """
    y_col = cfg.y_col
    out = Path(output_dir)

    if cfg.dge:
        df = df[df[y_col].notna()] if cfg.use_pvalue else df[df["logFC"].notna()]

    fig, ax = plt.subplots()

    _scatter_sig(ax, df, y_col, s=1, label="Not significant", color="grey")
    up, down = classify_regulation(df, y_col, cfg.lfc_thresh, cfg.fdr_thresh)
    _scatter_sig(ax, down, y_col, s=3, label="Down-regulated", color="blue")
    _scatter_sig(ax, up,   y_col, s=3, label="Up-regulated",   color="red")

    write_tsv(down, str(out / "down_regulated.tsv"))
    write_tsv(up,   str(out / "up_regulated.tsv"))

    ax.set_xlabel("logFC")
    ax.set_ylabel(cfg.y_label)
    _draw_thresholds(ax, cfg.lfc_thresh, cfg.fdr_thresh)
    ax.legend(fontsize=8)
    plot_title = cfg.title or cfg.outname
    ax.set_title(plot_title)

    # Save base (no highlights yet)
    fig.savefig(out / "base_volcanoPlot.pdf", format="pdf")
    fig.savefig(out / "base_volcanoPlot.png", format="png")

    # Overlay highlights
    colors = cfg.HIGHLIGHT_COLORS
    if highlights:
        for i, (hdf, hname) in enumerate(zip(highlights, cfg.region_names)):
            if hdf is None or hdf.empty:
                continue
            color = colors[i % len(colors)]
            _scatter_sig(ax, hdf, y_col, s=3, label=hname, color=color)

    # Labels
    if cfg.label_volcano or cfg.label_specific:
        hbg_override = hbg1_found = hbg2_found = None
        if cfg.override_hbg:
            hbg_override, hbg1_found, hbg2_found = check_hbg_label(df)

        idx_df = df.set_index("geneName") if "geneName" in df.columns else df

        if cfg.label_volcano:
            _add_top_labels(ax, idx_df, y_col, cfg.num_labels)

        if cfg.label_specific:
            _add_specific_labels(
                ax, idx_df, cfg.label_specific, y_col, cfg,
                hbg_override, hbg1_found, hbg2_found,
            )

    ax.legend(fontsize=8)
    fig.savefig(out / "full_volcanoPlot.pdf", format="pdf")
    fig.savefig(out / "full_volcanoPlot.png", format="png")

    xl, xr = ax.get_xlim()
    yl, yu = ax.get_ylim()
    plt.close(fig)
    return xl, xr, yl, yu


def plot_highlight(
    highlight_df: pd.DataFrame,
    name: str,
    color: str,
    xl: float, xr: float, yl: float, yu: float,
    cfg: VolcanoConfig,
    full_df: Optional[pd.DataFrame] = None,
    output_dir: str = ".",
) -> None:
    """
    Produce a standalone volcano plot for a single highlight region set,
    using the same axis limits as the full plot.

    Parameters
    ----------
    highlight_df : subset DataFrame from bed_intersect
    name         : region file stem (used in output filenames)
    color        : scatter colour for highlighted points
    xl/xr/yl/yu  : axis limits from plot_volcano()
    cfg          : VolcanoConfig
    full_df      : full DESeq2 DataFrame (needed for specific-label lookup)
    output_dir   : where to save output files
    """
    y_col = cfg.y_col
    out = Path(output_dir)

    up, down = classify_regulation(highlight_df, y_col, cfg.lfc_thresh, cfg.fdr_thresh)
    write_tsv(up,   str(out / f"up_highlighted_{name}.tsv"))
    write_tsv(down, str(out / f"down_highlighted_{name}.tsv"))
    up.to_csv(out / f"up_highlighted_{name}.bed",   sep="\t", index=False, header=False)
    down.to_csv(out / f"down_highlighted_{name}.bed", sep="\t", index=False, header=False)

    fig, ax = plt.subplots()
    _scatter_sig(ax, highlight_df, y_col, s=3, label=name, color=color)
    ax.set_title(cfg.title or f"{cfg.outname}_highlight")
    ax.set_xlabel("logFC")
    ax.set_ylabel(cfg.y_label)
    ax.set_xlim(xl, xr)
    ax.set_ylim(yl, yu)
    _draw_thresholds(ax, cfg.lfc_thresh, cfg.fdr_thresh)

    if cfg.label_volcano or cfg.label_specific:
        hbg_override = hbg1_found = hbg2_found = None
        if cfg.override_hbg and full_df is not None:
            hbg_override, hbg1_found, hbg2_found = check_hbg_label(full_df)

        src_df = full_df if full_df is not None else highlight_df
        idx_df = src_df.set_index("geneName") if "geneName" in src_df.columns else src_df

        if cfg.label_volcano:
            _add_top_labels(ax, idx_df, y_col, cfg.num_labels)
        if cfg.label_specific:
            _add_specific_labels(
                ax, idx_df, cfg.label_specific, y_col, cfg,
                hbg_override, hbg1_found, hbg2_found,
            )

    handles, labels_ = ax.get_legend_handles_labels()
    ax.legend(dict(zip(labels_, handles)).values(),
              dict(zip(labels_, handles)).keys(), fontsize=8)

    fig.savefig(out / f"highlighted_{name}.pdf", format="pdf")
    fig.savefig(out / f"highlighted_{name}.png", format="png")
    plt.close(fig)


# ---------------------------------------------------------------------------
# MA plot
# ---------------------------------------------------------------------------

def _render_ma(
    df: pd.DataFrame,
    padj_col: str,
    alpha: float,
    title: str,
    figsize: tuple,
    point_size: int,
    symmetric_ylim: bool,
    out_path_stem: Path,
) -> None:
    """
    Draw and save a single MA plot from a pre-computed DataFrame that already
    has 'A' (mean average) and 'M' (fold change) columns.

    Parameters
    ----------
    df             : DataFrame with 'A', 'M', and padj_col columns
    padj_col       : column name for adjusted p-value
    alpha          : significance cutoff
    title          : plot title
    figsize        : figure size tuple
    point_size     : scatter dot size
    symmetric_ylim : force symmetric y-axis around 0
    out_path_stem  : Path without extension; .pdf and .png are appended
    """
    plot_df = df.copy()
    plot_df["category"] = "ns"
    plot_df.loc[(plot_df[padj_col] < alpha) & (plot_df["M"] > 0), "category"] = "up"
    plot_df.loc[(plot_df[padj_col] < alpha) & (plot_df["M"] < 0), "category"] = "down"

    palette = {"ns": "grey", "up": "red", "down": "blue"}
    fig, ax = plt.subplots(figsize=figsize)
    for cat, color in palette.items():
        sub = plot_df[plot_df["category"] == cat]
        ax.scatter(sub["A"], sub["M"], s=point_size, color=color,
                   alpha=0.7 if cat == "ns" else 1.0)
    ax.axhline(0, linestyle="--", linewidth=1, color="black")
    if symmetric_ylim:
        lim = np.nanmax(np.abs(plot_df["M"]))
        ax.set_ylim(-lim, lim)
    ax.set_xlabel("A (mean log2 tag density)")
    ax.set_ylabel("M (log2 fold change)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(str(out_path_stem) + ".pdf", format="pdf")
    fig.savefig(str(out_path_stem) + ".png", format="png")
    plt.close(fig)


def plot_ma(
    df: pd.DataFrame,
    background_col: str = "background_avgTD",
    treat_col: str = "treat_avgTD",
    logfc_col: str = "logFC",
    padj_col: str = "adj_pVal",
    alpha: float = 0.05,
    figsize: tuple = (6, 6),
    point_size: int = 10,
    title: str = "MA Plot",
    symmetric_ylim: bool = True,
    iqr_fence: float = 3.0,
    output_dir: str = ".",
) -> None:
    """
    Generate and save two MA plots: a full plot and an IQR-filtered plot that
    excludes extreme A-axis outliers so the bulk of the data is visible.

    Parameters
    ----------
    df             : DESeq2 DataFrame
    background_col : column for mean background tag density
    treat_col      : column for mean treat tag density
    logfc_col      : column for log2 fold change
    padj_col       : column for adjusted p-value
    alpha          : significance cutoff
    figsize        : figure size tuple
    point_size     : scatter dot size
    title          : plot title
    symmetric_ylim : force symmetric y-axis around 0
    iqr_fence      : IQR multiplier for outlier exclusion on the A axis
                     (default 3.0 — points beyond Q3 + 3*IQR are removed)
    output_dir     : where to save output files
    """
    out = Path(output_dir)
    ma_df = df.copy()
    ma_df["A"] = pd.to_numeric(
        (ma_df[background_col] + ma_df[treat_col]) / 2, errors="coerce"
    )
    ma_df["M"] = pd.to_numeric(ma_df[logfc_col], errors="coerce")
    ma_df = ma_df.dropna(subset=["A", "M"])

    # --- Full MA plot ---
    _render_ma(ma_df, padj_col, alpha, title, figsize, point_size,
               symmetric_ylim, out / "MA_plot")

    # --- IQR-filtered MA plot ---
    q1 = ma_df["A"].quantile(0.25)
    q3 = ma_df["A"].quantile(0.75)
    iqr = q3 - q1
    upper_fence = q3 + iqr_fence * iqr
    filtered_df = ma_df[ma_df["A"] <= upper_fence]
    n_removed = len(ma_df) - len(filtered_df)
    logger.debug(
        "MA plot IQR filter: fence=%.2f, removed %d outlier(s) from A axis.",
        upper_fence, n_removed,
    )
    _render_ma(filtered_df, padj_col, alpha, f"{title} (IQR filtered)",
               figsize, point_size, symmetric_ylim, out / "MA_plot_filtered")


# ---------------------------------------------------------------------------
# Replicate correlation regression
# ---------------------------------------------------------------------------

def replicate_regression(
    df: pd.DataFrame,
    group_label: str,
    rep_cols: tuple[str, str],
    cfg: VolcanoConfig,
    output_dir: str = ".",
) -> None:
    """
    Linear + Spearman regression on the two replicates of one sample group
    (either test or background) and save a dot-dot plot, outlier file, and
    regression stats — all prefixed with the group name so outputs from both
    groups can coexist.

    Parameters
    ----------
    df          : standardised DESeq2 DataFrame (columns already renamed to
                  original input labels)
    group_label : human-readable name for the group (e.g. '14_treated'),
                  used in output filenames and plot titles
    rep_cols    : (rep1_col, rep2_col) — the actual DataFrame column names for
                  the two replicates of this group, after renaming
    cfg         : VolcanoConfig
    output_dir  : where to save output files
    """
    out = Path(output_dir)

    # Sanitise the label for use in filenames (spaces → underscores)
    file_stem = group_label.replace(" ", "_")

    x_col, y_col = rep_cols

    x_raw = df[x_col]
    y_raw = df[y_col]
    # For DGE: apply log10 here (not stored in df to keep TSV output clean)
    # For peak: columns are log2 tag densities; convert back to linear
    if cfg.dge:
        x_vals = x_raw.apply(lambda v: np.log10(v) if v >= 1 else 0.0)
        y_vals = y_raw.apply(lambda v: np.log10(v) if v >= 1 else 0.0)
    else:
        x_vals = 2 ** x_raw
        y_vals = 2 ** y_raw

    valid = ~pd.isnull(x_vals) & ~pd.isnull(y_vals) & np.isfinite(x_vals) & np.isfinite(y_vals)
    xv, yv = x_vals[valid], y_vals[valid]

    if len(xv) < 2:
        logger.warning("Fewer than 2 valid points for %s; skipping regression.", group_label)
        slope = intercept = r_value = p_value = std_err = spearman_corr = spearman_p = np.nan
    else:
        slope, intercept, r_value, p_value, std_err = linregress(xv, yv)
        spearman_corr, spearman_p = spearmanr(xv, yv)

    # Outlier detection (deviation from identity line)
    residuals = y_vals - x_vals
    valid_resid = ~pd.isnull(residuals) & np.isfinite(residuals)
    resid_std = residuals[valid_resid].std()
    mask_out = valid_resid & (residuals.abs() > 3 * resid_std)
    outliers = df[mask_out].copy()
    if not outliers.empty:
        outliers["x_val"] = x_vals[mask_out].values
        outliers["y_val"] = y_vals[mask_out].values

    outlier_file = out / f"{file_stem}_replicate_outliers.txt"
    outliers.to_csv(outlier_file, sep="\t", index=False)

    # Save regression stats
    stats_file = out / f"{file_stem}_linregress_stats.txt"
    with stats_file.open("w") as fh:
        fh.write(
            f"GROUP:\t{group_label}\n"
            f"SLOPE:\t{slope}\nINTERCEPT:\t{intercept}\n"
            f"PEARSON R-VALUE:\t{r_value}\nPEARSON P-VALUE:\t{p_value}\n"
            f"PEARSON STD ERR:\t{std_err}\n"
            f"SPEARMAN CORR:\t{spearman_corr}\nSPEARMAN P-VALUE:\t{spearman_p}\n"
        )

    # --- Plot helper ---
    xlabel = "Replicate 1 log10(RNA counts)" if cfg.dge else "Replicate 1 Tag Density"
    ylabel = "Replicate 2 log10(RNA counts)" if cfg.dge else "Replicate 2 Tag Density"
    if cfg.title:
        base_title = f"{cfg.title.split(':')[0]}:{group_label.upper()}_Replicate_Correlation"
    else:
        base_title = f"{group_label} Replicate Correlation"

    def _render_replicate_plot(
        xv_plot: pd.Series,
        yv_plot: pd.Series,
        title: str,
        out_stem: Path,
    ) -> None:
        sl, ic, rv, pv, se = linregress(xv_plot, yv_plot)
        sc, _ = spearmanr(xv_plot, yv_plot)
        fig, ax = plt.subplots()
        ax.scatter(xv_plot, yv_plot, alpha=0.5 if not cfg.dge else 0.6)
        ax.axline((0, ic), slope=sl, color="black")
        max_val = max(xv_plot.max(), yv_plot.max())
        padding = 0.05 * max_val
        ax.set_xlim(0, max_val + padding)
        ax.set_ylim(0, max_val + padding)
        ax.text(0.01, 0.99, f"R²={rv**2:.3f}", ha="left", va="top", transform=ax.transAxes)
        ax.text(0.01, 0.95, f"Spearman={sc:.3f}", ha="left", va="top", transform=ax.transAxes)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        fig.savefig(str(out_stem) + ".pdf", format="pdf")
        fig.savefig(str(out_stem) + ".png", format="png")
        plt.close(fig)

    # Full plot
    _render_replicate_plot(
        xv, yv, base_title,
        out / f"{file_stem}_replicate_correlation_plot",
    )

    # Conditional IQR-filtered plot: produce if BOTH conditions are met on
    # either axis:
    #   1. max / fence > 2.0        (axis stretched beyond double the fence
    #      — the bulk of data is being compressed)
    #   2. n_extreme / n_total < 0.01  (fewer than 1% of points lie beyond
    #      10× the fence — distinguishes a few truly extreme points from
    #      genuine heavy-tailed spread where many points are very high)
    iqr_fence = 3.0
    extreme_mult = 10.0   # multiplier for "truly extreme" point detection
    n_total = len(xv)
    produce_filtered = False
    for axis_vals in (xv, yv):
        q3    = axis_vals.quantile(0.75)
        iqr   = axis_vals.quantile(0.75) - axis_vals.quantile(0.25)
        fence = q3 + iqr_fence * iqr
        ratio = axis_vals.max() / fence if fence > 0 else 0
        n_extreme = (axis_vals > extreme_mult * fence).sum()
        frac_extreme = n_extreme / n_total
        logger.debug(
            "Replicate plot IQR check (%s): max/fence=%.2f, "
            "n_extreme(>%.0fx)=%d, frac_extreme=%.4f",
            group_label, ratio, extreme_mult, n_extreme, frac_extreme,
        )
        if ratio > 2.0 and frac_extreme < 0.01:
            produce_filtered = True
            break

    if produce_filtered:
        xq3  = xv.quantile(0.75); xiqr = xv.quantile(0.75) - xv.quantile(0.25)
        yq3  = yv.quantile(0.75); yiqr = yv.quantile(0.75) - yv.quantile(0.25)
        x_cut = xq3 + extreme_mult * xiqr
        y_cut = yq3 + extreme_mult * yiqr
        mask = (xv <= x_cut) & (yv <= y_cut)
        xv_filt, yv_filt = xv[mask], yv[mask]
        logger.debug(
            "Replicate plot IQR filter (%s): removed %d outlier(s) "
            "(x_cut=%.1f, y_cut=%.1f).",
            group_label, len(xv) - mask.sum(), x_cut, y_cut,
        )
        if len(xv_filt) >= 2:
            _render_replicate_plot(
                xv_filt, yv_filt,
                f"{base_title} (IQR filtered)",
                out / f"{file_stem}_replicate_correlation_plot_filtered",
            )


# ---------------------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------------------

def cleanup(
    output_dir: str,
    has_regions: bool = False,
    keep_temp: bool = False,
) -> None:
    """
    Move all outputs into *output_dir*. Optionally delete temporary files.

    Parameters
    ----------
    output_dir  : destination directory name
    has_regions : True if region-highlight files were produced
    keep_temp   : if True, leave temp* files in place (debug mode)
    """
    out = Path(output_dir)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    globs = ["*.pdf", "*.png", "down_regulated.tsv", "up_regulated.tsv",
             "*_replicate_outliers.txt", "*_linregress_stats.txt"]
    if has_regions:
        globs += ["up_highlighted*", "down_highlighted*", "region_intersect*.tsv"]

    for pattern in globs:
        for f in Path(".").glob(pattern):
            shutil.move(str(f), out / f.name)

    if not keep_temp:
        for f in Path(".").glob("temp*"):
            f.unlink(missing_ok=True)
        for f in ["hbg_promoters.bed", "hbg_override_check.sh"]:
            Path(f).unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    """Return the argument parser for the volcano CLI."""
    p = argparse.ArgumentParser(
        description="HATkit volcano — volcano/MA/replicate plots from DESeq2 output"
    )
    p.add_argument("-f", "--file", required=True, help="DESeq2 result file")
    p.add_argument(
        "-r", "--regions", nargs="*", default=None,
        help="Optional BED region file(s) of interest (up to 3)",
    )
    p.add_argument("--reverseFC", action="store_true", help="Swap background/treat assignment")
    p.add_argument("-p", "--usePvalue", action="store_true", help="Use raw p-value instead of FDR")
    p.add_argument("-lt", "--lfcThresh", type=float, default=2.0, help="LFC threshold (default 2)")
    p.add_argument("-pt", "--pThresh", type=float, default=0.01, help="p/FDR threshold (default 0.01)")
    p.add_argument("-t", "--title", default=None, help="Plot title")
    p.add_argument("--labelVolcano", action="store_true", help="Auto-label top DE genes")
    p.add_argument("--numberLabels", type=int, default=10, help="Labels per direction (default 10)")
    p.add_argument("--labelSpecific", nargs="*", default=None, help="Specific genes to label")
    p.add_argument("--labelPromoterOnly", action="store_true", help="Only label genes at promoters")
    p.add_argument("--promoterWindow", type=int, default=2000, choices=[1000, 2000, 5000],
                   help="Promoter window in bp (default 2000)")
    p.add_argument("--overrideHBG", action="store_true", help="Enable HBG1/HBG2 obfuscation check")
    p.add_argument("--dge", action="store_true", help="Input is differential gene expression output")
    p.add_argument("--genome", default="hg19", help="UCSC genome assembly (default hg19)")
    p.add_argument("--debug", action="store_true", help="Keep temp files; enable verbose logging")
    return p


def main(argv: Optional[list[str]] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format="%(levelname)s: %(message)s",
    )
    # matplotlib's font resolver is extremely verbose at DEBUG; keep it quiet
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("PIL").setLevel(logging.WARNING)

    cfg = VolcanoConfig(
        file=args.file,
        regions=args.regions or [],
        reverse_fc=args.reverseFC,
        use_pvalue=args.usePvalue,
        lfc_thresh=args.lfcThresh,
        fdr_thresh=args.pThresh,
        title=args.title,
        label_volcano=args.labelVolcano,
        num_labels=args.numberLabels,
        label_specific=args.labelSpecific or [],
        label_promoter_only=args.labelPromoterOnly,
        promoter_window=args.promoterWindow,
        override_hbg=args.overrideHBG,
        dge=args.dge,
        genome=args.genome,
        debug=args.debug,
    )

    logger.debug("Config: %s", cfg)

    # --- Parse input ---
    safe_copy(cfg.file, "temp.tsv")
    df, meta = parse_deseq2("temp.tsv", dge=cfg.dge)
    back, treat, back_cols, treat_cols = resolve_background_treat(df, meta, cfg.reverse_fc)
    if cfg.reverse_fc:
        df["logFC"] = -df["logFC"]

    logger.info("Background: %s | Test: %s", back, treat)

    # --- Rename internal columns to original input labels ---
    # Applied before any write_tsv calls (plot_volcano, plot_highlight).
    # plot_ma receives explicit col names so it is unaffected by the rename.
    _col_rename: dict[str, str] = {
        "backgroundrep1TagDensity": back_cols[0],
        "backgroundrep2TagDensity": back_cols[1],
        "treatrep1TagDensity":      treat_cols[0],
        "treatrep2TagDensity":      treat_cols[1],
        "treat_avgTD":              f"{_extract_sample_label(treat)}_avgTD",
        "background_avgTD":         f"{_extract_sample_label(back)}_avgTD",
    }
    df.rename(columns=_col_rename, inplace=True)

    # --- Region intersections ---
    highlights: list[pd.DataFrame] = []
    if cfg.regions and not cfg.dge:
        highlights = collect_region_highlights(df, cfg.regions)

    # --- Volcano ---
    xl, xr, yl, yu = plot_volcano(df, cfg, highlights=highlights or None)

    # --- Per-region highlight plots ---
    if highlights:
        colors = cfg.HIGHLIGHT_COLORS
        for i, (hdf, hname) in enumerate(zip(highlights, cfg.region_names)):
            if hdf is not None and not hdf.empty:
                plot_highlight(hdf, hname, colors[i % len(colors)],
                               xl, xr, yl, yu, cfg, full_df=df)

    # --- MA plot ---
    # Pass renamed avgTD column names explicitly so plot_ma works after rename.
    treat_avg_col = f"{_extract_sample_label(treat)}_avgTD"
    back_avg_col  = f"{_extract_sample_label(back)}_avgTD"
    ma_kwargs: dict = {
        "background_col": back_avg_col if not cfg.reverse_fc else treat_avg_col,
        "treat_col":      treat_avg_col if not cfg.reverse_fc else back_avg_col,
    }
    if cfg.use_pvalue:
        ma_kwargs["padj_col"] = "pVal"
    plot_ma(df, title=cfg.title or "MA Plot", **ma_kwargs)

    # --- Replicate regression (both groups) ---
    # rep_cols use the renamed column names now in df.
    is_test = Path(cfg.file).name.startswith("test")
    if not is_test:
        replicate_regression(df, treat, rep_cols=(treat_cols[0], treat_cols[1]), cfg=cfg)
        replicate_regression(df, back,  rep_cols=(back_cols[0],  back_cols[1]),  cfg=cfg)
    else:
        logger.warning("Test dataframe detected — skipping replicate regression.")

    # --- Cleanup ---
    cleanup(cfg.outname, has_regions=bool(highlights), keep_temp=cfg.debug)

    print(f"\nOutputs:    {cfg.outname}/")
    print(f"Background: {back}")
    print(f"Test:       {treat}")
    print("Re-run with --reverseFC to swap groups.\n")


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    main()
