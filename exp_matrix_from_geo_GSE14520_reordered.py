#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Download GSE14520 expression matrix from GEO and write a reordered matrix where
tumor samples appear first, then non-tumor samples, in the format:

#NAME    GSM...   GSM...
#CLASS   tumor    tumor   non-tumor ...

Output:
  data/external/GSE14520_expression_matrix_reordered.csv   (tab-delimited)

Dependencies:
  pip install GEOparse pandas numpy
"""

from __future__ import annotations

import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import GEOparse


# -----------------------------
# Repo paths
# -----------------------------
def find_repo_root(start: Path) -> Path:
    """
    Find a repo root by searching upwards for a 'data' directory.
    Falls back to the current working directory.
    """
    cur = start.resolve()
    for _ in range(6):
        if (cur / "data").exists():
            return cur
        cur = cur.parent
    return Path.cwd().resolve()


REPO_ROOT = find_repo_root(Path(__file__).resolve())
EXTERNAL_DIR = REPO_ROOT / "data" / "external"
CACHE_DIR = EXTERNAL_DIR / "geo_cache"
EXTERNAL_DIR.mkdir(parents=True, exist_ok=True)
CACHE_DIR.mkdir(parents=True, exist_ok=True)

OUT_PATH = EXTERNAL_DIR / "GSE14520_expression_matrix_reordered.csv"


# -----------------------------
# Robust extraction utilities
# -----------------------------
def extract_expression_series(gsm: GEOparse.GSM) -> Optional[pd.Series]:
    """
    Extract probe-level expression series for a GSM.
    Tries common GEO column conventions and falls back to the first numeric column.
    Returns a Series indexed by probe ID.
    """
    if not hasattr(gsm, "table") or gsm.table is None or gsm.table.empty:
        return None

    tab = gsm.table

    # Probe/feature ID column
    if "ID_REF" in tab.columns:
        id_col = "ID_REF"
    elif "ID" in tab.columns:
        id_col = "ID"
    else:
        id_col = tab.columns[0]

    # Expression value column
    expr_col = None
    for c in tab.columns:
        if str(c).upper() == "VALUE":
            expr_col = c
            break

    if expr_col is None:
        candidates = [c for c in tab.columns if c != id_col]
        numeric = [c for c in candidates if pd.api.types.is_numeric_dtype(tab[c])]
        expr_col = numeric[0] if numeric else (candidates[0] if candidates else id_col)

    s = tab.set_index(id_col)[expr_col]
    s = pd.to_numeric(s, errors="coerce")
    s.name = gsm.name
    return s


def extract_genes_from_string(gene_string: str) -> List[str]:
    gene_string = str(gene_string).replace('"', "").replace("'", "").strip()
    separators = ["///", ";", ",", "//", "|", "\t"]

    for sep in separators:
        if sep in gene_string:
            genes = [g.strip() for g in gene_string.split(sep)]
            genes = [g for g in genes if g and not g.startswith("---") and g not in {"--", "-"}]
            clean = []
            for g in genes:
                if "(" in g:
                    g = g.split("(")[0].strip()
                if "[" in g:
                    g = g.split("[")[0].strip()
                if g and len(g) > 1:
                    clean.append(g)
            return clean[:1]

    if gene_string and not gene_string.startswith("---") and gene_string not in {"--", "-"} and len(gene_string) > 1:
        if "(" in gene_string:
            gene_string = gene_string.split("(")[0].strip()
        return [gene_string]

    return []


def guess_gene_symbol_col(columns: List[str]) -> Optional[str]:
    """
    Identify a gene symbol column in a GPL annotation table.
    Uses exact matches and fuzzy regex matches.
    """
    # Exact matches first
    exact = [
        "Gene Symbol", "Gene_Symbol", "GENE_SYMBOL", "Symbol", "SYMBOL",
        "gene_assignment", "Gene name", "GENE_NAME", "Gene", "GENE",
        "HGNC symbol", "HGNC Symbol", "hgnc_symbol", "gene_symbol"
    ]
    for c in exact:
        if c in columns:
            return c

    # Fuzzy matches
    for c in columns:
        cl = str(c).strip().lower()
        if re.fullmatch(r"(gene\s*symbol|symbol|hgnc\s*symbol)", cl):
            return c
        if "gene symbol" in cl or (cl == "symbol"):
            return c

    return None


def create_gene_mapping(gpl: GEOparse.GPL) -> Dict[str, str]:
    """
    Create mapping from probe IDs to HGNC gene symbols for one platform.
    """
    probe_to_gene: Dict[str, str] = {}
    if not hasattr(gpl, "table") or gpl.table is None or gpl.table.empty:
        return probe_to_gene

    gpl_table = gpl.table
    cols = list(gpl_table.columns)

    # Probe ID column
    probe_col = "ID" if "ID" in cols else cols[0]
    gene_col = guess_gene_symbol_col(cols)
    if gene_col is None:
        return probe_to_gene

    for _, row in gpl_table.iterrows():
        probe_id = row.get(probe_col, None)
        gene_info = row.get(gene_col, None)
        if pd.isna(probe_id) or pd.isna(gene_info):
            continue
        genes = extract_genes_from_string(str(gene_info))
        if genes:
            probe_to_gene[str(probe_id)] = genes[0]

    return probe_to_gene


def map_probes_to_genes(expr_probe: pd.DataFrame, probe_to_gene: Dict[str, str]) -> pd.DataFrame:
    """
    Convert probe-level matrix (probe x sample) to gene-level (gene x sample)
    by averaging probes mapping to the same gene symbol.
    """
    row_mapping = {probe: probe_to_gene.get(probe, None) for probe in expr_probe.index}
    row_mapping = {p: g for p, g in row_mapping.items() if g and g not in {"", "---", "--"}}

    gene_groups = defaultdict(list)
    for probe_id, gene_name in row_mapping.items():
        if probe_id in expr_probe.index:
            gene_groups[gene_name].append(probe_id)

    gene_data = []
    gene_names = []
    for gene_name, probes in gene_groups.items():
        if len(probes) == 1:
            gene_data.append(expr_probe.loc[probes[0]])
        else:
            gene_data.append(expr_probe.loc[probes].mean(axis=0))
        gene_names.append(gene_name)

    if not gene_data:
        # fallback: keep probes if mapping fails
        return expr_probe.copy()

    out = pd.DataFrame(gene_data, index=gene_names)
    out = out[~out.index.duplicated(keep="first")]
    return out


# -----------------------------
# Tumor / non-tumor labeling
# -----------------------------
_NON_TUMOR_PAT = re.compile(
    r"(non[\s\-]*tumou?r|adjacent|para[\s\-]*carcinoma|paracarcinoma|normal(\s*liver)?|control|non[\s\-]*cancer|non[\s\-]*malignan)",
    re.IGNORECASE,
)
_TUMOR_PAT = re.compile(
    r"(tumou?r|hcc|carcinoma|cancer|hepatocellular)",
    re.IGNORECASE,
)


def classify_sample_gse14520(gsm: GEOparse.GSM) -> str:
    """
    Return 'tumor', 'non-tumor', or 'unknown' based on GSM metadata text.
    Precedence: non-tumor patterns override tumor patterns.
    """
    md = getattr(gsm, "metadata", {}) or {}
    parts: List[str] = []

    for k in ["title", "source_name_ch1", "description", "characteristics_ch1"]:
        v = md.get(k, [])
        if isinstance(v, list):
            parts.extend([str(x) for x in v])
        elif v is not None:
            parts.append(str(v))

    text = " ; ".join(parts)
    if _NON_TUMOR_PAT.search(text):
        return "non-tumor"
    if _TUMOR_PAT.search(text):
        return "tumor"
    return "unknown"


# -----------------------------
# Main: download + build matrix
# -----------------------------
def download_gse14520_gene_matrix() -> Tuple[pd.DataFrame, Dict[str, str]]:
    """
    Download GSE14520 and return (gene_matrix, sample_class), where:
      - gene_matrix: genes x samples (HGNC symbols as index)
      - sample_class: dict GSM -> {'tumor','non-tumor','unknown'}
    Handles multiple platforms by mapping each platform separately then merging.
    """
    print("Downloading GSE14520 via GEOparse...")
    gse = GEOparse.get_GEO(geo="GSE14520", destdir=str(CACHE_DIR))

    # Prepare per-platform probe matrices
    platform_to_samples: Dict[str, Dict[str, pd.Series]] = defaultdict(dict)
    sample_class: Dict[str, str] = {}

    for gsm_id, gsm in gse.gsms.items():
        s = extract_expression_series(gsm)
        if s is None:
            continue

        # Platform id for this GSM
        md = gsm.metadata or {}
        platform_id = None
        if "platform_id" in md and isinstance(md["platform_id"], list) and md["platform_id"]:
            platform_id = md["platform_id"][0]
        if platform_id is None:
            # fallback: best-effort pick the first GPL (rarely needed)
            platform_id = list(gse.gpls.keys())[0]

        platform_to_samples[platform_id][gsm_id] = s
        sample_class[gsm_id] = classify_sample_gse14520(gsm)

    if not platform_to_samples:
        raise RuntimeError("No GSM expression tables were found for GSE14520 (gsm.table empty).")

    gene_mats: List[pd.DataFrame] = []
    for platform_id, sample_dict in platform_to_samples.items():
        print(f"Processing platform {platform_id} with {len(sample_dict)} samples...")
        expr_probe = pd.DataFrame(sample_dict)

        gpl = gse.gpls.get(platform_id, None)
        if gpl is None:
            print(f"  WARNING: platform {platform_id} not present in GSE object; keeping probes.")
            gene_mats.append(expr_probe)
            continue

        probe_to_gene = create_gene_mapping(gpl)
        if not probe_to_gene:
            print(f"  WARNING: could not build probe->gene mapping for {platform_id}; keeping probes.")
            gene_mats.append(expr_probe)
            continue

        expr_gene = map_probes_to_genes(expr_probe, probe_to_gene)
        gene_mats.append(expr_gene)

    # Merge platforms at gene level (outer join)
    gene_mat = pd.concat(gene_mats, axis=1, join="outer")
    gene_mat.index = gene_mat.index.astype(str)

    # Ensure stable gene order
    gene_mat = gene_mat.sort_index()

    # Ensure stable sample order BEFORE reordering (alphabetical within class)
    gene_mat = gene_mat.reindex(sorted(gene_mat.columns), axis=1)

    return gene_mat, sample_class


def write_reordered_matrix_with_class_header(
    gene_mat: pd.DataFrame, sample_class: Dict[str, str], out_path: Path
) -> None:
    """
    Reorder columns: tumor first, then non-tumor, then unknown.
    Write tab-delimited file with two header rows (#NAME / #CLASS) and no extra column header.
    """
    cols = list(gene_mat.columns)

    tumor_cols = sorted([c for c in cols if sample_class.get(c, "unknown") == "tumor"])
    nont_cols = sorted([c for c in cols if sample_class.get(c, "unknown") == "non-tumor"])
    unk_cols = sorted([c for c in cols if sample_class.get(c, "unknown") == "unknown"])

    ordered = tumor_cols + nont_cols + unk_cols
    if not ordered:
        raise RuntimeError("No columns to write after classification.")

    print(f"Column order: tumor={len(tumor_cols)}, non-tumor={len(nont_cols)}, unknown={len(unk_cols)}")

    gene_mat = gene_mat[ordered]

    # Write file
    with out_path.open("w", encoding="utf-8", newline="\n") as f:
        f.write("#NAME\t" + "\t".join(ordered) + "\n")
        f.write("#CLASS\t" + "\t".join([sample_class.get(c, "unknown") for c in ordered]) + "\n")

        # gene rows
        for gene, row in gene_mat.iterrows():
            vals = row.values
            # Convert to strings, blanks for missing
            sval = ["" if (v is None or (isinstance(v, float) and np.isnan(v))) else str(v) for v in vals]
            f.write(str(gene) + "\t" + "\t".join(sval) + "\n")

    print(f"Wrote: {out_path}")


def main() -> None:
    gene_mat, sample_class = download_gse14520_gene_matrix()
    write_reordered_matrix_with_class_header(gene_mat, sample_class, OUT_PATH)


if __name__ == "__main__":
    main()
