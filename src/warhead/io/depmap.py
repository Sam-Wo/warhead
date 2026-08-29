"""DepMap 26Q1 loader (model metadata + doubling time, expression, CRISPR).

Doubling-time metadata is the input to G2b. DepMap distributes the release as a
set of CSVs from the portal download page (depmap.org/portal/download); drop them
into ``data/raw/depmap/``. Expected files and their tidy output schemas are
documented per function. Reuses the intent of the existing
``depmap_adc_payload_scoring.py`` loader.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from ..config import DATA_RAW
from .base import RawMissingError

DEPMAP_RAW = DATA_RAW / "depmap"

# Candidate column names DepMap has used for doubling time across releases.
_DOUBLING_CANDIDATES = [
    "doubling_time_hours",
    "DoublingTime",
    "doubling_time",
    "GrowthRate",
    "doubling_time_hrs",
]


def _require(path: Path, how: str) -> Path:
    if not path.exists():
        raise RawMissingError(path, "DepMap 26Q1", how)
    return path


def load_model_metadata(raw_dir: Path | str = DEPMAP_RAW) -> pd.DataFrame:
    """Tidy model metadata. Output columns:
        ModelID, CellLineName, StrippedCellLineName, OncotreeLineage,
        OncotreeCode, doubling_time_hours (may be NaN if not in the release).

    Expected file: ``Model.csv`` (portal 'Model' download).
    """
    raw_dir = Path(raw_dir)
    path = _require(raw_dir / "Model.csv", "Download 'Model.csv' from depmap.org/portal/download")
    df = pd.read_csv(path)
    dt_col = next((c for c in _DOUBLING_CANDIDATES if c in df.columns), None)
    out = pd.DataFrame({"ModelID": df["ModelID"]})
    for c in ("CellLineName", "StrippedCellLineName", "OncotreeLineage", "OncotreeCode"):
        out[c] = df[c] if c in df.columns else pd.NA
    out["doubling_time_hours"] = df[dt_col] if dt_col else pd.NA

    if dt_col is None:
        # Optional separate doubling-time table shipped alongside Model.csv.
        dt_path = raw_dir / "DoublingTime.csv"
        if dt_path.exists():
            dt = pd.read_csv(dt_path)
            key = next((c for c in _DOUBLING_CANDIDATES if c in dt.columns), None)
            if key and "ModelID" in dt.columns:
                out = out.drop(columns="doubling_time_hours").merge(
                    dt[["ModelID", key]].rename(columns={key: "doubling_time_hours"}),
                    on="ModelID", how="left",
                )
    return out


def load_growth_rate(raw_dir: Path | str = DEPMAP_RAW) -> pd.DataFrame:
    """Screen-inferred relative growth rate per model (a proliferation proxy).

    Expected file: ``CRISPRInferredModelGrowthRate.csv`` (columns per screen
    library). Returns ModelID + ``growth`` (mean over available libraries).
    """
    raw_dir = Path(raw_dir)
    path = _require(raw_dir / "CRISPRInferredModelGrowthRate.csv",
                    "Download 'CRISPRInferredModelGrowthRate.csv' from the DepMap release")
    g = pd.read_csv(path)
    lib_cols = [c for c in g.columns if c != "ModelID"]
    g["growth"] = g[lib_cols].mean(axis=1)
    return g[["ModelID", "growth"]].dropna(subset=["growth"])


def build_growth_model_meta(raw_dir: Path | str = DEPMAP_RAW) -> pd.DataFrame:
    """Model metadata joined to growth rate, with a doubling-time PROXY.

    ``doubling_time_hours`` here is ``1 / growth`` (arbitrary units): higher =
    slower-growing, so WARHEAD's G2b sign convention (positive slope = loses
    potency in slow lines = mitotic-dependent) applies unchanged. It is a proxy
    for real doubling time, not a measured value.
    Carries SangerModelID / COSMICID for joining GDSC.
    """
    raw_dir = Path(raw_dir)
    m = pd.read_csv(_require(raw_dir / "Model.csv", "Download 'Model.csv' from the DepMap release"))
    keep = [c for c in ["ModelID", "CellLineName", "StrippedCellLineName", "OncotreeCode",
                        "OncotreeLineage", "SangerModelID", "COSMICID"] if c in m.columns]
    meta = m[keep].merge(load_growth_rate(raw_dir), on="ModelID")
    meta["doubling_time_hours"] = 1.0 / meta["growth"]
    return meta
    """Long expression frame: ModelID, gene, expression (log2 TPM+1).

    Expected file: ``OmicsExpressionProteinCodingGenesTPMLogp1.csv`` (wide:
    rows = ModelID, columns = 'SYMBOL (ENTREZ)'). Pass ``genes`` (e.g.
    ['ABCB1', 'ABCG2']) to subset before melting.
    """
    raw_dir = Path(raw_dir)
    path = _require(
        raw_dir / "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
        "Download the Omics expression matrix from depmap.org/portal/download",
    )
    wide = pd.read_csv(path, index_col=0)
    wide.columns = [c.split(" (")[0] for c in wide.columns]
    if genes:
        keep = [g for g in genes if g in wide.columns]
        wide = wide[keep]
    long = wide.reset_index(names="ModelID").melt(
        id_vars="ModelID", var_name="gene", value_name="expression"
    )
    return long.dropna(subset=["expression"])


def load_crispr(raw_dir: Path | str = DEPMAP_RAW, genes: list[str] | None = None) -> pd.DataFrame:
    """Long CRISPR (Chronos) frame: ModelID, gene, chronos.

    Expected file: ``CRISPRGeneEffect.csv`` (wide: rows = ModelID, columns =
    'SYMBOL (ENTREZ)'). More negative = more essential.
    """
    raw_dir = Path(raw_dir)
    path = _require(
        raw_dir / "CRISPRGeneEffect.csv",
        "Download 'CRISPRGeneEffect.csv' from depmap.org/portal/download",
    )
    wide = pd.read_csv(path, index_col=0)
    wide.columns = [c.split(" (")[0] for c in wide.columns]
    if genes:
        keep = [g for g in genes if g in wide.columns]
        wide = wide[keep]
    long = wide.reset_index(names="ModelID").melt(
        id_vars="ModelID", var_name="gene", value_name="chronos"
    )
    return long.dropna(subset=["chronos"])
