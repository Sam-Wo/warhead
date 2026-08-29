"""Novartis PDX Encyclopedia (Gao et al. 2015) - in-vivo drug response.

A DIFFERENT modality from the cell-line screens: PDXE measures tumour-volume
response of patient-derived xenografts (1x1x1 design), so there is no IC50 / EC90
/ dose-response curve. The analog outputs are (a) which treatments shrink CRC PDX
tumours most, and (b) whether a treatment is more effective in CRC than in the
other tumour types. PDXE has NO liver/HCC arm (tumour types: BRCA, CRC, CM, GC,
NSCLC, PDAC), so HCC cannot be assessed here.

BestAvgResponse = best average % change in tumour volume (negative = shrinkage).
ResponseCategory ~ RECIST (CR/PR = objective response, SD, PD).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..config import DATA_INTERIM

METRICS_PKL = DATA_INTERIM / "pdxe_metrics.pkl"


def load_metrics(pkl: Path | str = METRICS_PKL) -> pd.DataFrame:
    return pd.read_pickle(pkl)


def _objective(cat: pd.Series) -> pd.Series:
    return cat.astype(str).str.startswith(("CR", "PR"))


def crc_response_ranking(m: pd.DataFrame, *, single_only: bool = True) -> pd.DataFrame:
    """CRC treatments ranked by median tumour response (most shrinkage first)."""
    crc = m[m["Tumor Type"] == "CRC"]
    if single_only:
        crc = crc[crc["Treatment type"] == "single"]

    def _agg(g: pd.DataFrame) -> pd.Series:
        return pd.Series({
            "target": g["Treatment target"].iloc[0],
            "n_models": len(g),
            "median_response": float(np.median(g["BestAvgResponse"])),
            "pct_objective_response": float(_objective(g["ResponseCategory"]).mean() * 100),
        })

    return (crc.groupby("Treatment").apply(_agg, include_groups=False)
            .reset_index().sort_values("median_response").reset_index(drop=True))


def crc_response_selectivity(m: pd.DataFrame, *, single_only: bool = True, min_models: int = 5) -> pd.DataFrame:
    """Per treatment: is CRC response better (more shrinkage) than other tumour
    types? ``delta`` = median response in other types - CRC (positive = CRC more
    responsive)."""
    work = m[m["Treatment type"] == "single"] if single_only else m
    rows = []
    for tx, g in work.groupby("Treatment"):
        crc = g.loc[g["Tumor Type"] == "CRC", "BestAvgResponse"].to_numpy()
        oth = g.loc[g["Tumor Type"] != "CRC", "BestAvgResponse"].to_numpy()
        if crc.size < min_models or oth.size < min_models:
            continue
        rows.append({
            "Treatment": tx, "target": g["Treatment target"].iloc[0],
            "n_crc": int(crc.size), "median_response_crc": float(np.median(crc)),
            "median_response_other": float(np.median(oth)),
            "delta": float(np.median(oth) - np.median(crc)),
            "pct_objective_crc": float(_objective(g.loc[g["Tumor Type"] == "CRC", "ResponseCategory"]).mean() * 100),
        })
    return pd.DataFrame(rows).sort_values("delta", ascending=False).reset_index(drop=True)
