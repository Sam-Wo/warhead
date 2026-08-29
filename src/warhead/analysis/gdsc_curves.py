"""Reconstruct real GDSC dose-response curves from the raw well-level data.

The fitted table only gives LN_IC50 + AUC; to SEE the measured curves (and where
IC50 / EC90 actually fall relative to the tested range) we normalise the raw
per-well fluorescence. Per plate (BARCODE):

    viability = (treated_intensity - blank) / (neg_control - blank)

where negative controls are the untreated/vehicle wells (DRUG_ID null, TAG ~ 'NC')
and blanks are the no-cell wells (DRUG_ID null, TAG 'B'). Treatment wells carry a
DRUG_ID + CONC. The 2 GB raw file is streamed in chunks; only the target drugs'
treatment wells are kept, plus running per-plate control/blank means.
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

_USECOLS = ["BARCODE", "CELL_LINE_NAME", "SANGER_MODEL_ID", "TAG", "DRUG_ID", "CONC", "INTENSITY"]


def _is_neg_control(tag: pd.Series) -> pd.Series:
    t = tag.astype(str).str.upper()
    return t.str.startswith("NC") | t.str.contains("DMSO")


def _is_blank(tag: pd.Series) -> pd.Series:
    t = tag.astype(str).str.upper()
    return (t == "B") | t.str.startswith("B-") | t.str.contains("BLANK")


def extract_raw_curves(raw_csv: str | Path, drug_ids: list[int], *,
                       chunksize: int = 3_000_000) -> pd.DataFrame:
    """Stream the raw CSV once; return normalised long curves for ``drug_ids``:
    columns [drug_id, cell_line, sanger_model_id, conc_uM, viability]."""
    raw_csv = Path(raw_csv)
    targets = set(int(x) for x in drug_ids)

    nc_sum: dict = defaultdict(float); nc_cnt: dict = defaultdict(int)
    b_sum: dict = defaultdict(float); b_cnt: dict = defaultdict(int)
    treat_parts = []

    for chunk in pd.read_csv(raw_csv, usecols=_USECOLS, chunksize=chunksize):
        drug = pd.to_numeric(chunk["DRUG_ID"], errors="coerce")
        is_treat = drug.isin(targets)
        no_drug = drug.isna()

        # accumulate per-plate negative-control and blank means (all plates)
        nc = chunk[no_drug & _is_neg_control(chunk["TAG"])]
        for bc, inten in zip(nc["BARCODE"].to_numpy(), nc["INTENSITY"].to_numpy()):
            nc_sum[bc] += inten; nc_cnt[bc] += 1
        bl = chunk[no_drug & _is_blank(chunk["TAG"])]
        for bc, inten in zip(bl["BARCODE"].to_numpy(), bl["INTENSITY"].to_numpy()):
            b_sum[bc] += inten; b_cnt[bc] += 1

        if is_treat.any():
            t = chunk[is_treat][["BARCODE", "CELL_LINE_NAME", "SANGER_MODEL_ID", "DRUG_ID", "CONC", "INTENSITY"]].copy()
            treat_parts.append(t)

    treat = pd.concat(treat_parts, ignore_index=True)
    nc_mean = {bc: nc_sum[bc] / nc_cnt[bc] for bc in nc_cnt}
    b_mean = {bc: b_sum[bc] / b_cnt[bc] for bc in b_cnt if b_cnt[bc]}

    treat["nc"] = treat["BARCODE"].map(nc_mean)
    treat["blank"] = treat["BARCODE"].map(b_mean).fillna(0.0)
    denom = treat["nc"] - treat["blank"]
    treat["viability"] = (treat["INTENSITY"] - treat["blank"]) / denom
    treat = treat[np.isfinite(treat["viability"]) & (denom > 0)]

    return treat.rename(columns={
        "CELL_LINE_NAME": "cell_line", "SANGER_MODEL_ID": "sanger_model_id",
        "DRUG_ID": "drug_id", "CONC": "conc_uM",
    })[["drug_id", "cell_line", "sanger_model_id", "conc_uM", "viability"]]


def pool_by_conc(curves: pd.DataFrame) -> pd.DataFrame:
    """Median viability + IQR across cell lines at each tested concentration."""
    g = curves.groupby(["drug_id", "conc_uM"])["viability"]
    out = g.agg(median="median",
                q1=lambda s: s.quantile(.25),
                q3=lambda s: s.quantile(.75),
                n="size").reset_index()
    return out.sort_values(["drug_id", "conc_uM"]).reset_index(drop=True)
