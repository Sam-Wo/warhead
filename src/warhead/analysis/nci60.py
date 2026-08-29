"""NCI-60 / DTP drug-screen analysis (CellMiner z-score activity).

Two hard constraints shape what NCI-60 can answer here:
  * the public processed file is a z-score matrix (activity standardised PER
    COMPOUND across the 60 lines), so it gives RELATIVE selectivity but NOT
    absolute GI50/EC90 - there is no potency ranking or dose-response curve
    without the raw -logGI50 dataset.
  * the panel has NO liver line (so no HCC) and only 7 colon lines (CRC
    selectivity is over the huge >25k-compound library but statistically thin).

So the NCI-60 output is a CRC (colon)-selectivity ranking: which compounds are
more active in the colon lines than the rest, higher activity z = more sensitive.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..config import DATA_RAW

NCI60_XLSX = DATA_RAW / "nci60" / "output" / "DTP_NCI60_ZSCORE.xlsx"


def load_crc_selectivity(xlsx: Path | str = NCI60_XLSX, *, min_colon: int = 5) -> pd.DataFrame:
    """Per-compound mean activity z-score in colon (CRC) lines vs the rest."""
    d = pd.read_excel(xlsx, sheet_name="all", header=8)
    d = d.rename(columns={d.columns[0]: "NSC", d.columns[1]: "drug",
                          d.columns[2]: "FDA", d.columns[3]: "MOA"})
    cell_cols = [c for c in d.columns if ":" in str(c)]
    colon = [c for c in cell_cols if str(c).startswith("CO:")]
    other = [c for c in cell_cols if c in cell_cols and not str(c).startswith("CO:")]
    zc = d[colon].apply(pd.to_numeric, errors="coerce")
    zo = d[other].apply(pd.to_numeric, errors="coerce")
    out = pd.DataFrame({
        "drug": d["drug"], "NSC": d["NSC"], "FDA": d["FDA"], "MOA": d["MOA"],
        "n_colon": zc.notna().sum(axis=1),
        "mean_z_colon": zc.mean(axis=1), "mean_z_other": zo.mean(axis=1),
    })
    out["delta_z"] = out["mean_z_colon"] - out["mean_z_other"]
    out = out[(out["n_colon"] >= min_colon)].dropna(subset=["delta_z"])
    return out.sort_values("delta_z", ascending=False).reset_index(drop=True)


def _annotated(df: pd.DataFrame) -> pd.DataFrame:
    """Rows with a known mechanism or clinical status (interpretable hits)."""
    known = (df["MOA"].astype(str).str.strip().ne("-") & df["MOA"].notna()) | \
            (df["FDA"].astype(str).str.strip().ne("-") & df["FDA"].notna())
    return df[known]
