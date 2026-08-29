"""PRISM Repurposing secondary screen loader (Corsello et al. 2020).

The secondary screen fits a 4-parameter dose-response per compound x line and
publishes the parameters, so EC90 comes straight from the fit - and, unlike
GDSC, PRISM fits a FREE lower asymptote (``lower_limit`` = Emax), so it captures
incomplete kill. It also ships ``target``, ``moa`` and clinical ``phase`` inline.

    viability(d) = lower + (upper - lower) / (1 + (d/ec50)^slope)
    EC90 = ec50 * 9^(1/|slope|)     (concentration for 90% of the max effect)

Files (figshare article 9393293):
  * secondary-screen-dose-response-curve-parameters.csv   (the fit parameters)
Concentrations are micromolar; the secondary screen tops out near 10 uM.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..config import DATA_RAW
from .base import RawMissingError

PRISM_RAW = DATA_RAW / "prism"
PRISM_MAX_CONC_UM = 10.0   # secondary screen max tested dose (~10 uM)

_PARAM_COLS = ["depmap_id", "ccle_name", "screen_id", "upper_limit", "lower_limit",
               "slope", "r2", "auc", "ec50", "ic50", "name", "moa", "target", "phase"]


def load_curve_params(raw_dir: Path | str = PRISM_RAW) -> pd.DataFrame:
    raw_dir = Path(raw_dir)
    path = raw_dir / "secondary-screen-dose-response-curve-parameters.csv"
    if not path.exists():
        raise RawMissingError(path, "PRISM Repurposing secondary",
                              "Download from figshare article 9393293 "
                              "(ndownloader.figshare.com/files/20237739)")
    return pd.read_csv(path, usecols=_PARAM_COLS)


def to_canonical(params: pd.DataFrame, indication_map: dict[str, str]) -> pd.DataFrame:
    """Params -> the canonical screen frame (see analysis.screen_potency).

    ``indication_map``: depmap_id -> 'CRC' / 'HCC' / 'other'.
    """
    df = params.copy()
    # keep interpretable fits: sane EC50, real slope, plausible asymptotes
    df = df[np.isfinite(df["ec50"]) & (df["ec50"] > 1e-5) & (df["ec50"] < 1e3)
            & np.isfinite(df["slope"]) & (df["slope"].abs() > 1e-3)
            & df["upper_limit"].between(0.5, 1.5)]
    df["ec90_uM"] = df["ec50"] * np.power(9.0, 1.0 / df["slope"].abs())
    # a near-zero slope makes 9^(1/slope) explode; cap at 10 mM (beyond any assay)
    df = df[np.isfinite(df["ec90_uM"]) & (df["ec90_uM"] < 1e4)]

    # aggregate the (rare) duplicate screens per compound x line
    agg = (df.groupby(["name", "depmap_id"], as_index=False)
             .agg(target=("target", "first"), moa=("moa", "first"),
                  clinical_phase=("phase", "first"),
                  ic50_uM=("ic50", "median"), ec90_uM=("ec90_uM", "median"),
                  emax=("lower_limit", "median")))
    agg["indication"] = agg["depmap_id"].map(indication_map).fillna("other")
    agg["ic50_nM"] = agg["ic50_uM"] * 1e3
    agg["ec90_nM"] = agg["ec90_uM"] * 1e3
    agg["emax"] = agg["emax"].clip(lower=0, upper=1.2)
    agg["ec90_extrapolated"] = agg["ec90_uM"] > PRISM_MAX_CONC_UM
    agg["source"] = "PRISM"
    agg = agg.rename(columns={"name": "compound", "depmap_id": "model_id"})
    return agg[["source", "compound", "target", "moa", "model_id", "indication",
                "ic50_nM", "ec90_nM", "emax", "ec90_extrapolated", "clinical_phase"]]


def indication_map_from_model(model_csv: Path | str) -> dict[str, str]:
    """depmap_id -> CRC/HCC/other from a DepMap Model.csv (OncotreeCode)."""
    m = pd.read_csv(model_csv, usecols=["ModelID", "OncotreeCode"])
    # Real DepMap Oncotree codes: hepatocellular = HCC (not the TCGA "LIHC").
    code_to_ind = {"COAD": "CRC", "READ": "CRC", "COADREAD": "CRC",
                   "HCC": "HCC", "HCCIHCH": "HCC"}
    return {r.ModelID: code_to_ind.get(r.OncotreeCode, "other") for r in m.itertuples()}


def load_canonical(raw_dir: Path | str = PRISM_RAW,
                   model_csv: Path | str = DATA_RAW / "depmap" / "Model.csv") -> pd.DataFrame:
    return to_canonical(load_curve_params(raw_dir), indication_map_from_model(model_csv))
