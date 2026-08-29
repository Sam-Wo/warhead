"""CTRP v2 loader (via the base-R export of the Zenodo PharmacoGx PharmacoSet).

CTRP v2 has no working flat-file download any more; we read the Zenodo PharmacoSet
(.rds) with base R (see scripts/ctrp_export.R - no PharmacoGx needed) into
``data/interim/ctrp_export.csv``, then map to the canonical screen frame here.

CTRP is the widest-window screen (16-point, to ~66 uM) and fits a free lower
asymptote (``E_inf`` = residual viability %), so EC90 = EC50 * 9^(1/HS) is rarely
extrapolated and carries a real Emax. Targets and FDA status come inline.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..config import DATA_INTERIM
from .base import RawMissingError

CTRP_EXPORT = DATA_INTERIM / "ctrp_export.csv"
_SITE_INDICATION = {"large_intestine": "CRC", "liver": "HCC"}


def load_canonical(export_csv: Path | str = CTRP_EXPORT) -> pd.DataFrame:
    export_csv = Path(export_csv)
    if not export_csv.exists():
        raise RawMissingError(export_csv, "CTRP v2",
                              "Run scripts/ctrp_export.R (base R) on the Zenodo "
                              "CTRPv2.rds (record 3905470) to produce ctrp_export.csv")
    d = pd.read_csv(export_csv, low_memory=False)
    # keep interpretable fits: sane EC50 (drop the 1e6 non-fit sentinel) + real slope
    d = d[np.isfinite(d["ec50_uM"]) & (d["ec50_uM"] > 1e-5) & (d["ec50_uM"] < 1e4)
          & np.isfinite(d["HS"]) & (d["HS"] > 1e-3)]
    d["ec90_uM"] = d["ec50_uM"] * np.power(9.0, 1.0 / d["HS"])
    # a near-zero Hill slope makes 9^(1/HS) explode; cap at 10 mM (beyond any assay)
    d = d[np.isfinite(d["ec90_uM"]) & (d["ec90_uM"] < 1e4)]
    d["emax"] = (d["E_inf"] / 100.0).clip(lower=0, upper=1.2)
    d["indication"] = d["primary_site"].map(_SITE_INDICATION).fillna("other")
    fda_bool = d["fda"].map(lambda v: v is True or str(v).strip().lower() in ("true", "1"))
    d["clinical_phase"] = np.where(fda_bool, "FDA approved", "")

    agg = (d.groupby(["drug", "cellid"], as_index=False)
             .agg(target=("target", "first"), moa=("moa", "first"),
                  clinical_phase=("clinical_phase", "first"), indication=("indication", "first"),
                  ic50_uM=("ic50_uM", "median"), ec90_uM=("ec90_uM", "median"),
                  emax=("emax", "median"), max_conc_uM=("max_conc_uM", "median")))
    agg["source"] = "CTRP v2"
    agg["ic50_nM"] = agg["ic50_uM"] * 1e3
    agg["ec90_nM"] = agg["ec90_uM"] * 1e3
    agg["ec90_extrapolated"] = agg["ec90_uM"] > agg["max_conc_uM"]
    agg = agg.rename(columns={"drug": "compound", "cellid": "model_id"})
    return agg[["source", "compound", "target", "moa", "model_id", "indication",
                "ic50_nM", "ec90_nM", "emax", "ec90_extrapolated", "clinical_phase"]]
