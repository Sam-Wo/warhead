"""Refit QC. Flags dubious curves; never silently drops without a reason."""
from __future__ import annotations

import numpy as np
import pandas as pd


def qc_flags(fits: pd.DataFrame, *, max_rmse: float = 0.2, max_sigma: float = 0.25) -> pd.DataFrame:
    """Annotate a refit frame with QC flags and a single ``qc_pass`` boolean.

    Flags (any True -> qc_pass False):
      * ``qc_not_converged`` - optimiser did not report success
      * ``qc_high_rmse``     - poor fit to the observed points
      * ``qc_high_sigma``    - noisy curve
      * ``qc_flat``          - top and emax within noise (no real response)
    """
    out = fits.copy()
    out["qc_not_converged"] = ~out["converged"].astype(bool)
    out["qc_high_rmse"] = out["rmse"] > max_rmse
    out["qc_high_sigma"] = out["sigma"] > max_sigma
    out["qc_flat"] = (out["top"] - out["emax"]).abs() < (2.0 * out["sigma"].clip(lower=1e-6))
    flag_cols = ["qc_not_converged", "qc_high_rmse", "qc_high_sigma", "qc_flat"]
    out["qc_pass"] = ~out[flag_cols].any(axis=1)
    return out


def summarise_qc(fits_qc: pd.DataFrame) -> dict[str, float]:
    n = len(fits_qc)
    return {
        "n_curves": n,
        "qc_pass": int(fits_qc["qc_pass"].sum()),
        "qc_fail": int((~fits_qc["qc_pass"]).sum()),
        "left_censored": int((fits_qc["ic50_censoring"] == "left").sum()),
        "right_censored": int((fits_qc["ic50_censoring"] == "right").sum()),
        "pass_frac": float(fits_qc["qc_pass"].mean()) if n else 0.0,
    }
