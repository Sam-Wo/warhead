"""GDSC (Genomics of Drug Sensitivity in Cancer) loader - release 8.5.

Source: cancerrxgene.org bulk download (Sanger COG). We use the fitted
dose-response table (LN_IC50 + AUC per drug x cell line) plus tissue and compound
annotation. GDSC1/GDSC2 must NOT be pooled naively for overlapping compounds
(WARHEAD.md sec 2) - load one dataset at a time.

Expected files in ``data/raw/gdsc/``:
  * GDSC2_fitted_dose_response.xlsx   (or GDSC1_...)
  * Cell_Lines_Details.xlsx           (optional; TCGA_DESC is already in fitted)
  * screened_compounds.csv            (optional; target/pathway already in fitted)

EC90 note (WARHEAD.md G1a): GDSC's public fit is a 2-parameter sigmoid in
ln-concentration with the bottom fixed at 0 (it assumes complete kill). We recover
the slope ``scal`` from (AUC, LN_IC50, tested range) and report
``ec90_uM = IC50 * 9**scal``. Two flags travel with it: ``ec90_range`` (within the
tested range vs extrapolated beyond MAX_CONC) and ``ec90_confidence`` (low when the
IC50 sits mid-range, where AUC ~ 0.5 for any slope so scal is weakly identified).
The reliable quantities are the fitted IC50 and the potency ORDERING; the EC90
magnitude is a model extrapolation. For an Emax-aware EC90 refit the raw
well-level data (the 2 GB file).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq

from ..config import DATA_RAW
from .base import RawMissingError

GDSC_RAW = DATA_RAW / "gdsc"

# WARHEAD indication <- GDSC TCGA_DESC.
INDICATION_TCGA = {"CRC": "COREAD", "HCC": "LIHC"}

_COLS = ["DATASET", "COSMIC_ID", "CELL_LINE_NAME", "SANGER_MODEL_ID", "TCGA_DESC",
         "DRUG_ID", "DRUG_NAME", "PUTATIVE_TARGET", "PATHWAY_NAME",
         "MIN_CONC", "MAX_CONC", "LN_IC50", "AUC", "RMSE", "Z_SCORE"]


def load_fitted(raw_dir: Path | str = GDSC_RAW, dataset: str = "GDSC2") -> pd.DataFrame:
    """Load the fitted dose-response table as a tidy frame (one row per curve).

    Adds ``ic50_uM`` = exp(LN_IC50). Concentrations are micromolar.
    """
    raw_dir = Path(raw_dir)
    path = raw_dir / f"{dataset}_fitted_dose_response.xlsx"
    if not path.exists():
        raise RawMissingError(
            path, f"GDSC {dataset}",
            "Download from cancerrxgene.org bulk download "
            "(cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/)",
        )
    raw = pd.read_excel(path)
    cols = [c for c in _COLS if c in raw.columns]
    df = raw[cols].rename(columns={
        "CELL_LINE_NAME": "cell_line", "SANGER_MODEL_ID": "sanger_model_id",
        "TCGA_DESC": "tcga_desc", "DRUG_NAME": "drug_name", "DRUG_ID": "drug_id",
        "PUTATIVE_TARGET": "target", "PATHWAY_NAME": "pathway",
        "MIN_CONC": "min_conc_uM", "MAX_CONC": "max_conc_uM",
        "LN_IC50": "ln_ic50", "AUC": "auc", "RMSE": "rmse", "Z_SCORE": "z_score",
    })
    df["ic50_uM"] = np.exp(df["ln_ic50"])
    return df


# ---------------------------------------------------------------------------
# EC90 recovery from the GDSC 2-parameter fit
# ---------------------------------------------------------------------------
def _Fint(u: np.ndarray) -> np.ndarray:
    """Stable antiderivative of 1/(1+e^u): u - log(1+e^u)."""
    u = np.asarray(u, float)
    return np.where(u > 0, -np.log1p(np.exp(-np.abs(u))), u - np.log1p(np.exp(-np.abs(u))))


def _mean_viability(scal: float, xmid: float, xlo: float, xhi: float) -> float:
    """Mean fitted viability over [xlo, xhi] in ln-concentration."""
    umax = (xhi - xmid) / scal
    umin = (xlo - xmid) / scal
    return float(scal * (_Fint(umax) - _Fint(umin)) / (xhi - xlo))


def _solve_scal(auc: float, xmid: float, xlo: float, xhi: float) -> float:
    """Recover the slope that reproduces the reported AUC (mean viability)."""
    if not np.isfinite(auc) or xhi <= xlo:
        return np.nan
    f = lambda s: _mean_viability(s, xmid, xlo, xhi) - auc
    try:
        if f(0.02) * f(20.0) > 0:
            return np.nan
        return brentq(f, 0.02, 20.0, maxiter=200)
    except Exception:
        return np.nan


def add_ec90(df: pd.DataFrame) -> pd.DataFrame:
    """Add ``scal``, ``ec90_uM`` and ``ec90_range`` (interpolated/extrapolated).

    EC90 = IC50 * 9**scal (concentration for 90% of the fitted maximal effect).
    ``ec90_range`` is 'within' when EC90 <= MAX_CONC (observed), else 'extrapolated'.
    """
    out = df.copy()
    xlo = np.log(out["min_conc_uM"].to_numpy())
    xhi = np.log(out["max_conc_uM"].to_numpy())
    xmid = out["ln_ic50"].to_numpy()
    auc = out["auc"].to_numpy()
    scal = np.array([_solve_scal(auc[i], xmid[i], xlo[i], xhi[i]) for i in range(len(out))])
    out["scal"] = scal
    out["ec90_uM"] = out["ic50_uM"] * np.power(9.0, scal)
    ln_ec90 = np.log(out["ec90_uM"].to_numpy())
    out["ec90_range"] = np.where(ln_ec90 <= xhi, "within", "extrapolated")
    # IC50 outside the tested range is itself censored - record it too.
    out["ic50_range"] = np.where(xmid > xhi, "above_max",
                                 np.where(xmid < xlo, "below_min", "within"))
    # Identifiability: when IC50 sits mid-range, AUC ~ 0.5 for ANY slope, so scal
    # (hence EC90 magnitude) is weakly constrained. Flag those + boundary solves.
    at_bound = (scal <= 0.05) | (scal >= 19.0)
    centred = (out["ic50_range"].to_numpy() == "within") & (np.abs(auc - 0.5) < 0.04)
    out["ec90_confidence"] = np.where(at_bound | centred | ~np.isfinite(scal), "low", "ok")
    return out


def load_with_ec90(raw_dir: Path | str = GDSC_RAW, dataset: str = "GDSC2") -> pd.DataFrame:
    return add_ec90(load_fitted(raw_dir, dataset))
