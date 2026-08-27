"""4-parameter-logistic refit with an interval-censored (Tobit) likelihood.

Why not AUC / IC50 off the shelf (WARHEAD.md G1a):
  * AUC and IC50 conflate potency with completeness of kill. For a payload those
    are not interchangeable: an ADC cannot out-dose a persister fraction, so Emax
    matters independently of IC50. We therefore extract IC50, Emax and Hill as
    SEPARATE features and never collapse them.
  * Screens bottom out at 1-10 nM, exactly where payload-relevant activity
    begins. A reading pinned at the assay floor ("full response at lowest dose")
    or ceiling ("no response at highest dose") is a CENSORED observation, not a
    point estimate. Imputing the lowest tested dose systematically flattens the
    ranking. We fit those points with a Tobit likelihood (Normal density for
    interior points, Normal CDF mass for censored ones) and, when the fitted
    IC50 falls outside the tested range, we REPORT it as left/right censored
    rather than clamping it.

Model (viability normalised so 1 = untreated, 0 = complete kill):

    v(d) = emax + (top - emax) / (1 + (d / ic50) ** hill)

with hill > 0, so v -> top as d -> 0 and v -> emax as d -> inf.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Sequence

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import log_ndtr

_LOG_SQRT_2PI = 0.5 * np.log(2.0 * np.pi)


@dataclass
class FourPLFit:
    """Result of a single-curve refit. Concentrations in MOLAR."""

    ic50_M: float
    log10_ic50_M: float
    emax: float
    top: float
    hill: float
    sigma: float
    n_points: int
    converged: bool
    ic50_censoring: str  # 'left' | 'none' | 'right'
    rmse: float
    neg_loglik: float

    def as_dict(self) -> dict:
        return asdict(self)


def four_pl(dose_M: np.ndarray, top: float, emax: float, log10_ic50: float, hill: float) -> np.ndarray:
    """Vectorised 4PL. ``dose_M`` in molar; ``log10_ic50`` is log10(IC50 [M])."""
    dose_M = np.asarray(dose_M, dtype=float)
    log10_d = np.log10(np.clip(dose_M, 1e-30, None))
    # (d/ic50)^hill == 10 ** (hill * (log10 d - log10 ic50)); stable in log space.
    ratio = np.power(10.0, hill * (log10_d - log10_ic50))
    return emax + (top - emax) / (1.0 + ratio)


def _initial_guess(log10_d: np.ndarray, y: np.ndarray) -> tuple[float, float, float, float]:
    top0 = float(np.clip(np.nanmax(y), 0.6, 1.2))
    emax0 = float(np.clip(np.nanmin(y), -0.1, top0 - 0.05))
    mid = 0.5 * (top0 + emax0)
    # First dose whose response drops below the midpoint -> rough IC50.
    order = np.argsort(log10_d)
    ld, ys = log10_d[order], y[order]
    below = np.where(ys <= mid)[0]
    log_ic50_0 = float(ld[below[0]]) if below.size else float(np.median(ld))
    return top0, emax0, log_ic50_0, 1.0


def refit_curve(
    dose_M: Sequence[float],
    viability: Sequence[float],
    *,
    response_floor: float = 0.02,
    response_ceiling: float = 1.05,
    hill_bounds: tuple[float, float] = (0.3, 8.0),
    emax_bounds: tuple[float, float] = (-0.2, 1.2),
    top_bounds: tuple[float, float] = (0.6, 1.3),
    ic50_pad_logs: float = 3.0,
) -> FourPLFit:
    """Refit one dose-response curve with the Tobit 4PL likelihood.

    Points with ``viability <= response_floor`` are treated as left-censored in
    the readout (true value somewhere at/below the floor) and points with
    ``viability >= response_ceiling`` as right-censored. Interior points use a
    Normal density. The IC50 bound is padded ``ic50_pad_logs`` beyond the tested
    dose range so the optimiser can place IC50 outside the range; when it does,
    we flag it censored instead of pretending we measured it.
    """
    dose_M = np.asarray(dose_M, dtype=float)
    y = np.asarray(viability, dtype=float)
    ok = np.isfinite(dose_M) & (dose_M > 0) & np.isfinite(y)
    dose_M, y = dose_M[ok], y[ok]
    n = y.size
    if n == 0:
        return FourPLFit(np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 0, False, "none", np.nan, np.nan)

    log10_d = np.log10(dose_M)
    dmin, dmax = float(log10_d.min()), float(log10_d.max())

    left = y <= response_floor
    right = y >= response_ceiling
    interior = ~(left | right)

    top0, emax0, lic0, hill0 = _initial_guess(log10_d, y)
    resid0 = y[interior] - four_pl(dose_M[interior], top0, emax0, lic0, hill0) if interior.any() else np.array([0.05])
    logsig0 = float(np.log(np.clip(np.std(resid0) if resid0.size else 0.08, 0.02, 0.5)))
    x0 = np.array([top0, emax0, lic0, hill0, logsig0])

    lic_lo, lic_hi = dmin - ic50_pad_logs, dmax + ic50_pad_logs
    bounds = [top_bounds, emax_bounds, (lic_lo, lic_hi), hill_bounds, (np.log(1e-3), np.log(1.0))]

    def nll(theta: np.ndarray) -> float:
        top, emax, lic, hill, logsig = theta
        sigma = np.exp(logsig)
        v = four_pl(dose_M, top, emax, lic, hill)
        total = 0.0
        if interior.any():
            r = (y[interior] - v[interior]) / sigma
            total += np.sum(_LOG_SQRT_2PI + logsig + 0.5 * r * r)
        if left.any():
            total += -np.sum(log_ndtr((response_floor - v[left]) / sigma))
        if right.any():
            total += -np.sum(log_ndtr((v[right] - response_ceiling) / sigma))
        return float(total)

    res = minimize(nll, x0, method="L-BFGS-B", bounds=bounds)
    top, emax, lic, hill, logsig = res.x
    sigma = float(np.exp(logsig))

    # Censoring verdict on the DERIVED IC50: report, do not clamp. IC50 is only
    # identified when the tested dose range brackets the half-maximal response.
    # When it does not (or the curve is flat because there is no transition
    # in-range), the potency is censored - which direction is set by whether the
    # cells were mostly killed (left) or mostly untouched (right).
    resp_range = top - emax
    halfmax = 0.5 * (top + emax)
    v_lo = float(four_pl(np.array([10.0 ** dmin]), top, emax, lic, hill)[0])
    v_hi = float(four_pl(np.array([10.0 ** dmax]), top, emax, lic, hill)[0])
    if resp_range < 0.15:
        censoring = "right" if float(np.median(y)) > 0.5 else "left"
    elif v_hi > halfmax + 0.05 or lic >= dmax:
        censoring = "right"          # half-effect never reached by the top dose
    elif v_lo < halfmax - 0.05 or lic <= dmin:
        censoring = "left"           # already past half-effect at the bottom dose
    else:
        censoring = "none"

    v_all = four_pl(dose_M, top, emax, lic, hill)
    rmse = float(np.sqrt(np.mean((y - v_all) ** 2)))

    return FourPLFit(
        ic50_M=float(10.0 ** lic),
        log10_ic50_M=float(lic),
        emax=float(emax),
        top=float(top),
        hill=float(hill),
        sigma=sigma,
        n_points=int(n),
        converged=bool(res.success),
        ic50_censoring=censoring,
        rmse=rmse,
        neg_loglik=float(res.fun),
    )


def refit_frame(
    long_df: pd.DataFrame,
    *,
    compound_col: str = "compound_id",
    model_col: str = "ModelID",
    dose_col: str = "dose_M",
    viability_col: str = "viability",
    min_points: int = 5,
    **refit_kwargs,
) -> pd.DataFrame:
    """Refit every (compound x cell line) curve in a tidy long frame.

    Returns one row per curve with IC50/Emax/Hill and censoring flags. Curves
    with fewer than ``min_points`` usable points are dropped (recorded in QC),
    never silently imputed.
    """
    # 'model' is release metadata in gates.yaml's g1.refit block, not a curve arg.
    refit_kwargs.pop("model", None)
    rows = []
    keys = [compound_col, model_col]
    for (comp, model), g in long_df.groupby(keys, sort=False):
        d = g[dose_col].to_numpy(dtype=float)
        v = g[viability_col].to_numpy(dtype=float)
        usable = np.isfinite(d) & (d > 0) & np.isfinite(v)
        if int(usable.sum()) < min_points:
            continue
        fit = refit_curve(d[usable], v[usable], **refit_kwargs)
        rows.append({compound_col: comp, model_col: model, **fit.as_dict()})
    cols = [
        compound_col, model_col, "ic50_M", "log10_ic50_M", "emax", "top", "hill",
        "sigma", "n_points", "converged", "ic50_censoring", "rmse", "neg_loglik",
    ]
    return pd.DataFrame(rows, columns=cols)
