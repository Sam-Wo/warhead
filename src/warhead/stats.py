"""Small statistics helpers implemented on numpy/scipy.

Kept dependency-light on purpose: the core cascade must not require statsmodels.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats as sstats


@dataclass
class SlopeFit:
    slope: float
    intercept: float
    slope_se: float
    t: float
    p: float          # two-sided, H0: slope == 0
    n: int
    r: float          # weighted Pearson correlation of x, y
    std_slope: float  # slope on standardised x, y (comparable across compounds)


def weighted_linregress(x: np.ndarray, y: np.ndarray, w: np.ndarray | None = None) -> SlopeFit:
    """Weighted OLS of y on x with a t-test on the slope.

    Weights are analytic (reliability) weights. Returns a standardised slope
    (weighted-z-scored x and y) so magnitudes are comparable across compounds
    with different sensitivity scales.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    w = np.ones_like(x) if w is None else np.asarray(w, dtype=float)[m]
    n = x.size
    if n < 3 or np.allclose(x, x[0]):
        return SlopeFit(np.nan, np.nan, np.nan, np.nan, np.nan, n, np.nan, np.nan)

    W = w.sum()
    xm = np.sum(w * x) / W
    ym = np.sum(w * y) / W
    dx = x - xm
    dy = y - ym
    Sxx = np.sum(w * dx * dx)
    Sxy = np.sum(w * dx * dy)
    Syy = np.sum(w * dy * dy)
    slope = Sxy / Sxx
    intercept = ym - slope * xm

    resid = y - (intercept + slope * x)
    dof = n - 2
    # Weighted residual variance.
    sigma2 = np.sum(w * resid * resid) / (dof if dof > 0 else 1)
    slope_se = np.sqrt(sigma2 / Sxx)
    t = slope / slope_se if slope_se > 0 else np.nan
    p = float(2.0 * sstats.t.sf(np.abs(t), dof)) if np.isfinite(t) and dof > 0 else np.nan
    r = Sxy / np.sqrt(Sxx * Syy) if Sxx > 0 and Syy > 0 else np.nan

    sx = np.sqrt(Sxx / W)
    sy = np.sqrt(Syy / W)
    std_slope = slope * (sx / sy) if sy > 0 else np.nan
    return SlopeFit(float(slope), float(intercept), float(slope_se), float(t),
                    p, int(n), float(r), float(std_slope))


def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (q-values). NaNs pass through."""
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    finite = np.isfinite(p)
    pv = p[finite]
    m = pv.size
    if m == 0:
        return out
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * m / (np.arange(1, m + 1))
    # Enforce monotonicity from the largest p downward.
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    adj = np.empty(m)
    adj[order] = q
    out[finite] = adj
    return out


def balance_weights(values: np.ndarray, n_bins: int = 5, cap: float = 5.0) -> np.ndarray:
    """Weights that equalise the influence of equal-WIDTH bins of ``values``.

    Corrects pooled-PRISM under-representation of slow-growing lines
    (WARHEAD.md G2b): binning doubling time into equal-width intervals means the
    sparsely populated slow bins get up-weighted, so the regression is not
    dominated by the over-sampled fast lines. (Equal-COUNT / quantile bins would
    not do this - by construction they hold equal numbers, so they cannot
    up-weight a rare region.)

    Weights are trimmed to ``[1/cap, cap]`` before renormalising. Uncapped
    inverse-frequency weights let a near-empty bin's one or two lines dominate the
    regression, which explodes the slope variance; capping keeps the correction
    while bounding any single line's leverage.
    """
    v = np.asarray(values, dtype=float)
    w = np.ones_like(v)
    finite = np.isfinite(v)
    if finite.sum() < n_bins:
        return w
    lo, hi = float(np.min(v[finite])), float(np.max(v[finite]))
    if hi <= lo:
        return w
    edges = np.linspace(lo, hi, n_bins + 1)
    edges[-1] = np.inf  # include the max in the top bin
    binidx = np.digitize(v, edges[1:-1])  # 0 .. n_bins-1
    occupied = np.unique(binidx[finite])
    target = finite.sum() / len(occupied)
    for b in occupied:
        sel = finite & (binidx == b)
        cnt = int(sel.sum())
        if cnt > 0:
            w[sel] = target / cnt
    # Normalise to mean 1, trim leverage, renormalise.
    w[finite] *= finite.sum() / w[finite].sum()
    w[finite] = np.clip(w[finite], 1.0 / cap, cap)
    w[finite] *= finite.sum() / w[finite].sum()
    return w
