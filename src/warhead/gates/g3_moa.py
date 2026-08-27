"""G3 - MOA novelty and orthogonality. (Deferred: build order steps 5-6.)

G3a  Joint embedding of LINCS L1000 + Tahoe pseudobulk + JUMP morphology;
     score distance to nearest known-payload centroid CONDITIONAL on passing G1
     (novelty without potency is worthless). Normalise JUMP by dose first, or
     "novel MOA" collapses to "tested at a different concentration".
G3b  Orthogonal-resistance search - the exatecan-partner question. Regress out
     the Top1i component (or use SLFN11 as covariate) and rank compounds by
     potency in the RESIDUAL space, stratified by ABCB1 status so the result is
     not just efflux escape.
"""
from __future__ import annotations

import pandas as pd

from ..config import load_gates
from ..stats import weighted_linregress


def orthogonal_resistance_rank(
    sensitivity_matrix: pd.DataFrame,
    slfn11: pd.Series,
    *,
    config: dict | None = None,
) -> pd.DataFrame:
    """G3b (skeleton). ``sensitivity_matrix``: lines x compounds. ``slfn11``:
    SLFN11 expression indexed by line. Regress SLFN11 out of each compound's
    sensitivity and rank by mean residual potency on the SLFN11-low (Top1i-
    resistant) lines.

    Wired but not yet validated against ADCdb Top1i controls; kept minimal until
    real PRISM/CTRP data is local (build order step 5).
    """
    cfg = (config or load_gates())["g3"]["orthogonal_resistance"]
    common = sensitivity_matrix.index.intersection(slfn11.index)
    x = slfn11.loc[common].to_numpy(float)
    out = []
    for comp in sensitivity_matrix.columns:
        y = sensitivity_matrix.loc[common, comp].to_numpy(float)
        fit = weighted_linregress(x, y)
        resid_mean = float((y - (fit.intercept + fit.slope * x)).mean()) if fit.n else float("nan")
        out.append({"compound_id": comp, "slfn11_slope": fit.slope, "residual_potency": resid_mean})
    return (
        pd.DataFrame(out)
        .sort_values("residual_potency", ascending=False)
        .head(cfg["top_n"])
        .reset_index(drop=True)
    )
