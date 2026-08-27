"""G3 - MOA novelty and orthogonality.

G3a  Joint embedding of LINCS L1000 + Tahoe pseudobulk + JUMP morphology; score
     distance to nearest known-payload centroid CONDITIONAL on passing G1
     (novelty without potency is worthless). Normalise JUMP by dose first, or
     "novel MOA" collapses to "tested at a different concentration". [DEFERRED -
     needs Tahoe/JUMP local; build order step 6.]

G3b  Orthogonal-resistance search - the exatecan-partner question. Implemented
     below. Regress out the Top1i axis (SLFN11) and rank compounds by potency in
     the RESIDUAL space, i.e. on the lines exatecan does not handle, controlling
     for ABCB1 so the result is not just efflux escape. Current lead partner is
     ATRi; this generates partner candidates empirically to corroborate that or
     surface alternatives.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from ..config import load_gates
from ..stats import benjamini_hochberg, weighted_linregress


def sensitivity_matrix_from_fits(
    sensitivity: pd.DataFrame,
    *,
    compound_col: str = "compound_id",
    model_col: str = "ModelID",
    value_col: str = "sensitivity",
) -> pd.DataFrame:
    """Pivot the long sensitivity frame to a lines x compounds matrix."""
    return sensitivity.pivot_table(index=model_col, columns=compound_col, values=value_col)


def _covariate_series(expression: pd.DataFrame, gene: str,
                      gene_col: str = "gene", model_col: str = "ModelID",
                      expr_col: str = "expression") -> pd.Series:
    sub = expression[expression[gene_col] == gene]
    return sub.set_index(model_col)[expr_col]


def exatecan_partner_search(
    sensitivity_matrix: pd.DataFrame,
    slfn11: pd.Series,
    abcb1: pd.Series,
    *,
    config: dict | None = None,
    resistant_quantile: float = 1.0 / 3.0,
    top_n: int | None = None,
) -> pd.DataFrame:
    """Rank compounds as orthogonal exatecan partners.

    For each compound (columns of ``sensitivity_matrix``; higher value = more
    potent, i.e. -log10 IC50):

      1. Remove the ABCB1 (efflux) component of its sensitivity pattern, so we do
         not merely rediscover efflux escape.
      2. In that efflux-controlled residual space, regress on SLFN11. A positive
         ``slfn11_slope`` means the compound needs SLFN11 like a Top1i does (a
         BAD partner); ~0 or negative means its potency is orthogonal to the
         Top1i axis.
      3. Score ``orthogonality`` = mean efflux-controlled residual potency on the
         SLFN11-low (Top1i-resistant) lines - how much it over-performs exactly
         where exatecan fails.

    Ranked by ``orthogonality`` descending. ``resistant_potency`` reports the raw
    potency on the resistant subset (a partner must also be potent in absolute
    terms - feed G1 passers only).
    """
    cfg = (config or load_gates())["g3"]["orthogonal_resistance"]
    top_n = top_n or cfg["top_n"]

    idx = sensitivity_matrix.index
    slfn = slfn11.reindex(idx).astype(float)
    abc = abcb1.reindex(idx).astype(float)
    resistant_thr = slfn.quantile(resistant_quantile)

    rows = []
    for comp in sensitivity_matrix.columns:
        y = sensitivity_matrix[comp].astype(float)
        m = y.notna() & slfn.notna() & abc.notna()
        if int(m.sum()) < 10:
            continue
        yy = y[m].to_numpy()
        sz = (slfn[m] - slfn[m].mean()) / (slfn[m].std() + 1e-9)
        az = (abc[m] - abc[m].mean()) / (abc[m].std() + 1e-9)
        resistant = (slfn[m] <= resistant_thr).to_numpy()
        if int(resistant.sum()) < 3:
            continue

        # 1. strip efflux (ABCB1) component
        fa = weighted_linregress(az.to_numpy(), yy)
        resid = yy - (fa.intercept + fa.slope * az.to_numpy())
        # 2. SLFN11 dependence in the efflux-controlled space
        fs = weighted_linregress(sz.to_numpy(), resid)
        # 3. orthogonal residual potency on the Top1i-resistant lines
        orthogonality = float(resid[resistant].mean())
        resistant_potency = float(yy[resistant].mean())

        rows.append(
            {
                "compound_id": comp,
                "slfn11_slope": fs.std_slope,   # >0 = Top1i-like (bad partner)
                "slfn11_p": fs.p,
                "abcb1_slope": fa.std_slope,
                "orthogonality": orthogonality,
                "resistant_potency": resistant_potency,
                "n_lines": int(m.sum()),
                "n_resistant": int(resistant.sum()),
            }
        )

    out = pd.DataFrame(rows)
    if len(out):
        out["slfn11_q"] = benjamini_hochberg(out["slfn11_p"].to_numpy())
        # A partner: over-performs on resistant lines AND not SLFN11-dependent.
        out["is_partner_candidate"] = (out["orthogonality"] > 0) & (
            out["slfn11_slope"] <= cfg.get("partner_slfn11_slope_max", 0.1)
        )
        out = out.sort_values("orthogonality", ascending=False).reset_index(drop=True)
        out["rank"] = np.arange(1, len(out) + 1)
    return out.head(top_n)


def run_partner_search_from_expression(
    sensitivity: pd.DataFrame,
    expression: pd.DataFrame,
    *,
    config: dict | None = None,
    **kwargs,
) -> pd.DataFrame:
    """Convenience wrapper: build the sensitivity matrix and pull SLFN11 / ABCB1
    from a long expression frame (real or synthetic - same schema)."""
    matrix = sensitivity_matrix_from_fits(sensitivity)
    slfn11 = _covariate_series(expression, "SLFN11")
    abcb1 = _covariate_series(expression, "ABCB1")
    return exatecan_partner_search(matrix, slfn11, abcb1, config=config, **kwargs)
