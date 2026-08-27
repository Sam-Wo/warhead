"""reports/exatecan_partner.pdf - the standalone G3b deliverable for the
dual-payload program (WARHEAD.md sec 7): which payload class is potent on the
lines exatecan cannot handle, orthogonal to the Top1i/SLFN11 axis?"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from ..config import load_gates  # noqa: E402
from ..stats import weighted_linregress  # noqa: E402

RRB_MAROON = "#6E1426"
_PARTNER = "#2C7FB8"
_NOT = "#B5651D"


def _slfn11_series(expression: pd.DataFrame) -> pd.Series:
    s = expression[expression["gene"] == "SLFN11"].set_index("ModelID")["expression"]
    return s


def render_exatecan_report(
    g3b: pd.DataFrame,
    sensitivity: pd.DataFrame,
    expression: pd.DataFrame,
    *,
    out_path: str | Path,
    top1i_ref: str = "exatecan_like",
    config: dict | None = None,
) -> Path:
    cfg = config or load_gates()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    ranked = g3b.sort_values("orthogonality").reset_index(drop=True)
    slfn11 = _slfn11_series(expression)
    top_partner = g3b.sort_values("orthogonality", ascending=False).iloc[0]["compound_id"]

    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("WARHEAD - G3b Exatecan Partner Search", color=RRB_MAROON,
                 fontsize=16, fontweight="bold", x=0.5, y=0.97)
    fig.text(0.5, 0.935,
             "Potency orthogonal to the Top1i / SLFN11 axis  -  who kills where exatecan cannot",
             ha="center", fontsize=10, color="#444")

    gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 1.0], hspace=0.42, wspace=0.30,
                          left=0.16, right=0.96, top=0.90, bottom=0.09)

    # Panel A: orthogonality ranking.
    axA = fig.add_subplot(gs[0, :])
    colors = [_PARTNER if p else _NOT for p in ranked["is_partner_candidate"]]
    ypos = np.arange(len(ranked))
    axA.barh(ypos, ranked["orthogonality"], color=colors, edgecolor="#222", linewidth=0.4)
    axA.axvline(0, color="#222", lw=0.8)
    axA.set_yticks(ypos)
    axA.set_yticklabels(ranked["compound_id"], fontsize=8)
    axA.set_xlabel("orthogonality  =  efflux-controlled residual potency on SLFN11-low (Top1i-resistant) lines")
    axA.set_title("A.  Partner ranking (blue = orthogonal partner candidate; higher = better)",
                  fontsize=9, loc="left")

    # Panel B: the two axes - Top1i dependence vs orthogonal potency.
    axB = fig.add_subplot(gs[1, 0])
    axB.scatter(g3b["slfn11_slope"], g3b["orthogonality"], s=28,
                color=[_PARTNER if p else _NOT for p in g3b["is_partner_candidate"]],
                edgecolor="#222", linewidth=0.4, zorder=3)
    axB.axvline(0, color="#888", lw=0.7, ls=":")
    axB.axhline(0, color="#888", lw=0.7, ls=":")
    for _, r in g3b.iterrows():
        if r["compound_id"] in (top_partner, top1i_ref):
            axB.annotate(r["compound_id"], (r["slfn11_slope"], r["orthogonality"]),
                         fontsize=8, color=RRB_MAROON, xytext=(4, 4),
                         textcoords="offset points")
    axB.set_xlabel("SLFN11 dependence (>0 = Top1i-like)")
    axB.set_ylabel("orthogonality")
    axB.set_title("B.  Orthogonal potency vs Top1i dependence", fontsize=9, loc="left")

    # Panel C: sensitivity vs SLFN11 for the top partner and the Top1i reference.
    axC = fig.add_subplot(gs[1, 1])
    for comp, col in [(top_partner, _PARTNER), (top1i_ref, RRB_MAROON)]:
        sub = sensitivity[sensitivity["compound_id"] == comp].merge(
            slfn11.rename("slfn11"), left_on="ModelID", right_index=True, how="inner"
        )
        if not len(sub):
            continue
        x = sub["slfn11"].to_numpy(float)
        y = sub["sensitivity"].to_numpy(float)
        axC.scatter(x, y, s=14, color=col, alpha=0.6, edgecolor="none")
        if len(x) > 3:
            f = weighted_linregress(x, y)
            xs = np.linspace(x.min(), x.max(), 40)
            axC.plot(xs, f.intercept + f.slope * xs, color=col, lw=1.8,
                     label=f"{comp} (slope {f.std_slope:+.2f})")
    axC.set_xlabel("SLFN11 expression")
    axC.set_ylabel("potency  (-log10 IC50)")
    axC.set_title("C.  Partner is potent where SLFN11 is low", fontsize=9, loc="left")
    axC.legend(fontsize=7, loc="best", frameon=False)

    n_partner = int(g3b["is_partner_candidate"].sum())
    fig.text(0.16, 0.02,
             f"candidates = {len(g3b)}   |   partner candidates = {n_partner}"
             f"   |   top partner = {top_partner}   |   Top1i ref = {top1i_ref}",
             fontsize=8, color="#555")

    fmt = out_path.suffix.lstrip(".").lower() or "pdf"
    fig.savefig(out_path, format=fmt, dpi=150)
    plt.close(fig)
    return out_path
