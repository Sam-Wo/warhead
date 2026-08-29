"""reports/pdxe_crc_response.pdf - Novartis PDXE in-vivo CRC drug response.

Not EC90/curves (PDXE has no dose-response); the analog is tumour-volume response.
Panel A ranks CRC treatments by median response (negative = shrinkage); Panel B
shows which treatments are more effective in CRC than in the other tumour types.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.transforms as mtransforms  # noqa: E402
import numpy as np  # noqa: E402

RRB_MAROON = "#6E1426"
_SHRINK = "#6E1426"
_GROW = "#B7BAC2"


def render_pdxe_report(ranking, selectivity, *, out_path, n_crc_models=None):
    out_path = Path(out_path); out_path.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(12, 10.5))
    fig.suptitle("WARHEAD - Novartis PDXE: in-vivo CRC drug response",
                 color=RRB_MAROON, fontsize=15, fontweight="bold", y=.98)
    fig.text(.5, .955, f"patient-derived xenografts ({n_crc_models or '?'} CRC models).  "
             "response = best avg % tumour-volume change (negative = shrinkage).  "
             "no dose-response / EC90; no HCC arm in PDXE.", ha="center", fontsize=9, color="#555")

    gs = fig.add_gridspec(2, 1, height_ratios=[1.25, 1.0], hspace=.34, left=.22, right=.95, top=.925, bottom=.07)

    # Panel A: response ranking (waterfall)
    axA = fig.add_subplot(gs[0])
    r = ranking.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(r))
    colors = [_SHRINK if v < 0 else _GROW for v in r["median_response"]]
    axA.barh(y, r["median_response"], color=colors, edgecolor="#222", linewidth=.4)
    axA.axvline(0, color="#222", lw=.8)
    axA.axvline(-30, color=RRB_MAROON, lw=.9, ls="--")  # ~PR threshold
    axA.text(-30, len(r) - .2, "PR (-30%)", color=RRB_MAROON, fontsize=7, ha="center", va="bottom")
    axA.set_yticks(y); axA.set_yticklabels(r["Treatment"], fontsize=8)
    axA.set_xlabel("median best-avg tumour response (%)  -  negative = shrinkage")
    axA.set_title("A.  CRC single-agent treatments ranked by response "
                  "(maroon = shrinkage)", fontsize=10, loc="left")
    tr = mtransforms.blended_transform_factory(axA.transAxes, axA.transData)
    for yi, (_, row) in zip(y, r.iterrows()):
        axA.text(.992, yi, f"{str(row['target'])[:20]}  ORR {row['pct_objective_response']:.0f}%",
                 transform=tr, ha="right", va="center", fontsize=6.6, color="#555", family="monospace")

    # Panel B: CRC-vs-other selectivity
    axB = fig.add_subplot(gs[1])
    if len(selectivity):
        s = selectivity.head(12).iloc[::-1].reset_index(drop=True)
        yb = np.arange(len(s))
        axB.barh(yb, s["delta"], color=[RRB_MAROON if d > 0 else _GROW for d in s["delta"]],
                 edgecolor="#222", linewidth=.4)
        axB.axvline(0, color="#222", lw=.8)
        axB.set_yticks(yb); axB.set_yticklabels(s["Treatment"], fontsize=8)
        axB.set_xlabel("CRC selectivity  =  median response(other types) - median response(CRC)   (>0 = CRC more responsive)")
        axB.set_title("B.  Response selectivity: CRC vs other tumour types", fontsize=10, loc="left")
        trb = mtransforms.blended_transform_factory(axB.transAxes, axB.transData)
        for yi, (_, row) in zip(yb, s.iterrows()):
            axB.text(.992, yi, str(row["target"])[:22], transform=trb, ha="right", va="center",
                     fontsize=6.6, color="#555", family="monospace")
    fig.savefig(out_path, format=out_path.suffix.lstrip(".").lower() or "pdf", dpi=150)
    plt.close(fig)
    return out_path
