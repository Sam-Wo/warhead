"""reports/nci60_crc_selectivity.pdf - NCI-60 CRC (colon)-selectivity.

NCI-60 constraints are stated on the figure: z-scored (relative) activity so no
absolute EC90/curve; 7 colon lines; NO liver (no HCC). The value is breadth - the
selectivity is computed over the >25k-compound DTP library.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.transforms as mtransforms  # noqa: E402
import numpy as np  # noqa: E402

RRB_MAROON = "#6E1426"


def render_nci60_report(selectivity, annotated, *, out_path, n_total, n_colon):
    out_path = Path(out_path); out_path.parent.mkdir(parents=True, exist_ok=True)
    top = annotated.sort_values("delta_z", ascending=False).head(20).iloc[::-1].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(12, 0.42 * len(top) + 2.2))
    fig.suptitle("WARHEAD - NCI-60 CRC (colon)-selective compounds", color=RRB_MAROON,
                 fontsize=14, fontweight="bold", y=.99)
    fig.text(.5, .95, f"activity z-score, colon vs rest, over {n_total:,} DTP compounds ({n_colon} colon lines).  "
             "RELATIVE activity only (no absolute EC90/IC50/curves); NO liver line -> no HCC.  "
             "annotated (known MoA/clinical) compounds shown.", ha="center", fontsize=8.4, color="#555")
    y = np.arange(len(top))
    fda = top["FDA"].astype(str).str.contains("FDA|approved", case=False, na=False)
    clin = top["FDA"].astype(str).str.contains("Clinical", case=False, na=False)
    colors = [RRB_MAROON if a else ("#C06A7C" if c else "#9AA0A6") for a, c in zip(fda, clin)]
    ax.barh(y, top["delta_z"], color=colors, edgecolor="#222", linewidth=.4)
    ax.axvline(0, color="#222", lw=.8)
    ax.set_yticks(y); ax.set_yticklabels(top["drug"].astype(str).str.slice(0, 34), fontsize=8)
    ax.set_xlabel("CRC selectivity  =  mean activity z(colon) - z(rest)   (>0 = more active in colon)")
    tr = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
    for yi, (_, r) in zip(y, top.iterrows()):
        moa = str(r["MOA"])[:28] if str(r["MOA"]).strip() != "-" else ""
        fdas = str(r["FDA"])[:16] if str(r["FDA"]).strip() != "-" else ""
        ax.text(.992, yi, (moa + ("  " + fdas if fdas else "")), transform=tr, ha="right",
                va="center", fontsize=6.5, color="#555", family="monospace")
    fig.tight_layout(rect=[0, 0, 1, .93])
    fig.savefig(out_path, format=out_path.suffix.lstrip(".").lower() or "pdf", dpi=150)
    plt.close(fig)
    return out_path
