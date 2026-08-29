"""reports/proliferation_independence_real.pdf - real-data G2b (GDSC2 x DepMap).

With ~286 compounds the per-compound bar is unreadable, so the headline view is
proliferation dependence aggregated by MECHANISM (pathway): which mechanisms lose
potency in slow-growing lines, and which stay flat - the HCC payload question.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from ..config import load_gates  # noqa: E402
from ..stats import balance_weights, weighted_linregress  # noqa: E402

RRB_MAROON = "#6E1426"
_FLAT = "#2C7FB8"
_MITO = "#B5651D"


def _scatter(ax, merged, comp, dt_col, letter, tag):
    sub = merged[merged["compound_id"] == comp]
    x = sub[dt_col].to_numpy(float)
    y = sub["sensitivity"].to_numpy(float)
    ax.scatter(x, y, s=12, color=RRB_MAROON, alpha=.55, edgecolor="none")
    if len(x) > 3:
        w = balance_weights(x)
        f = weighted_linregress(x, y, w)
        xs = np.linspace(x.min(), x.max(), 40)
        ax.plot(xs, f.intercept + f.slope * xs, color="#222", lw=1.7)
        ax.text(.04, .06, f"std slope = {f.std_slope:+.2f}", transform=ax.transAxes,
                fontsize=8, color="#222")
    ax.set_title(f"{letter}.  {comp}\n({tag})", fontsize=9, loc="left")
    ax.set_xlabel("1 / growth rate  (proxy doubling time)")
    ax.set_ylabel("log10(IC50 [uM])")


def render_real_g2b_report(stats, sensitivity, model_meta, *, out_path, config=None,
                           min_pathway_n=3):
    cfg = config or load_gates()
    pcfg = cfg["g2"]["proliferation"]
    dt_col = pcfg["doubling_time_col"]
    smax = pcfg["std_slope_max"]
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    merged = sensitivity.merge(model_meta[["ModelID", dt_col]], on="ModelID", how="inner")

    # pathway-level aggregation
    pw = (stats.dropna(subset=["pathway"]).groupby("pathway")
          .agg(median_slope=("std_slope", "median"), n=("std_slope", "size"))
          .reset_index())
    pw = pw[pw["n"] >= min_pathway_n].sort_values("median_slope")

    fig = plt.figure(figsize=(12, 9))
    fig.suptitle("WARHEAD - G2b Proliferation Independence  (REAL: GDSC2 x DepMap)",
                 color=RRB_MAROON, fontsize=15, fontweight="bold", x=.5, y=.975)
    fig.text(.5, .945, f"log10(IC50) vs DepMap growth-rate proxy across {merged['ModelID'].nunique()} "
             "lines  -  which mechanisms need active proliferation?", ha="center", fontsize=9.5, color="#555")

    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1.0], hspace=.42, wspace=.28,
                          left=.30, right=.965, top=.90, bottom=.09)

    axA = fig.add_subplot(gs[0, :])
    y = np.arange(len(pw))
    colors = [_FLAT if abs(m) < smax else _MITO for m in pw["median_slope"]]
    axA.barh(y, pw["median_slope"], color=colors, edgecolor="#222", linewidth=.4)
    axA.axvspan(-smax, smax, color=RRB_MAROON, alpha=.08)
    axA.axvline(0, color="#222", lw=.8)
    axA.set_yticks(y); axA.set_yticklabels(pw["pathway"], fontsize=8)
    for yi, (_, row) in zip(y, pw.iterrows()):
        axA.text(row["median_slope"] + (.004 if row["median_slope"] >= 0 else -.004), yi,
                 f"n={int(row['n'])}", va="center",
                 ha="left" if row["median_slope"] >= 0 else "right", fontsize=6.5, color="#666")
    axA.set_xlabel("median standardised slope of log10(IC50) on proxy doubling time")
    axA.set_title("A.  Proliferation dependence by MECHANISM  "
                  "(shaded = independent band; orange = proliferation-gated)", fontsize=9, loc="left")

    st = stats.sort_values("std_slope")
    steep = st.iloc[-1]["compound_id"]
    flat = st.iloc[st["std_slope"].abs().to_numpy().argmin()]["compound_id"]
    _scatter(fig.add_subplot(gs[1, 0]), merged, steep, dt_col, "B", "most proliferation-dependent")
    _scatter(fig.add_subplot(gs[1, 1]), merged, flat, dt_col, "C", "proliferation-independent")

    npass = int((stats["q"] > pcfg["fdr_alpha"]).sum())
    fig.text(.30, .02, f"compounds = {len(stats)}   |   proliferation-independent (pass G2b) = {npass}"
             "   |   proxy doubling time = 1 / DepMap screen-inferred growth rate",
             fontsize=8, color="#666")
    fig.savefig(out_path, format=out_path.suffix.lstrip(".").lower() or "pdf", dpi=150)
    plt.close(fig)
    return out_path
