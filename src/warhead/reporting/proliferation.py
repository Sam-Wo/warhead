"""reports/proliferation_independence.pdf - the standalone G2b deliverable for
the HCC / MASH ADC argument (WARHEAD.md sec 7)."""
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


def _fit_line(x, y):
    w = balance_weights(x)
    f = weighted_linregress(x, y, w)
    xs = np.linspace(np.min(x), np.max(x), 50)
    return xs, f.intercept + f.slope * xs, f


def render_proliferation_report(
    stats: pd.DataFrame,
    sensitivity: pd.DataFrame,
    model_meta: pd.DataFrame,
    *,
    out_path: str | Path,
    compound_col: str = "compound_id",
    config: dict | None = None,
) -> Path:
    cfg = config or load_gates()
    pcfg = cfg["g2"]["proliferation"]
    dt_col = pcfg["doubling_time_col"]
    smax = pcfg["std_slope_max"]
    alpha = pcfg["fdr_alpha"]
    min_lines = pcfg["min_lines"]
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    merged = sensitivity.merge(model_meta[["ModelID", dt_col]], on="ModelID", how="inner")
    st = stats.sort_values("std_slope").reset_index(drop=True)
    # Pass = slope not significantly different from zero (the G2b criterion).
    st["passes_g2b"] = (st["q"] > alpha) & (st["n_lines"] >= min_lines)

    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("WARHEAD - G2b Proliferation Independence", color=RRB_MAROON,
                 fontsize=16, fontweight="bold", x=0.5, y=0.97)
    fig.text(0.5, 0.935, "Potency loss (log10 IC50) vs DepMap doubling time  -  keep the flat ones",
             ha="center", fontsize=10, color="#444")

    gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 1.0], hspace=0.42, wspace=0.28,
                          left=0.16, right=0.96, top=0.90, bottom=0.09)

    # Panel A: standardised slope per compound.
    axA = fig.add_subplot(gs[0, :])
    colors = [_FLAT if p else _MITO for p in st["passes_g2b"]]
    ypos = np.arange(len(st))
    axA.barh(ypos, st["std_slope"], color=colors, edgecolor="#222", linewidth=0.4)
    axA.axvspan(-smax, smax, color=RRB_MAROON, alpha=0.08)
    axA.axvline(0, color="#222", lw=0.8)
    axA.axvline(smax, color=RRB_MAROON, ls="--", lw=0.9)
    axA.axvline(-smax, color=RRB_MAROON, ls="--", lw=0.9)
    axA.set_yticks(ypos)
    axA.set_yticklabels(st[compound_col], fontsize=8)
    axA.set_xlabel("standardised slope of log10(IC50) on doubling time")
    axA.set_title("A.  Proliferation dependence per compound "
                  "(shaded band = independent; blue = passes G2b)", fontsize=9, loc="left")

    # Panels B/C: example scatters - steepest (mitotic) and flattest (independent).
    steep = st.iloc[-1][compound_col]
    flat = st.iloc[(st["std_slope"].abs()).idxmin()][compound_col] if len(st) else None
    for ax_i, (comp, tag) in enumerate([(steep, "mitotic-dependent"), (flat, "proliferation-independent")]):
        ax = fig.add_subplot(gs[1, ax_i])
        sub = merged[merged[compound_col] == comp]
        x = sub[dt_col].to_numpy(float)
        y = sub["sensitivity"].to_numpy(float)
        ax.scatter(x, y, s=18, color=RRB_MAROON, alpha=0.7, edgecolor="none")
        if len(x) > 3:
            xs, ys, f = _fit_line(x, y)
            ax.plot(xs, ys, color="#222", lw=1.6)
            ax.text(0.04, 0.06, f"std slope = {f.std_slope:+.2f}", transform=ax.transAxes,
                    fontsize=8, color="#222")
        ax.set_title(f"{'B' if ax_i == 0 else 'C'}.  {comp}\n({tag})", fontsize=9, loc="left")
        ax.set_xlabel("doubling time (h)")
        ax.set_ylabel("log10(IC50 [M])")

    n_pass = int(st["passes_g2b"].sum())
    fig.text(0.16, 0.02,
             f"n compounds = {len(st)}   |   pass G2b (proliferation-independent) = {n_pass}"
             f"   |   band = |std slope| < {smax}",
             fontsize=8, color="#555")

    fmt = out_path.suffix.lstrip(".").lower() or "pdf"
    fig.savefig(out_path, format=fmt, dpi=150)
    plt.close(fig)
    return out_path
