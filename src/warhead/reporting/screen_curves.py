"""Top-N dose-response curves for any screen.

- render_fitted_curves: draws the fitted 4PL with a FREE lower asymptote
  viability(d) = lower + (upper-lower)/(1+(d/ec50)^slope), i.e. it shows the real
  Emax plateau (used for PRISM / CTRP fits, unlike GDSC's bottom=0).
- render_measured_curves: draws measured median viability + IQR (used for CTRP
  raw wells), with an optional fitted overlay.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

RRB_MAROON = "#6E1426"
_FIT = "#6E1426"
_EXTRAP = "#EDE7E9"
_TESTED = "#2E7D6B"


def _grid(n):
    ncol = 3
    return ncol, int(np.ceil(n / ncol))


def _fourpl_free(d_uM, upper, lower, ec50, slope):
    return lower + (upper - lower) / (1.0 + np.power(d_uM / ec50, abs(slope)))


def render_fitted_curves(summary: pd.DataFrame, *, source, out_path,
                         min_conc_uM=6e-4, max_conc_uM=10.0):
    """summary rows: compound, target, ic50_uM, ec50_uM, slope, upper, lower,
    ec90_uM (+ optional per-row min_conc_uM/max_conc_uM)."""
    out_path = Path(out_path); out_path.parent.mkdir(parents=True, exist_ok=True)
    s = summary.reset_index(drop=True); n = len(s)
    ncol, nrow = _grid(n)
    fig, axes = plt.subplots(nrow, ncol, figsize=(14, 3.05 * nrow))
    axes = np.atleast_1d(axes).ravel()
    fig.suptitle(f"WARHEAD - {source} top-{n} most potent: fitted dose-response (with Emax)",
                 color=RRB_MAROON, fontsize=15, fontweight="bold", y=.997)
    fig.text(.5, .972, "4PL with a free lower asymptote (Emax).  green = tested window;  grey = beyond "
             "tested max (EC90 extrapolated).  dashed = IC50 & EC90;  dotted line = Emax plateau.",
             ha="center", fontsize=8.6, color="#555")

    for i, r in s.iterrows():
        ax = axes[i]
        ic50_nM = r["ic50_uM"] * 1e3; ec90_nM = r["ec90_uM"] * 1e3
        mint = r.get("min_conc_uM", min_conc_uM) * 1e3
        maxt = r.get("max_conc_uM", max_conc_uM) * 1e3
        emax = r["lower"]
        xlo = min(mint, ic50_nM if np.isfinite(ic50_nM) else mint) / 3
        xhi = max(ec90_nM, maxt) * 2.5
        ax.axvspan(mint, maxt, color=_TESTED, alpha=.10, zorder=0)
        ax.axvspan(maxt, xhi, color=_EXTRAP, zorder=0)
        xx = np.logspace(np.log10(xlo), np.log10(xhi), 240)
        ax.plot(xx, _fourpl_free(xx / 1e3, r["upper"], r["lower"], r["ec50_uM"], r["slope"]),
                color=_FIT, lw=1.8, zorder=3)
        ax.axhline(.5, color="#ccc", lw=.7, ls=":")
        ax.axhline(emax, color=RRB_MAROON, lw=.8, ls=":", alpha=.6)
        if np.isfinite(ic50_nM):
            ax.axvline(ic50_nM, color="#222", lw=1.0, ls="--", zorder=4)
            ax.text(ic50_nM, 1.08, "IC50", color="#222", fontsize=6.5, ha="center")
        ax.axvline(ec90_nM, color=RRB_MAROON, lw=1.0, ls="--", zorder=4)
        ax.text(ec90_nM, 1.08, "EC90", color=RRB_MAROON, fontsize=6.5, ha="center")
        ax.set_xscale("log"); ax.set_xlim(xlo, xhi); ax.set_ylim(-.05, 1.2)
        ax.set_title(f"{r['compound']}  ({str(r['target'])[:22] if pd.notna(r['target']) else 'n/a'})",
                     fontsize=8.3, loc="left")
        ax.tick_params(labelsize=7)
        ax.text(.03, .06, f"IC50 {ic50_nM:.0f} nM · EC90 {ec90_nM:.0f} nM\nEmax {emax:.2f}",
                transform=ax.transAxes, fontsize=6.5, color="#444", va="bottom")
        if i % ncol == 0:
            ax.set_ylabel("viability (fitted)", fontsize=8)
        if i >= n - ncol:
            ax.set_xlabel("concentration (nM, log)", fontsize=8)
    for j in range(n, len(axes)):
        axes[j].set_visible(False)
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    fig.legend(handles=[Line2D([0], [0], color=_FIT, lw=2, label="fitted 4PL (with Emax)"),
                        Patch(facecolor=_TESTED, alpha=.10, label="tested window"),
                        Patch(facecolor=_EXTRAP, label="beyond tested max")],
               loc="lower center", ncol=3, fontsize=8, frameon=False, bbox_to_anchor=(.5, .0005))
    fig.tight_layout(rect=[0, .02, 1, .958])
    fig.savefig(out_path, format=out_path.suffix.lstrip(".").lower() or "pdf", dpi=150)
    plt.close(fig)
    return out_path


def pool_measured_curves(raw: pd.DataFrame, *, nbins=14, min_bin_n=8, range_frac=0.5):
    """Pool raw per-well CTRP viability into a clean median+IQR curve per compound.

    CTRP doses each line on a compound-specific grid, and different lines get
    different grids (jittered, variable density, sometimes a different range).
    Plotting the raw pooled points therefore zig-zags. This collapses each
    compound onto a common set of log-spaced concentration bins:

      raw: compound, conc_uM, viability  (viability already a FRACTION, not %),
           optional model_id to identify a "line".

    1. drop lines whose max tested dose < range_frac x the compound's median
       per-line max (removes shallow/off-range lines that cause the zig-zag);
    2. cut each compound's log10(conc) into `nbins` equal-width bins, take the
       median viability + IQR per bin;
    3. drop bins with fewer than `min_bin_n` wells.

    Returns pooled: compound, conc_uM, median, q1, q3, n.
    """
    d = raw.copy()
    d = d[np.isfinite(d["conc_uM"]) & (d["conc_uM"] > 0) & np.isfinite(d["viability"])]
    line_col = "model_id" if "model_id" in d.columns else ("cellid" if "cellid" in d.columns else None)
    out = []
    for cmp, g in d.groupby("compound"):
        if line_col is not None:
            line_max = g.groupby(line_col)["conc_uM"].max()
            keep = line_max[line_max >= range_frac * line_max.median()].index
            g = g[g[line_col].isin(keep)]
        if len(g) < min_bin_n:
            continue
        lg = np.log10(g["conc_uM"].to_numpy())
        edges = np.linspace(lg.min(), lg.max(), nbins + 1)
        idx = np.clip(np.digitize(lg, edges[1:-1]), 0, nbins - 1)
        centres = 10 ** ((edges[:-1] + edges[1:]) / 2)
        v = g["viability"].to_numpy()
        for b in range(nbins):
            m = idx == b
            if m.sum() < min_bin_n:
                continue
            vb = v[m]
            out.append({"compound": cmp, "conc_uM": float(centres[b]),
                        "median": float(np.median(vb)), "q1": float(np.percentile(vb, 25)),
                        "q3": float(np.percentile(vb, 75)), "n": int(m.sum())})
    return pd.DataFrame(out)


def render_measured_curves(pooled: pd.DataFrame, summary: pd.DataFrame, *, source, out_path):
    """pooled: compound, conc_uM, median, q1, q3.  summary: compound, target,
    ic50_uM, ec90_uM, max_conc_uM (for the markers)."""
    out_path = Path(out_path); out_path.parent.mkdir(parents=True, exist_ok=True)
    s = summary.reset_index(drop=True); n = len(s)
    ncol, nrow = _grid(n)
    fig, axes = plt.subplots(nrow, ncol, figsize=(14, 3.05 * nrow))
    axes = np.atleast_1d(axes).ravel()
    fig.suptitle(f"WARHEAD - {source} top-{n} most potent: measured dose-response",
                 color=RRB_MAROON, fontsize=15, fontweight="bold", y=.997)
    fig.text(.5, .972, "points = measured median viability across lines (IQR band).  "
             "dashed = IC50 & EC90.", ha="center", fontsize=8.6, color="#555")
    for i, r in s.iterrows():
        ax = axes[i]
        cur = pooled[pooled["compound"] == r["compound"]].sort_values("conc_uM")
        if not len(cur):
            ax.set_visible(False); continue
        c_nM = cur["conc_uM"].to_numpy() * 1e3
        ax.fill_between(c_nM, cur["q1"], cur["q3"], color=_FIT, alpha=.15, zorder=1)
        ax.plot(c_nM, cur["median"], "-o", color=_FIT, ms=3.2, lw=1.3, zorder=3)
        ax.axhline(.5, color="#ccc", lw=.7, ls=":"); ax.axhline(.1, color="#ccc", lw=.7, ls=":")
        for xv, col in [(r["ic50_uM"] * 1e3, "#222"), (r["ec90_uM"] * 1e3, RRB_MAROON)]:
            if np.isfinite(xv):
                ax.axvline(xv, color=col, lw=1.0, ls="--", zorder=4)
        ax.set_xscale("log"); ax.set_ylim(-.05, 1.2)
        ax.set_title(f"{r['compound']}  ({str(r['target'])[:22] if pd.notna(r['target']) else 'n/a'})",
                     fontsize=8.3, loc="left")
        ax.tick_params(labelsize=7)
        ax.text(.03, .06, f"IC50 {r['ic50_uM']*1e3:.0f} nM · EC90 {r['ec90_uM']*1e3:.0f} nM",
                transform=ax.transAxes, fontsize=6.5, color="#444", va="bottom")
        if i % ncol == 0:
            ax.set_ylabel("viability", fontsize=8)
        if i >= n - ncol:
            ax.set_xlabel("concentration (nM, log)", fontsize=8)
    for j in range(n, len(axes)):
        axes[j].set_visible(False)
    fig.tight_layout(rect=[0, .01, 1, .958])
    fig.savefig(out_path, format=out_path.suffix.lstrip(".").lower() or "pdf", dpi=150)
    plt.close(fig)
    return out_path
