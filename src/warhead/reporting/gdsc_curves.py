"""reports/gdsc_top15_curves.pdf - measured dose-response for the 15 most potent
GDSC2 compounds, showing where IC50 and (extrapolated) EC90 fall vs the doses
actually tested. This is the visual answer to "why is EC90 >> IC50": these
compounds were screened over low, narrow windows that barely reach IC50, so 90%
kill is never observed and EC90 is a model extrapolation.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

RRB_MAROON = "#6E1426"
_MEAS = "#6E1426"
_FIT = "#2C7FB8"
_EXTRAP = "#EDE7E9"


def _gdsc_curve(conc_uM, ic50_uM, scal):
    # y = 1 / (1 + (c/IC50)^(1/scal))  (GDSC 2-param sigmoid, bottom = 0)
    return 1.0 / (1.0 + np.power(conc_uM / ic50_uM, 1.0 / scal))


def render_top15_fitted_curves(summary: pd.DataFrame, *, out_path: str | Path) -> Path:
    """Plot each top compound's GDSC FITTED dose-response (IC50 + recovered slope,
    bottom=0) with the tested concentration window shaded, so it is clear that the
    doses tested barely reach IC50 and that EC90 is an extrapolation.

    `summary`: drug_name, target, median_ic50_uM, median_scal, median_ec90_uM,
    min_conc_uM, max_conc_uM.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    summary = summary.reset_index(drop=True)
    n = len(summary)
    ncol, nrow = 3, int(np.ceil(n / 3))

    fig, axes = plt.subplots(nrow, ncol, figsize=(14, 3.05 * nrow))
    axes = np.atleast_1d(axes).ravel()
    fig.suptitle("WARHEAD - GDSC2 top-15 most potent compounds: fitted dose-response",
                 color=RRB_MAROON, fontsize=15, fontweight="bold", y=.997)
    fig.text(.5, .972,
             "GDSC fitted 2-param sigmoid (bottom=0).  green band = concentrations actually tested;  "
             "grey = beyond tested max (EC90 extrapolated).  dashed = IC50 & EC90;  dotted = 50% / 90% effect.",
             ha="center", fontsize=8.6, color="#555")

    for i, (_, s) in enumerate(summary.iterrows()):
        ax = axes[i]
        ic50_nM = s["median_ic50_uM"] * 1e3
        ec90_nM = s["median_ec90_uM"] * 1e3
        mint_nM = s["min_conc_uM"] * 1e3
        maxt_nM = s["max_conc_uM"] * 1e3
        v_at_max = float(_gdsc_curve(s["max_conc_uM"], s["median_ic50_uM"], s["median_scal"]))

        xlo = min(mint_nM, ic50_nM) / 3
        xhi = max(ec90_nM, maxt_nM) * 2.5
        ax.axvspan(mint_nM, maxt_nM, color="#2E7D6B", alpha=.10, zorder=0)   # tested window
        ax.axvspan(maxt_nM, xhi, color=_EXTRAP, zorder=0)                     # extrapolation

        xx = np.logspace(np.log10(xlo), np.log10(xhi), 240)
        ax.plot(xx, _gdsc_curve(xx / 1e3, s["median_ic50_uM"], s["median_scal"]),
                color=_MEAS, lw=1.8, zorder=3)
        ax.axhline(.5, color="#bbb", lw=.7, ls=":")
        ax.axhline(.1, color="#bbb", lw=.7, ls=":")
        ax.axvline(ic50_nM, color="#222", lw=1.0, ls="--", zorder=4)
        ax.axvline(ec90_nM, color=RRB_MAROON, lw=1.0, ls="--", zorder=4)
        ax.text(ic50_nM, 1.08, "IC50", color="#222", fontsize=6.5, ha="center")
        ax.text(ec90_nM, 1.08, "EC90", color=RRB_MAROON, fontsize=6.5, ha="center")

        ax.set_xscale("log")
        ax.set_xlim(xlo, xhi); ax.set_ylim(-0.05, 1.2)
        ax.set_title(f"{s['drug_name']}  ({s['target'] if isinstance(s['target'], str) else 'n/a'})",
                     fontsize=8.5, loc="left")
        ax.tick_params(labelsize=7)
        ax.text(.03, .06,
                f"IC50 {ic50_nM:.0f} nM · EC90 {ec90_nM:.0f} nM\n"
                f"tested {mint_nM:.2g}-{maxt_nM:.0f} nM · viab@max {v_at_max:.2f}",
                transform=ax.transAxes, fontsize=6.5, color="#444", va="bottom")
        if i % ncol == 0:
            ax.set_ylabel("viability (fitted)", fontsize=8)
        if i >= n - ncol:
            ax.set_xlabel("concentration (nM, log)", fontsize=8)

    for j in range(n, len(axes)):
        axes[j].set_visible(False)

    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    fig.legend(handles=[
        Line2D([0], [0], color=_MEAS, lw=2, label="GDSC fitted curve"),
        Patch(facecolor="#2E7D6B", alpha=.10, label="tested concentration window"),
        Patch(facecolor=_EXTRAP, label="beyond tested max (EC90 extrapolated)"),
    ], loc="lower center", ncol=3, fontsize=8, frameon=False, bbox_to_anchor=(.5, .0005))
    fig.tight_layout(rect=[0, .02, 1, .958])
    fig.savefig(out_path, format=out_path.suffix.lstrip(".").lower() or "pdf", dpi=150)
    plt.close(fig)
    return out_path


def render_top15_curves(pooled: pd.DataFrame, summary: pd.DataFrame, *, out_path: str | Path) -> Path:
    """`pooled`: drug_id, conc_uM, median, q1, q3 (from analysis.gdsc_curves).
    `summary`: one row per drug_id with drug_name, target, median_ic50_uM,
    median_scal, median_ec90_uM, max_conc_uM."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    summary = summary.reset_index(drop=True)
    n = len(summary)
    ncol, nrow = 3, int(np.ceil(n / 3))

    fig, axes = plt.subplots(nrow, ncol, figsize=(14, 3.05 * nrow))
    axes = np.atleast_1d(axes).ravel()
    fig.suptitle("WARHEAD - GDSC2 top-15 most potent compounds: measured dose-response",
                 color=RRB_MAROON, fontsize=15, fontweight="bold", y=.997)
    fig.text(.5, .983,
             "points = measured median viability across lines (IQR band).  grey = beyond the "
             "tested max dose (EC90 extrapolated).  dashed = IC50 (fitted) and EC90.",
             ha="center", fontsize=9, color="#555")

    for i, (_, s) in enumerate(summary.iterrows()):
        ax = axes[i]
        cur = pooled[pooled["drug_id"] == s["drug_id"]].sort_values("conc_uM")
        if not len(cur):
            ax.set_visible(False); continue
        c_nM = cur["conc_uM"].to_numpy() * 1e3
        ic50_nM = s["median_ic50_uM"] * 1e3
        ec90_nM = s["median_ec90_uM"] * 1e3
        maxt_nM = s["max_conc_uM"] * 1e3
        v_at_max = float(cur.sort_values("conc_uM")["median"].iloc[-1])

        xlo, xhi = c_nM.min() / 2, max(ec90_nM, maxt_nM) * 2.2
        ax.axvspan(maxt_nM, xhi, color=_EXTRAP, zorder=0)  # extrapolation zone

        # measured curve + IQR
        ax.fill_between(c_nM, cur["q1"], cur["q3"], color=_MEAS, alpha=.16, zorder=1)
        ax.plot(c_nM, cur["median"], "-o", color=_MEAS, ms=3.5, lw=1.4, zorder=3)
        # GDSC fitted model across full x
        xx = np.logspace(np.log10(xlo), np.log10(xhi), 200)
        ax.plot(xx, _gdsc_curve(xx / 1e3, s["median_ic50_uM"], s["median_scal"]),
                color=_FIT, lw=1.4, ls="-", alpha=.9, zorder=2)

        ax.axhline(.5, color="#bbb", lw=.7, ls=":")
        ax.axhline(.1, color="#bbb", lw=.7, ls=":")
        for xv, lab, col in [(ic50_nM, "IC50", "#222"), (ec90_nM, "EC90", RRB_MAROON)]:
            ax.axvline(xv, color=col, lw=1.0, ls="--", zorder=4)
        ax.set_xscale("log")
        ax.set_xlim(xlo, xhi); ax.set_ylim(-0.05, 1.15)
        ax.set_title(f"{s['drug_name']}  ({s['target'] if isinstance(s['target'], str) else 'n/a'})",
                     fontsize=8.5, loc="left")
        ax.tick_params(labelsize=7)
        ax.text(.03, .06,
                f"IC50 {ic50_nM:.0f} nM · EC90 {ec90_nM:.0f} nM\n"
                f"max dose {maxt_nM:.0f} nM · viab@max {v_at_max:.2f}",
                transform=ax.transAxes, fontsize=6.6, color="#444", va="bottom")
        if i % ncol == 0:
            ax.set_ylabel("viability", fontsize=8)
        if i >= n - ncol:
            ax.set_xlabel("concentration (nM, log)", fontsize=8)

    for j in range(n, len(axes)):
        axes[j].set_visible(False)

    # legend proxies
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    fig.legend(handles=[
        Line2D([0], [0], color=_MEAS, marker="o", ms=4, label="measured median (IQR band)"),
        Line2D([0], [0], color=_FIT, label="GDSC fitted model"),
        Patch(facecolor=_EXTRAP, label="beyond tested max dose"),
    ], loc="lower center", ncol=3, fontsize=8, frameon=False, bbox_to_anchor=(.5, .0005))

    fig.tight_layout(rect=[0, .02, 1, .975])
    fig.savefig(out_path, format=out_path.suffix.lstrip(".").lower() or "pdf", dpi=150)
    plt.close(fig)
    return out_path
