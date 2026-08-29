"""reports/gdsc_ec90_selectivity.pdf - GDSC EC90 potency and HCC/CRC selectivity.

Rows 1-2: which compounds reach 90% effect at the lowest concentration in each
indication (top 20, nM), with target annotated and extrapolated EC90 hatched.
Row 3: potency vs tissue-selectivity - compounds in the upper right are both potent
and selectively more active in the indication (the mechanistic-advantage quadrant).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.transforms as mtransforms  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter  # noqa: E402

RRB_MAROON = "#6E1426"
_HILITE = "#6E1426"
_GREY = "#B7BAC2"
_SEL_WEAK = "#C98A2E"


def _short(target: str, n: int = 34) -> str:
    if not isinstance(target, str) or not target or target.lower() == "nan":
        return "(target n/a)"
    return target if len(target) <= n else target[: n - 1] + "…"


def _rank_panel(ax, rank: pd.DataFrame, title: str, n: int = 20):
    r = rank.head(n).iloc[::-1].reset_index(drop=True)
    y = np.arange(len(r))
    for yi, (_, row) in zip(y, r.iterrows()):
        extra = bool(row["median_extrapolated"])
        ax.barh(yi, row["median_ec90_nM"], color=RRB_MAROON, edgecolor="#3a0a13",
                linewidth=.5, height=.74, hatch="////" if extra else None)
    ax.set_yticks(y)
    ax.set_yticklabels(r["drug_name"], fontsize=8.5)
    ax.set_xscale("log")
    ax.set_ylim(-0.7, len(r) - 0.3)

    # headroom so the target column (right) never collides with the bars
    vmax = float(r["median_ec90_nM"].max())
    vmin = float(r["median_ec90_nM"].min())
    ax.set_xlim(vmin / 2.5, vmax * 30)
    ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5), numticks=12))
    ax.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * .1, numticks=40))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.tick_params(axis="x", labelsize=8)
    ax.set_xlabel("median EC90 in indication (nM, log scale)", fontsize=9)
    ax.set_title(title, fontsize=10, loc="left")

    # right-aligned target column (axes-x, data-y) so targets read as a list
    tr = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
    for yi, (_, row) in zip(y, r.iterrows()):
        ax.text(0.992, yi, _short(row["target"]), transform=tr, ha="right", va="center",
                fontsize=6.8, color="#555", family="monospace")


def _sel_panel(ax, sel: pd.DataFrame, title: str, annotate: int = 6):
    if not len(sel):
        ax.text(.5, .5, "no compounds", ha="center", va="center", transform=ax.transAxes,
                color="#888"); ax.set_title(title, fontsize=10, loc="left"); return
    pot = sel["potency_in"].to_numpy()
    dp = sel["delta_potency"].to_numpy()
    is_sp = sel["selective_potent"].to_numpy()
    is_s = sel["selective"].to_numpy() & ~is_sp
    colors = np.where(is_sp, _HILITE, np.where(is_s, _SEL_WEAK, _GREY))
    sizes = np.where(is_sp, 48, np.where(is_s, 30, 15))
    ax.axhline(0, color="#888", lw=.8, ls=":")
    ax.scatter(pot, dp, s=sizes, c=colors, edgecolor="#222", linewidth=.4, zorder=3)
    ax.margins(x=.12, y=.12)
    xlo, xhi = ax.get_xlim()
    xmid = 0.5 * (xlo + xhi)
    ann = sel[sel["selective_potent"]].head(annotate)
    if not len(ann):
        ann = sel.head(annotate)
    for _, row in ann.iterrows():
        x, yv = row["potency_in"], row["delta_potency"]
        right = x > xmid
        col = RRB_MAROON if row.get("selective_potent") else "#777"
        ax.annotate(f"{row['drug_name']} [{_short(row['target'], 12)}]", (x, yv),
                    fontsize=6.6, color=col,
                    xytext=(-4 if right else 4, 3), textcoords="offset points",
                    ha="right" if right else "left")
    ax.set_xlabel("potency in indication  =  -log10(median IC50 / nM)", fontsize=9)
    ax.set_ylabel("selectivity  =  log-potency(in) - (rest)", fontsize=9)
    ax.set_title(title, fontsize=10, loc="left")


def render_gdsc_report(rankings: dict, selectivities: dict, *, out_path: str | Path,
                       counts: dict | None = None) -> Path:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    counts = counts or {}
    tc = counts.get("total_compounds", "?")
    tl = counts.get("total_lines", "?")

    fig = plt.figure(figsize=(13, 16))
    fig.suptitle("WARHEAD - GDSC2 EC90 Potency & HCC / CRC Selectivity",
                 color=RRB_MAROON, fontsize=16, fontweight="bold", x=.5, y=.985)
    fig.text(.5, .966,
             f"GDSC2 dataset: {tc} compounds x {tl} cell lines.  "
             "EC90 = concentration for 90% of fitted max effect (2-param fit, bottom=0);  "
             "hatched bar = median EC90 extrapolated beyond the tested max dose.",
             ha="center", fontsize=9, color="#555")

    gs = fig.add_gridspec(3, 2, height_ratios=[1.35, 1.35, 1.0], hspace=.30, wspace=.26,
                          left=.15, right=.965, top=.94, bottom=.055)

    crc_n = counts.get("CRC_lines", "?")
    hcc_n = counts.get("HCC_lines", "?")
    axA = fig.add_subplot(gs[0, :])
    _rank_panel(axA, rankings["CRC"], f"A.  CRC (COREAD) - top 20 lowest EC90   ({crc_n} cell lines)")
    axA.legend(handles=[Patch(facecolor=RRB_MAROON, hatch="////", edgecolor="#3a0a13",
                              label="EC90 extrapolated beyond tested max")],
               loc="lower right", fontsize=7.5, framealpha=.9)
    _rank_panel(fig.add_subplot(gs[1, :]), rankings["HCC"],
                f"B.  HCC (LIHC) - top 20 lowest EC90   ({hcc_n} cell lines)")
    _sel_panel(fig.add_subplot(gs[2, 0]), selectivities["CRC"],
               "C.  CRC selectivity (maroon = selective & potent; amber = selective, weak)")
    _sel_panel(fig.add_subplot(gs[2, 1]), selectivities["HCC"], "D.  HCC selectivity")

    n_crc = int(selectivities["CRC"]["selective_potent"].sum()) if len(selectivities["CRC"]) else 0
    n_hcc = int(selectivities["HCC"]["selective_potent"].sum()) if len(selectivities["HCC"]) else 0
    fig.text(.15, .012,
             f"selective & potent (q<0.1, IC50<=1000 nM):  CRC = {n_crc}   HCC = {n_hcc}"
             "     |  selectivity anchored on fitted IC50; EC90 drives the ranking.  "
             "hypothesis-generating: tissue selectivity can reflect a co-enriched dependency.",
             fontsize=7.5, color="#666")

    fig.savefig(out_path, format=out_path.suffix.lstrip(".").lower() or "pdf", dpi=150)
    plt.close(fig)
    return out_path
