"""Generic per-source report: EC90 potency ranking + HCC/CRC selectivity.

Works from the canonical rank_potency / selectivity frames (analysis.screen_potency),
so GDSC / PRISM / CTRP all render identically. Ranking bars are coloured by
clinical phase (validation), hatched where EC90 is extrapolated, and annotated
with each compound's target and Emax (residual viability).
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
_GREY = "#B7BAC2"
_SEL_WEAK = "#C98A2E"

# clinical-phase palette (validation)
_PHASE_COLOR = {"Launched": "#6E1426", "Phase 3": "#9E3A50", "Phase 2/Phase 3": "#9E3A50",
                "Phase 2": "#C06A7C", "Phase 1/Phase 2": "#C06A7C", "Phase 1": "#D9A6B0",
                "Preclinical": "#9AA0A6", "Withdrawn": "#5A5A5A"}


def _phase_color(p):
    return _PHASE_COLOR.get(str(p), "#C7B8BC")


def _short(t, n=30):
    if not isinstance(t, str) or not t or t.lower() == "none" or t.lower() == "nan":
        return "target n/a"
    return t if len(t) <= n else t[: n - 1] + "…"


def _rank_panel(ax, rank, title, n=20):
    r = rank.head(n).iloc[::-1].reset_index(drop=True)
    y = np.arange(len(r))
    for yi, (_, row) in zip(y, r.iterrows()):
        extra = (row.get("frac_ec90_extrapolated") or 0) > 0.5
        ax.barh(yi, row["median_ec90_nM"], color=_phase_color(row.get("clinical_phase")),
                edgecolor="#3a0a13", linewidth=.5, height=.74, hatch="////" if extra else None)
    ax.set_yticks(y); ax.set_yticklabels(r["compound"], fontsize=8.5)
    ax.set_xscale("log"); ax.set_ylim(-.7, len(r) - .3)
    vmax = float(r["median_ec90_nM"].max()); vmin = float(r["median_ec90_nM"].min())
    ax.set_xlim(vmin / 2.5, vmax * 60)
    ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5), numticks=14))
    ax.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * .1, numticks=40))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.tick_params(axis="x", labelsize=8)
    ax.set_xlabel("median EC90 in indication (nM, log scale)", fontsize=9)
    ax.set_title(title, fontsize=10, loc="left")
    tr = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
    for yi, (_, row) in zip(y, r.iterrows()):
        em = row.get("median_emax")
        emx = f"  Emax {em:.2f}" if pd.notna(em) else ""
        ax.text(.992, yi, _short(row["target"]) + emx, transform=tr, ha="right", va="center",
                fontsize=6.6, color="#555", family="monospace")


def _sel_panel(ax, sel, title, annotate=6):
    if not len(sel):
        ax.text(.5, .5, "no compounds", ha="center", va="center", transform=ax.transAxes, color="#888")
        ax.set_title(title, fontsize=10, loc="left"); return
    is_sp = sel["selective_potent"].to_numpy(); is_s = sel["selective"].to_numpy() & ~is_sp
    colors = np.where(is_sp, RRB_MAROON, np.where(is_s, _SEL_WEAK, _GREY))
    sizes = np.where(is_sp, 48, np.where(is_s, 30, 15))
    ax.axhline(0, color="#888", lw=.8, ls=":")
    ax.scatter(sel["potency_in"], sel["delta_potency"], s=sizes, c=colors,
               edgecolor="#222", linewidth=.4, zorder=3)
    ax.margins(x=.12, y=.12)
    xlo, xhi = ax.get_xlim(); xmid = .5 * (xlo + xhi)
    ann = sel[sel["selective_potent"]].head(annotate)
    if not len(ann):
        ann = sel.head(3)  # nothing selective -> only flag the top few, avoid clutter
    for _, row in ann.iterrows():
        right = row["potency_in"] > xmid
        col = RRB_MAROON if row.get("selective_potent") else "#777"
        ax.annotate(f"{row['compound']} [{_short(row['target'], 12)}]",
                    (row["potency_in"], row["delta_potency"]), fontsize=6.6, color=col,
                    xytext=(-4 if right else 4, 3), textcoords="offset points",
                    ha="right" if right else "left")
    ax.set_xlabel("potency in indication  =  -log10(median IC50 / nM)", fontsize=9)
    ax.set_ylabel("selectivity  =  log-potency(in) - (rest)", fontsize=9)
    ax.set_title(title, fontsize=10, loc="left")


def render_screen_report(source, rankings, selectivities, *, out_path, counts=None, subtitle=""):
    out_path = Path(out_path); out_path.parent.mkdir(parents=True, exist_ok=True)
    counts = counts or {}
    fig = plt.figure(figsize=(13, 16.2))
    fig.suptitle(f"WARHEAD - {source} EC90 Potency & HCC / CRC Selectivity",
                 color=RRB_MAROON, fontsize=16, fontweight="bold", x=.5, y=.988)
    head = (f"{source}: {counts.get('total_compounds','?')} compounds x "
            f"{counts.get('total_lines','?')} lines.  bar colour = clinical phase; "
            f"hatched = EC90 extrapolated; Emax = residual viability (lower = deeper kill).")
    fig.text(.5, .967, subtitle or head, ha="center", fontsize=8.8, color="#555")

    gs = fig.add_gridspec(3, 2, height_ratios=[1.35, 1.35, 1.0], hspace=.30, wspace=.26,
                          left=.16, right=.965, top=.945, bottom=.06)
    _rank_panel(fig.add_subplot(gs[0, :]), rankings["CRC"],
                f"A.  CRC - top 20 lowest EC90   ({counts.get('CRC_lines','?')} lines)")
    _rank_panel(fig.add_subplot(gs[1, :]), rankings["HCC"],
                f"B.  HCC - top 20 lowest EC90   ({counts.get('HCC_lines','?')} lines)")
    _sel_panel(fig.add_subplot(gs[2, 0]), selectivities["CRC"],
               "C.  CRC selectivity (maroon = selective & potent)")
    _sel_panel(fig.add_subplot(gs[2, 1]), selectivities["HCC"], "D.  HCC selectivity")

    handles = [Patch(facecolor=_PHASE_COLOR["Launched"], label="Launched"),
               Patch(facecolor=_PHASE_COLOR["Phase 3"], label="Phase 2/3"),
               Patch(facecolor=_PHASE_COLOR["Phase 1"], label="Phase 1"),
               Patch(facecolor=_PHASE_COLOR["Preclinical"], label="Preclinical"),
               Patch(facecolor="white", edgecolor="#3a0a13", hatch="////", label="EC90 extrapolated")]
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=8, frameon=False,
               bbox_to_anchor=(.5, .012))
    fig.text(.16, .002, "selectivity anchored on fitted IC50; ranking on EC90 (Emax-filtered). "
             "hypothesis-generating: tissue selectivity can reflect a co-enriched dependency.",
             fontsize=7.5, color="#666")
    fig.savefig(out_path, format=out_path.suffix.lstrip(".").lower() or "pdf", dpi=150)
    plt.close(fig)
    return out_path
