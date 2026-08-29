"""reports/gdsc_ec90_selectivity.pdf - GDSC EC90 potency and HCC/CRC selectivity.

Top row: which compounds reach 90% effect at the lowest concentration in each
indication. Bottom row: potency vs tissue-selectivity - compounds in the upper
right are both potent and selectively more active in the indication (the
mechanistic-advantage quadrant).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

RRB_MAROON = "#6E1426"
_HILITE = "#6E1426"
_GREY = "#B7BAC2"
_SEL_WEAK = "#C98A2E"


def _rank_panel(ax, rank: pd.DataFrame, title: str, n: int = 12):
    r = rank.head(n).iloc[::-1]
    y = np.arange(len(r))
    ax.barh(y, r["median_ec90_uM"], color=RRB_MAROON, edgecolor="#222", linewidth=.4, height=.72)
    ax.set_yticks(y)
    ax.set_yticklabels(r["drug_name"], fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("median EC90 in indication (uM, log)")
    ax.set_title(title, fontsize=9.5, loc="left")
    # mark how much of the EC90 is observed vs extrapolated
    for yi, (_, row) in zip(y, r.iterrows()):
        ax.text(row["median_ec90_uM"] * 1.15, yi,
                f"{row['frac_ec90_within_range']*100:.0f}% obs",
                va="center", fontsize=6.5, color="#666")
    ax.margins(x=.18)


def _sel_panel(ax, sel: pd.DataFrame, title: str, annotate: int = 6):
    if not len(sel):
        ax.text(.5, .5, "no compounds", ha="center", va="center", transform=ax.transAxes,
                color="#888"); ax.set_title(title, fontsize=9.5, loc="left"); return
    pot = sel["potency_in"].to_numpy()
    dp = sel["delta_potency"].to_numpy()
    is_sp = sel["selective_potent"].to_numpy()
    is_s = sel["selective"].to_numpy() & ~is_sp
    colors = np.where(is_sp, _HILITE, np.where(is_s, _SEL_WEAK, _GREY))
    sizes = np.where(is_sp, 46, np.where(is_s, 30, 16))
    ax.axhline(0, color="#888", lw=.8, ls=":")
    ax.axvline(0, color="#ccc", lw=.7, ls=":")  # EC90 = 1 uM
    ax.scatter(pot, dp, s=sizes, c=colors, edgecolor="#222", linewidth=.4, zorder=3)
    # annotate the strongest selective + potent hits
    ann = sel[sel["selective_potent"]].head(annotate)
    if not len(ann):
        ann = sel.head(annotate)
    for _, row in ann.iterrows():
        ax.annotate(row["drug_name"], (row["potency_in"], row["delta_potency"]),
                    fontsize=7, color=RRB_MAROON, xytext=(4, 3), textcoords="offset points")
    ax.set_xlabel("potency in indication  =  -log10(median IC50 / uM)")
    ax.set_ylabel("selectivity  =  log-potency(in) - (rest)")
    ax.set_title(title, fontsize=9.5, loc="left")


def render_gdsc_report(rankings: dict, selectivities: dict, *, out_path: str | Path) -> Path:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(13.5, 9.2))
    fig.suptitle("WARHEAD - GDSC2 EC90 Potency & HCC / CRC Selectivity",
                 color=RRB_MAROON, fontsize=15, fontweight="bold", x=.5, y=.975)
    fig.text(.5, .945,
             "EC90 = concentration for 90% of fitted max effect (GDSC 2-param fit, bottom=0; "
             "% obs = curves reaching EC90 within the tested range)",
             ha="center", fontsize=8.5, color="#555")

    gs = fig.add_gridspec(2, 2, hspace=.34, wspace=.28, left=.11, right=.965, top=.90, bottom=.09)
    _rank_panel(fig.add_subplot(gs[0, 0]), rankings["CRC"], "A.  CRC (COREAD) - lowest EC90")
    _rank_panel(fig.add_subplot(gs[0, 1]), rankings["HCC"], "B.  HCC (LIHC) - lowest EC90")
    _sel_panel(fig.add_subplot(gs[1, 0]), selectivities["CRC"],
               "C.  CRC selectivity (maroon = selective & potent; amber = selective, weak)")
    _sel_panel(fig.add_subplot(gs[1, 1]), selectivities["HCC"], "D.  HCC selectivity")

    n_crc = int(selectivities["CRC"]["selective_potent"].sum()) if len(selectivities["CRC"]) else 0
    n_hcc = int(selectivities["HCC"]["selective_potent"].sum()) if len(selectivities["HCC"]) else 0
    fig.text(.11, .02,
             f"selective & potent (q<0.1, IC50<=1uM):  CRC = {n_crc}   HCC = {n_hcc}"
             "     |  selectivity on IC50 (fitted); ranking on EC90.  hypothesis-generating: "
             "cell-line tissue selectivity can reflect a co-enriched dependency",
             fontsize=7.5, color="#666")

    fig.savefig(out_path, format=out_path.suffix.lstrip(".").lower() or "pdf", dpi=150)
    plt.close(fig)
    return out_path
