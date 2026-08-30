"""reports/conjugation_scorecard_*.pdf - Phase-A conjugation-suitability scorecard.

Re-scores the top free-drug hits on the questions an ADC payload must answer with
data already in hand: the G1 sub-nM payload-potency bar (+ the size of the potency
gap), pan-panel G2b proliferation-independence, and the curated known-payload / DLT
annotation. The chemical gates (G4/G5) and G2a efflux are Phase B (need structures /
expression) and are shown as pending columns so the reader sees the whole ledger.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402

RRB = "#6E1426"
_PASS = "#2E7D6B"
_FAIL = "#B7BAC2"
_MITO = "#C98A2E"


def _fmt_ic50(v):
    if pd.isna(v):
        return "-"
    return f"{v/1000:.1f}µM" if v >= 1000 else f"{v:.0f}nM" if v >= 10 else f"{v:.1f}nM"


def render_scorecard(df: pd.DataFrame, *, out_path, indication="CRC", screen="CTRP v2") -> Path:
    out_path = Path(out_path); out_path.parent.mkdir(parents=True, exist_ok=True)
    n = len(df)
    fig = plt.figure(figsize=(15.5, 0.42 * n + 3.0))
    ax = fig.add_axes([0.01, 0.02, 0.98, 0.90]); ax.axis("off")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    fig.suptitle(f"WARHEAD - conjugation-suitability scorecard ({screen}, top {indication})",
                 color=RRB, fontsize=15, fontweight="bold", y=.985)
    fig.text(0.01, 0.945,
             "Phase A - delivery survivability from data in hand. The ADC payload bar (G1): IC50 ≤ 1 nM in "
             "≥20% of lines AND median Emax < 0.15 (complete kill).",
             fontsize=8.7, color="#333")
    fig.text(0.01, 0.928,
             "'gap' = log10 units the median IC50 sits ABOVE 1 nM = the potency a medchem campaign must close.  "
             "G5 = conjugatable handle present (necessary, not SAR-verified);  G4 = physchem in the bystander window.",
             fontsize=8.7, color="#333")

    # columns
    xc = {"compound": 0.0, "target": 0.155, "ic50": 0.34, "gap": 0.40, "emax": 0.46,
          "g1": 0.52, "prolif": 0.585, "efflux": 0.70, "handle": 0.76, "bystd": 0.82,
          "payload": 0.885}
    hy = 0.90
    hdrs = [("compound", "compound"), ("target", "target"), ("ic50", "med IC50"),
            ("gap", "gap↑"), ("emax", "Emax"), ("g1", "G1 bar"),
            ("prolif", "G2b prolif"), ("efflux", "G2a"), ("handle", "G5"),
            ("bystd", "G4"), ("payload", "known payload?")]
    has_chem = "g5_handle" in df.columns
    grey_keys = ("efflux",) if has_chem else ("efflux", "handle", "bystd")
    for key, lab in hdrs:
        ax.text(xc[key], hy, lab, fontsize=7.6, fontweight="bold",
                color="#bbb" if key in grey_keys else RRB)
    if has_chem:
        ax.text(xc["efflux"], hy - 0.028, "G2a pending", fontsize=6.4, color="#bbb")
    else:
        ax.text(xc["efflux"], hy - 0.028, "― Phase B (needs structures / expression) ―",
                fontsize=6.6, color="#bbb")

    gnorm = Normalize(vmin=0, vmax=2.5)
    row_h = (hy - 0.05) / max(n, 1)
    for i, r in df.reset_index(drop=True).iterrows():
        y = hy - 0.05 - (i + 0.5) * row_h
        ax.text(xc["compound"], y, str(r["compound"])[:24], fontsize=7.6, va="center", fontweight="bold")
        ax.text(xc["target"], y, str(r.get("target") or "")[:26], fontsize=6.2, va="center",
                color="#555", family="monospace")
        ax.text(xc["ic50"], y, _fmt_ic50(r["median_ic50_nM"]), fontsize=7, va="center", ha="center")
        # potency gap - heat (bigger gap = redder background)
        gap = r.get("potency_gap_log10")
        if pd.notna(gap):
            col = plt.cm.OrRd(gnorm(max(gap, 0)))
            ax.add_patch(plt.Rectangle((xc["gap"] - 0.026, y - row_h * .42), 0.052, row_h * .84,
                                       facecolor=col, edgecolor="#eee", lw=.3, zorder=1))
            lum = 0.299 * col[0] + 0.587 * col[1] + 0.114 * col[2]
            ax.text(xc["gap"], y, f"{gap:.1f}", fontsize=6.6, va="center", ha="center",
                    color="white" if lum < 0.5 else "#222", zorder=2)
        em = r.get("median_emax")
        ax.text(xc["emax"], y, f"{em:.2f}" if pd.notna(em) else "-", fontsize=6.8, va="center", ha="center",
                color=_PASS if (pd.notna(em) and em < 0.15) else "#444")
        ax.text(xc["g1"], y, "✓ pass" if r["g1_potency_pass"] else "✗ fail",
                fontsize=7, va="center", ha="center", fontweight="bold",
                color=_PASS if r["g1_potency_pass"] else _FAIL)
        pc = str(r.get("prolif_class"))
        pcol = {"independent": _PASS, "mitotic-dependent": _MITO}.get(pc, "#999")
        ax.text(xc["prolif"], y, pc, fontsize=6.6, va="center", ha="center", color=pcol,
                fontweight="bold" if pc == "mitotic-dependent" else "normal")
        ax.text(xc["efflux"], y, "·", fontsize=9, va="center", ha="center", color="#ccc")
        if has_chem:
            for key, col in (("handle", "g5_handle"), ("bystd", "g4_bystander")):
                v = r.get(col)
                if pd.isna(v):
                    ax.text(xc[key], y, "n/a", fontsize=6, va="center", ha="center", color="#ccc")
                else:
                    ax.text(xc[key], y, "✓" if v else "✗", fontsize=8, va="center", ha="center",
                            fontweight="bold", color=_PASS if v else _FAIL)
        else:
            for key in ("handle", "bystd"):
                ax.text(xc[key], y, "·", fontsize=9, va="center", ha="center", color="#ccc")
        pay = str(r.get("adc_payload_status") or "")
        is_known = pay and "not a payload" not in pay.lower()
        ax.text(xc["payload"], y, ("● " if is_known else "") + (pay[:26] if pay else "–"),
                fontsize=6.3, va="center", color=RRB if is_known else "#999")

    fig.text(0.01, 0.018,
             "G2b prolif = pan-panel Spearman(log10 IC50, DepMap CRISPR growth rate), BH q<0.05; indicative "
             "(growth-rate proxy; efflux substrates e.g. taxanes are confounded - see G2a).", fontsize=6.9, color="#777")
    fig.text(0.01, 0.006,
             "Authoritative G2b = `warhead g2b-real`. No hit clears G1 as a free drug -> the output is a "
             "chemotype for a potency campaign, not a molecule to conjugate as-is.", fontsize=6.9, color="#777")
    fig.savefig(out_path, format=out_path.suffix.lstrip(".").lower() or "pdf", dpi=150)
    plt.close(fig)
    return out_path
