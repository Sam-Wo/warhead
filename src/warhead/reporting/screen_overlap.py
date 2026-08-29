"""Compound-space overlap across the three dose-response screens.

Different screens name the same drug differently (salts, hyphenation, synonyms,
casing), so overlap is computed on a normalised name (lower-cased, punctuation and
parentheticals stripped). That recovers most true matches but is still a LOWER
BOUND - a synonym pair (e.g. an INN vs a code name) will not match. Identity by
structure (InChIKey) would be exact, but the public curve tables ship names only.

Public: `norm_name`, `tested_sets`, `three_set_counts`, `render_overlap_venn`.
"""
from __future__ import annotations

import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402

RRB_MAROON = "#6E1426"
# the three comparable EC90/IC50 dose-response screens (PDXE in-vivo / NCI-60
# z-scored live in a different space and are reported separately)
DR_SOURCES = ("GDSC2", "PRISM Repurposing (secondary)", "CTRP v2")
_VENN_COLORS = {"GDSC2": "#6E1426", "PRISM Repurposing (secondary)": "#2E7D6B",
                "CTRP v2": "#C98A2E"}


def norm_name(name) -> str:
    s = re.sub(r"\s*\(.*?\)\s*", "", str(name).lower())   # drop parentheticals
    return re.sub(r"[^a-z0-9]+", "", s)                   # strip spaces/punct


def tested_sets(cans: dict, sources=DR_SOURCES) -> dict:
    """{full_source_name: set(normalised compound names ever assayed in that screen)}."""
    return {s: {norm_name(c) for c in cans[s]["compound"].dropna().unique()}
            for s in sources if s in cans}


def three_set_counts(a: set, b: set, c: set) -> dict:
    """7 disjoint region sizes for a 3-set Venn (a=A only, ab=A&B only, ...)."""
    return {"a": len(a - b - c), "b": len(b - a - c), "c": len(c - a - b),
            "ab": len((a & b) - c), "ac": len((a & c) - b), "bc": len((b & c) - a),
            "abc": len(a & b & c)}


def render_overlap_venn(cans: dict, *, out_path, sources=DR_SOURCES) -> Path:
    """Hand-drawn 3-circle Venn of the compound libraries + a triple-overlap list."""
    out_path = Path(out_path); out_path.parent.mkdir(parents=True, exist_ok=True)
    ts = tested_sets(cans, sources)
    names = list(ts)
    A, B, C = (ts[names[0]], ts[names[1]], ts[names[2]])
    k = three_set_counts(A, B, C)
    tot = {names[0]: len(A), names[1]: len(B), names[2]: len(C)}

    fig = plt.figure(figsize=(12.5, 7.4))
    ax = fig.add_axes([0.0, 0.0, 0.62, 1.0]); ax.set_xlim(-2, 2); ax.set_ylim(-1.9, 2.0)
    ax.set_aspect("equal"); ax.axis("off")
    fig.suptitle("WARHEAD - compound-library overlap across the dose-response screens",
                 color=RRB_MAROON, fontsize=14, fontweight="bold", x=0.5, y=0.98)

    r = 1.05
    ctr = {names[0]: (-0.52, 0.34), names[1]: (0.52, 0.34), names[2]: (0.0, -0.55)}
    for s in names:
        ax.add_patch(Circle(ctr[s], r, facecolor=_VENN_COLORS[s], edgecolor=_VENN_COLORS[s],
                            alpha=0.28, lw=1.6, zorder=1))
    # screen labels + totals, pulled outside the circles
    ax.text(-1.35, 1.45, f"{names[0].split()[0]}\nn={tot[names[0]]}", color=_VENN_COLORS[names[0]],
            fontsize=11, fontweight="bold", ha="center")
    ax.text(1.35, 1.45, f"PRISM\nn={tot[names[1]]}", color=_VENN_COLORS[names[1]],
            fontsize=11, fontweight="bold", ha="center")
    ax.text(0.0, -1.72, f"CTRP v2   n={tot[names[2]]}", color=_VENN_COLORS[names[2]],
            fontsize=11, fontweight="bold", ha="center")

    def cnt(x, y, v):
        ax.text(x, y, str(v), fontsize=14, fontweight="bold", ha="center", va="center", color="#222")
    cnt(-0.95, 0.62, k["a"]); cnt(0.95, 0.62, k["b"]); cnt(0.0, -1.02, k["c"])
    cnt(0.0, 0.82, k["ab"]); cnt(-0.72, -0.32, k["ac"]); cnt(0.72, -0.32, k["bc"])
    cnt(0.0, 0.02, k["abc"])

    # ---- right column: what the overlap means + triple-shared compounds ----
    disp = {}                                            # normalised -> a real display name
    for s in names:
        for c in cans[s]["compound"].dropna().unique():
            disp.setdefault(norm_name(c), c)
    shared = sorted((disp[x] for x in (A & B & C)), key=str.lower)
    tx = 0.65
    fig.text(tx, 0.90, f"Tested in all three screens: {k['abc']}", fontsize=11,
             fontweight="bold", color=RRB_MAROON)
    fig.text(tx, 0.865, "(normalised-name match - a lower bound; synonyms/salts split some pairs)",
             fontsize=7.5, color="#777")
    col = shared[:34]
    half = (len(col) + 1) // 2
    for j, chunk in enumerate((col[:half], col[half:])):
        fig.text(tx + j * 0.17, 0.83, "\n".join(f"- {c[:22]}" for c in chunk),
                 fontsize=7.3, color="#333", va="top", family="monospace")
    if len(shared) > 34:
        fig.text(tx, 0.10, f"... and {len(shared) - 34} more", fontsize=7.5, color="#777")
    fig.text(tx, 0.045, "Overlap is what lets a hit REPLICATE across screens; a compound in only\n"
             "one circle has no independent confirmation (and cannot be cross-checked).",
             fontsize=7.6, color="#555", va="bottom")

    fig.savefig(out_path, format=out_path.suffix.lstrip(".").lower() or "pdf", dpi=150)
    plt.close(fig)
    return out_path
