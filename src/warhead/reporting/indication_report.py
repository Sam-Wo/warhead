"""Generic single-indication potency + selectivity report, one screen per row.

Reuses the rank / selectivity panels from reporting.screen (the CRC/HCC report) but
for an arbitrary indication across the screens that cover it - e.g. AML on CTRP and
GDSC (PRISM's secondary subset has no blood lines).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from .screen import _PHASE_COLOR, RRB_MAROON, _rank_panel, _sel_panel  # noqa: E402


def render_indication_report(indication: str, panels: list, *, out_path, caveat: str = "") -> Path:
    """panels: list of dicts {label, rank, sel, n_lines, n_total}. Renders a rank
    panel per screen (top-20 lowest EC90 in the indication) then a row of selectivity
    panels (indication vs the rest of that screen's panel). `caveat` overrides the
    bottom methodology note."""
    out_path = Path(out_path); out_path.parent.mkdir(parents=True, exist_ok=True)
    n = len(panels)
    fig = plt.figure(figsize=(13, 4.9 * n + 4.6))
    fig.suptitle(f"WARHEAD - {indication} EC90 potency & selectivity across public screens",
                 color=RRB_MAROON, fontsize=16, fontweight="bold", x=.5, y=.99)
    covered = ", ".join(f"{p['label']} ({p['n_lines']} {indication} lines / {p['n_total']})" for p in panels)
    fig.text(.5, .972, f"Screens with {indication} coverage: {covered}.  "
             "bar colour = clinical phase; hatched = EC90 extrapolated; Emax = residual viability.",
             ha="center", fontsize=8.6, color="#555")

    gs = fig.add_gridspec(n + 1, n, height_ratios=[1.35] * n + [1.0], hspace=.34, wspace=.26,
                          left=.16, right=.965, top=.945, bottom=.065)
    for i, p in enumerate(panels):
        _rank_panel(fig.add_subplot(gs[i, :]), p["rank"],
                    f"{chr(65 + i)}.  {p['label']} - top 20 lowest EC90 in {indication}   "
                    f"({p['n_lines']} lines)")
    for i, p in enumerate(panels):
        _sel_panel(fig.add_subplot(gs[n, i]), p["sel"],
                   f"{chr(65 + n + i)}.  {p['label']} {indication} selectivity"
                   + ("  (maroon = selective & potent)" if i == 0 else ""))

    handles = [Patch(facecolor=_PHASE_COLOR["Launched"], label="Launched"),
               Patch(facecolor=_PHASE_COLOR["Phase 3"], label="Phase 2/3"),
               Patch(facecolor=_PHASE_COLOR["Phase 1"], label="Phase 1"),
               Patch(facecolor=_PHASE_COLOR["Preclinical"], label="Preclinical"),
               Patch(facecolor="white", edgecolor="#3a0a13", hatch="////", label="EC90 extrapolated")]
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=8, frameon=False,
               bbox_to_anchor=(.5, .012))
    default = (f"selectivity anchored on fitted IC50 ({indication} vs the rest of each panel); "
               "ranking on EC90 (Emax-filtered). Hypothesis-generating: tissue selectivity can reflect a "
               "co-enriched dependency, not a payload-ready ADC.")
    fig.text(.16, .002, caveat or default, fontsize=7.5, color="#666", wrap=True)
    fig.savefig(out_path, format=out_path.suffix.lstrip(".").lower() or "pdf", dpi=150)
    plt.close(fig)
    return out_path
