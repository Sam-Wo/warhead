"""Interactive (Plotly) HCC/CRC selectivity plots as a self-contained HTML.

One page, a source dropdown (GDSC / PRISM / CTRP), CRC + HCC panels. Each point is
a compound; hover shows target, clinical status, IC50, selectivity and stats. The
plotly.js library is inlined so the file works offline in any browser.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

RRB_MAROON = "#6E1426"
_SEL_WEAK = "#C98A2E"
_GREY = "#B7BAC2"


def _colors(sel: pd.DataFrame):
    c, s = [], []
    for _, r in sel.iterrows():
        if r.get("selective_potent"):
            c.append(RRB_MAROON); s.append(13)
        elif r.get("selective"):
            c.append(_SEL_WEAK); s.append(9)
        else:
            c.append(_GREY); s.append(6)
    return c, s


def _customdata(sel: pd.DataFrame):
    def g(col):
        return sel[col] if col in sel else pd.Series([""] * len(sel))
    return np.column_stack([
        g("compound").astype(str), g("target").astype(str).str.slice(0, 40),
        g("clinical_phase").astype(str), g("median_ic50_in_nM").round(1).astype(str),
        g("delta_potency").round(2).astype(str), g("cliffs_delta").round(2).astype(str),
        g("q").map(lambda v: f"{v:.2g}"),
    ])


_HOVER = ("<b>%{customdata[0]}</b><br>target: %{customdata[1]}<br>"
          "clinical: %{customdata[2]}<br>median IC50(in): %{customdata[3]} nM<br>"
          "Δ potency: %{customdata[4]}   Cliff's δ: %{customdata[5]}   q: %{customdata[6]}"
          "<extra></extra>")


def render_selectivity_html(sel_by_source: dict, *, out_path, title="WARHEAD - screen selectivity (interactive)"):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    out_path = Path(out_path); out_path.parent.mkdir(parents=True, exist_ok=True)
    sources = list(sel_by_source)
    fig = make_subplots(rows=1, cols=2, subplot_titles=("CRC selectivity", "HCC selectivity"),
                        horizontal_spacing=.09)

    trace_source = []   # which source each trace belongs to
    for si, src in enumerate(sources):
        for col, ind in ((1, "CRC"), (2, "HCC")):
            sel = sel_by_source[src].get(ind)
            if sel is None or not len(sel):
                fig.add_trace(go.Scatter(x=[], y=[], mode="markers", showlegend=False), row=1, col=col)
            else:
                c, s = _colors(sel)
                fig.add_trace(go.Scatter(
                    x=sel["potency_in"], y=sel["delta_potency"], mode="markers",
                    marker=dict(color=c, size=s, line=dict(width=.5, color="#333")),
                    customdata=_customdata(sel), hovertemplate=_HOVER, showlegend=False,
                    visible=(si == 0)), row=1, col=col)
            trace_source.append(si)

    for col in (1, 2):
        fig.add_hline(y=0, line=dict(color="#999", width=1, dash="dot"), row=1, col=col)
        fig.update_xaxes(title_text="potency in indication = -log10(median IC50 / nM)", row=1, col=col)
    fig.update_yaxes(title_text="selectivity = log-potency(in) - (rest)", row=1, col=1)

    buttons = []
    for si, src in enumerate(sources):
        vis = [ts == si for ts in trace_source]
        buttons.append(dict(label=src, method="update",
                            args=[{"visible": vis}, {"title.text": f"{title}  -  {src}"}]))
    fig.update_layout(
        title=dict(text=f"{title}  -  {sources[0]}", font=dict(color=RRB_MAROON, size=18)),
        updatemenus=[dict(buttons=buttons, direction="down", x=0.0, xanchor="left",
                          y=1.16, yanchor="top", showactive=True)],
        template="plotly_white", height=560, margin=dict(t=110),
        annotations=list(fig.layout.annotations) + [dict(
            text="maroon = selective &amp; potent · amber = selective · grey = rest.  "
                 "hover for compound / target / clinical status.",
            showarrow=False, xref="paper", yref="paper", x=0, y=-0.22, font=dict(size=11, color="#666"))],
    )
    fig.write_html(str(out_path), include_plotlyjs="inline", full_html=True)
    return out_path
