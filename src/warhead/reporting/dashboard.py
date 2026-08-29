"""reports/screens_dashboard.html - one interactive page, a tab per screen.

Each tab: the screen's metadata, an interactive CRC/HCC selectivity plot (hover for
compound/target/clinical), and the top-20 ranking table. plotly.js is inlined once
so the file is fully offline. Dose-response screens (GDSC/PRISM/CTRP) show EC90/IC50
selectivity + ranking; PDXE shows in-vivo CRC response; NCI-60 shows CRC-selectivity.
"""
from __future__ import annotations

import base64
import html
from pathlib import Path

import numpy as np
import pandas as pd

from warhead.reporting.screen_overlap import norm_name

RRB = "#6E1426"
_WEAK = "#C98A2E"
_GREY = "#B7BAC2"

# short letter per dose-response screen, for the cross-screen coverage badge
_COV_ORDER = [("GDSC2", "G"), ("PRISM Repurposing (secondary)", "P"), ("CTRP v2", "C")]


def _coverage_badges(compound, tested):
    """Small G/P/C badges: bold-maroon where the compound is in that screen's
    library, faint grey where it was not tested there."""
    if not tested:
        return ""
    k = norm_name(compound)
    spans = []
    for src, letter in _COV_ORDER:
        on = k in tested.get(src, set())
        style = ("color:#6E1426;font-weight:700" if on
                 else "color:#cfcfcf;text-decoration:line-through")
        spans.append(f'<span style="{style}" title="{html.escape(src)}: '
                     f'{"tested" if on else "not tested"}">{letter}</span>')
    return "&nbsp;".join(spans)


def _img_div(png_bytes, alt):
    b64 = base64.b64encode(png_bytes).decode("ascii")
    return (f'<img src="data:image/png;base64,{b64}" alt="{html.escape(alt)}" '
            f'style="max-width:100%;height:auto;border:1px solid #eadfe2;border-radius:8px">')


def _sel_traces(sel, ind):
    import plotly.graph_objects as go
    if sel is None or not len(sel):
        return go.Scatter(x=[], y=[], mode="markers", name=ind, showlegend=False)
    col, size = [], []
    for _, r in sel.iterrows():
        if r.get("selective_potent"): col.append(RRB); size.append(12)
        elif r.get("selective"): col.append(_WEAK); size.append(8)
        else: col.append(_GREY); size.append(6)
    cd = np.column_stack([sel["compound"].astype(str),
                          sel.get("target", pd.Series([""] * len(sel))).astype(str).str.slice(0, 40),
                          sel.get("clinical_phase", pd.Series([""] * len(sel))).astype(str),
                          sel["median_ic50_in_nM"].round(1).astype(str),
                          sel["delta_potency"].round(2).astype(str),
                          sel["q"].map(lambda v: f"{v:.2g}")])
    return go.Scatter(
        x=sel["potency_in"], y=sel["delta_potency"], mode="markers",
        marker=dict(color=col, size=size, line=dict(width=.5, color="#333")),
        customdata=cd, showlegend=False,
        hovertemplate=("<b>%{customdata[0]}</b><br>target: %{customdata[1]}<br>"
                       "clinical: %{customdata[2]}<br>median IC50: %{customdata[3]} nM<br>"
                       "Δpotency %{customdata[4]}  q %{customdata[5]}<extra></extra>"))


def _selectivity_div(sel_crc, sel_hcc, label):
    from plotly.subplots import make_subplots
    fig = make_subplots(rows=1, cols=2, subplot_titles=("CRC selectivity", "HCC selectivity"),
                        horizontal_spacing=.1)
    fig.add_trace(_sel_traces(sel_crc, "CRC"), row=1, col=1)
    fig.add_trace(_sel_traces(sel_hcc, "HCC"), row=1, col=2)
    for c in (1, 2):
        fig.add_hline(y=0, line=dict(color="#999", width=1, dash="dot"), row=1, col=c)
        fig.update_xaxes(title_text="-log10(median IC50 / nM)", row=1, col=c)
    fig.update_yaxes(title_text="selectivity = log-potency(in)-(rest)", row=1, col=1)
    fig.update_layout(template="plotly_white", height=430, margin=dict(t=40, b=50))
    return fig.to_html(include_plotlyjs=False, full_html=False, div_id=f"sel_{label}",
                       default_width="100%")


def _vline(fig, x, color, row, col):
    """Vertical marker as a 2-point trace (robust on a log x-axis, unlike shapes)."""
    import plotly.graph_objects as go
    if x is None or not np.isfinite(x):
        return
    fig.add_trace(go.Scatter(x=[x, x], y=[-0.05, 1.2], mode="lines",
                             line=dict(color=color, width=1, dash="dash"),
                             hoverinfo="skip", showlegend=False), row=row, col=col)


def _curve_titles(summary):
    return [f"{r['compound'][:20]} ({str(r['target'])[:14] if pd.notna(r['target']) else 'n/a'})"
            for _, r in summary.reset_index(drop=True).iterrows()]


def _curves_band_div(pooled, summary, div_id, ncol=4, markers=True):
    """Interactive small-multiples of median + IQR-band dose-response curves.

    markers=True draws lines+markers (CTRP: real measured bins); markers=False
    draws a smooth line (PRISM: cross-line median of the per-line 4PL fits). Both
    get the shaded IQR band and IC50/EC90 dashed markers, so the two screens read
    the same."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    s = summary.reset_index(drop=True); n = len(s)
    nrow = int(np.ceil(n / ncol))
    mode = "lines+markers" if markers else "lines"
    fig = make_subplots(rows=nrow, cols=ncol, subplot_titles=_curve_titles(s),
                        vertical_spacing=0.055, horizontal_spacing=0.045)
    for i, r in s.iterrows():
        rr, cc = i // ncol + 1, i % ncol + 1
        cur = pooled[pooled["compound"] == r["compound"]].sort_values("conc_uM")
        x = (cur["conc_uM"] * 1e3).tolist()
        fig.add_trace(go.Scatter(x=x + x[::-1], y=cur["q3"].tolist() + cur["q1"].tolist()[::-1],
                                 fill="toself", fillcolor="rgba(110,20,38,0.12)", line=dict(width=0),
                                 hoverinfo="skip", showlegend=False), row=rr, col=cc)
        fig.add_trace(go.Scatter(x=x, y=cur["median"], mode=mode,
                                 line=dict(color=RRB, width=1.6),
                                 marker=dict(size=4, color=RRB) if markers else None,
                                 hovertemplate="%{x:.0f} nM<br>viab %{y:.2f}<extra></extra>",
                                 showlegend=False), row=rr, col=cc)
        _vline(fig, r["ic50_uM"] * 1e3, "#222", rr, cc)
        _vline(fig, r["ec90_uM"] * 1e3, RRB, rr, cc)
        fig.update_xaxes(type="log", row=rr, col=cc)
        fig.update_yaxes(range=[-0.05, 1.2], row=rr, col=cc)
    fig.update_annotations(font_size=10)
    fig.update_xaxes(title_text="", tickfont_size=8)
    fig.update_yaxes(tickfont_size=8)
    fig.update_layout(template="plotly_white", height=210 * nrow, margin=dict(t=34, b=20, l=30, r=10))
    return fig.to_html(include_plotlyjs=False, full_html=False, div_id=div_id, default_width="100%")


def _bar_div(labels, values, colors, title, xtitle, div_id, height=460):
    import plotly.graph_objects as go
    fig = go.Figure(go.Bar(x=values, y=labels, orientation="h",
                           marker_color=colors, hovertemplate="%{y}: %{x}<extra></extra>"))
    fig.update_layout(template="plotly_white", height=height, margin=dict(t=40, l=180),
                      title=title, xaxis_title=xtitle)
    return fig.to_html(include_plotlyjs=False, full_html=False, div_id=div_id, default_width="100%")


def _table(df, cols, headers, fmts):
    th = "".join(f"<th>{html.escape(h)}</th>" for h in headers)
    rows = []
    for _, r in df.iterrows():
        tds = []
        for c, f in zip(cols, fmts):
            v = r.get(c)
            tds.append(f"<td>{html.escape(f(v))}</td>")
        rows.append("<tr>" + "".join(tds) + "</tr>")
    return (f'<table class="rank"><thead><tr>{th}</tr></thead>'
            f'<tbody>{"".join(rows)}</tbody></table>')


def _table_cov(r, tested):
    """Top-N dose-response ranking table with an optional cross-screen 'also in'
    coverage column (G/P/C badges)."""
    head = ["compound", "target", "IC50 nM", "EC90 nM", "Emax", "clinical"]
    if tested:
        head.append("also in")
    th = "".join(f"<th>{html.escape(h)}</th>" for h in head)
    rows = []
    for _, row in r.iterrows():
        em = row.get("median_emax")
        em_s = "-" if pd.isna(em) else f"{em:.2f}"
        ph = row.get("clinical_phase")
        ph_s = "-" if (ph is None or pd.isna(ph) or str(ph).lower() in ("nan", "<na>")) else str(ph)
        tds = [f'<td>{html.escape(str(row.get("compound")))}</td>',
               f'<td>{html.escape(str(row.get("target"))[:40])}</td>',
               f'<td>{_num(row.get("median_ic50_nM"))}</td>',
               f'<td>{_num(row.get("median_ec90_nM"))}</td>',
               f'<td>{em_s}</td>',
               f'<td>{html.escape(ph_s)}</td>']
        if tested:
            tds.append(f'<td style="font-family:monospace;letter-spacing:1px">'
                       f'{_coverage_badges(row.get("compound"), tested)}</td>')
        rows.append("<tr>" + "".join(tds) + "</tr>")
    return (f'<table class="rank"><thead><tr>{th}</tr></thead>'
            f'<tbody>{"".join(rows)}</tbody></table>')


def _num(v):
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return "-"
    v = float(v)
    return f"{v/1000:.1f}k" if v >= 1000 else (f"{v:.0f}" if v >= 10 else f"{v:.2g}")


def _meta_card(m):
    fields = [("compounds", m["compounds"]), ("cell lines", m["cell_lines"]),
              ("CRC lines", m["CRC_lines"]), ("HCC lines", m["HCC_lines"]),
              ("dose range", m["dose_range"]), ("# doses", m["n_doses"]),
              ("metrics", m["metrics"]), ("assay", m["assay"])]
    items = "".join(f'<div class="mi"><span class="k">{html.escape(str(k))}</span>'
                    f'<span class="v">{html.escape(str(v))}</span></div>' for k, v in fields)
    return f'<div class="metacard">{items}</div>'


def _overlap_table(counts, totals):
    rows = [("GDSC2 only", counts["a"]), ("PRISM only", counts["b"]), ("CTRP v2 only", counts["c"]),
            ("GDSC2 &amp; PRISM", counts["ab"]), ("GDSC2 &amp; CTRP", counts["ac"]),
            ("PRISM &amp; CTRP", counts["bc"]), ("all three", counts["abc"])]
    body = "".join(f"<tr><td>{r}</td><td>{v}</td></tr>" for r, v in rows)
    return (f'<table class="rank" style="max-width:360px"><thead><tr><th>region</th>'
            f'<th>compounds</th></tr></thead><tbody>{body}</tbody></table>')


def render_dashboard(screens, *, out_path, tested=None, venn_png=None, overlap_counts=None,
                     overlap_totals=None):
    from plotly.offline import get_plotlyjs
    out_path = Path(out_path); out_path.parent.mkdir(parents=True, exist_ok=True)
    tabs_btn, tabs_div, first = [], [], None
    cov_h = ['<th title="also tested in: G=GDSC2 P=PRISM C=CTRP v2">also in</th>'] if tested else []
    for sc in screens:
        lab = sc["label"]; sid = "".join(ch for ch in lab if ch.isalnum())
        first = first or sid
        tabs_btn.append(f'<button class="tab-btn" onclick="showTab(\'{sid}\')" id="btn_{sid}">{html.escape(lab)}</button>')
        body = [f"<h2>{html.escape(lab)}</h2>"]
        if sc.get("meta"):
            body.append(_meta_card(sc["meta"]))
        if sc["type"] == "dr":
            body.append('<h3>HCC / CRC selectivity (hover for details)</h3>')
            body.append(_selectivity_div(sc["sel"]["CRC"], sc["sel"]["HCC"], sid))
            for ind in ("CRC", "HCC"):
                r = sc["rank"][ind].head(20)
                body.append(f"<h3>Top 20 by EC90 - {ind}</h3>")
                body.append(_table_cov(r, tested))
            cv = sc.get("curves")
            if cv is not None:
                body.append(f'<h3>Top-20 dose-response curves - {cv["note"]}</h3>')
                body.append('<p style="color:#666;font-size:12px;margin:0 0 6px">'
                            'median (line) + IQR across cell lines (shaded); hover for viability at '
                            'each dose; dashed = IC50 (black) &amp; EC90 (maroon).</p>')
                body.append(_curves_band_div(cv["pooled"], cv["summary"], f"curves_{sid}",
                                             markers=cv.get("markers", True)))
        elif sc["type"] == "overlap":
            body.append('<h3>Compound-library overlap across the dose-response screens</h3>')
            if venn_png is not None:
                body.append(_img_div(venn_png, "compound overlap Venn"))
            if overlap_counts is not None:
                body.append("<h3>Region sizes</h3>")
                body.append(_overlap_table(overlap_counts, overlap_totals or {}))
            body.append('<p style="color:#666;font-size:12.5px;max-width:760px">Overlap is computed on a '
                        'normalised compound name (a lower bound - synonyms/salts split some true pairs). '
                        'A compound tested in only one screen has no independent confirmation; the '
                        '"also in" column on each screen tab flags which of the other two libraries '
                        'contain each ranked compound.</p>')
        elif sc["type"] == "pdxe":
            rk = sc["ranking"]
            colors = [RRB if v < 0 else _GREY for v in rk["median_response"]]
            body.append(_bar_div(rk["Treatment"].tolist()[::-1], rk["median_response"].tolist()[::-1],
                                 colors[::-1], "CRC single-agent response (neg = shrinkage)",
                                 "median best-avg tumour response (%)", f"bar_{sid}", height=560))
            body.append("<h3>CRC treatments</h3>")
            body.append(_table(rk, ["Treatment", "target", "n_models", "median_response", "pct_objective_response"],
                               ["treatment", "target", "n models", "median response %", "ORR %"],
                               [str, lambda v: str(v)[:34], lambda v: str(int(v)),
                                lambda v: f"{v:.1f}", lambda v: f"{v:.0f}"]))
        elif sc["type"] == "nci60":
            an = sc["selectivity"].sort_values("delta_z", ascending=False).head(25).copy()
            # most NCI-60 compounds are unnamed research entries - fall back to the
            # NSC accession so every bar has a unique, identifiable label (two "-"
            # names otherwise collide onto one y-tick and read as a double bar)
            def _nci_label(row):
                d = str(row["drug"]).strip()
                return d[:30] if d and d.lower() not in ("-", "nan", "none") else f"NSC {row['NSC']}"
            an["label"] = an.apply(_nci_label, axis=1)
            body.append(_bar_div(an["label"].tolist()[::-1],
                                 an["delta_z"].tolist()[::-1], [RRB] * len(an),
                                 "CRC(colon)-selective (z colon - rest)", "delta z", f"bar_{sid}", height=620))
            body.append("<h3>Top CRC-selective (annotated)</h3>")
            body.append(_table(an, ["label", "MOA", "FDA", "delta_z"],
                               ["compound", "MoA", "clinical", "Δz colon"],
                               [str, lambda v: str(v)[:34], str, lambda v: f"{v:.2f}"]))
        tabs_div.append(f'<div id="{sid}" class="tab-pane">{"".join(body)}</div>')

    css = """
    body{font-family:'IBM Plex Sans',system-ui,sans-serif;margin:0;color:#17181d;background:#faf9fa}
    h1{color:#6E1426;padding:18px 26px 6px;margin:0}
    .sub{padding:0 26px 10px;color:#666;font-size:13px}
    .tabbar{display:flex;gap:4px;padding:0 20px;border-bottom:2px solid #eadfe2;flex-wrap:wrap}
    .tab-btn{border:none;background:none;padding:11px 18px;font-size:14px;font-weight:600;color:#888;cursor:pointer;border-bottom:3px solid transparent}
    .tab-btn.active{color:#6E1426;border-bottom-color:#6E1426}
    .tab-pane{display:none;padding:16px 26px 40px}
    .tab-pane.active{display:block}
    h2{color:#6E1426;margin:8px 0 12px} h3{margin:22px 0 8px;font-size:15px;color:#333}
    .metacard{display:flex;flex-wrap:wrap;gap:10px 26px;background:#fff;border:1px solid #eadfe2;border-radius:10px;padding:14px 18px;box-shadow:0 1px 2px rgba(0,0,0,.04)}
    .mi{display:flex;flex-direction:column} .mi .k{font-size:11px;letter-spacing:.05em;text-transform:uppercase;color:#999}
    .mi .v{font-size:14px;font-weight:600;color:#17181d}
    table.rank{border-collapse:collapse;font-size:12.5px;margin:6px 0 14px;width:100%;max-width:1100px}
    table.rank th{background:#6E1426;color:#fff;text-align:left;padding:6px 10px;position:sticky;top:0;cursor:pointer}
    table.rank td{padding:4px 10px;border-bottom:1px solid #eee}
    table.rank tbody tr:nth-child(even){background:#faf5f6}
    """
    js = """
    function showTab(id){
      document.querySelectorAll('.tab-pane').forEach(p=>p.classList.remove('active'));
      document.querySelectorAll('.tab-btn').forEach(b=>b.classList.remove('active'));
      document.getElementById(id).classList.add('active');
      document.getElementById('btn_'+id).classList.add('active');
      document.querySelectorAll('#'+id+' .js-plotly-plot').forEach(d=>{try{Plotly.Plots.resize(d);}catch(e){}});
    }
    document.querySelectorAll('table.rank th').forEach((th,ci)=>th.addEventListener('click',()=>{
      const tb=th.closest('table').querySelector('tbody');
      const rows=[...tb.rows];
      const num=v=>{v=v.replace('k','e3').replace(/[^0-9.eE-]/g,'');const f=parseFloat(v);return isNaN(f)?null:f;};
      th._asc=!th._asc;
      rows.sort((a,b)=>{let x=a.cells[ci].innerText,y=b.cells[ci].innerText;const nx=num(x),ny=num(y);
        if(nx!==null&&ny!==null)return th._asc?nx-ny:ny-nx; return th._asc?x.localeCompare(y):y.localeCompare(x);});
      rows.forEach(r=>tb.appendChild(r));
    }));
    """
    doc = (f"<!doctype html><html><head><meta charset='utf-8'>"
           f"<link rel='preconnect' href='https://fonts.googleapis.com'>"
           f"<link rel='stylesheet' href='https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;600&display=swap'>"
           f"<style>{css}</style><script>{get_plotlyjs()}</script></head><body>"
           f"<h1>WARHEAD - public drug-screen dashboard</h1>"
           f"<div class='sub'>EC90/IC50 potency, HCC/CRC selectivity and clinical/ADC context across five public screens. "
           f"Click a tab; hover points; click a table header to sort.</div>"
           f"<div class='tabbar'>{''.join(tabs_btn)}</div>{''.join(tabs_div)}"
           f"<script>{js}\nshowTab('{first}');</script></body></html>")
    out_path.write_text(doc, encoding="utf-8")
    return out_path
