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
_GREEN = "#2E7D6B"

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


def _sel_pair_div(pairs, div_id, height=430):
    """Selectivity scatter for an arbitrary set of (subplot_title, sel_frame) pairs -
    e.g. one screen per column for a single indication."""
    from plotly.subplots import make_subplots
    fig = make_subplots(rows=1, cols=len(pairs), subplot_titles=[t for t, _ in pairs],
                        horizontal_spacing=.1)
    for i, (t, sel) in enumerate(pairs):
        fig.add_trace(_sel_traces(sel, t), row=1, col=i + 1)
        fig.add_hline(y=0, line=dict(color="#999", width=1, dash="dot"), row=1, col=i + 1)
        fig.update_xaxes(title_text="-log10(median IC50 / nM)", row=1, col=i + 1)
    fig.update_yaxes(title_text="selectivity = log-potency(in)-(rest)", row=1, col=1)
    fig.update_layout(template="plotly_white", height=height, margin=dict(t=40, b=50))
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


def _flag_cell(v, *, good_is_true=True):
    """✓/✗/n-a cell. good_is_true=False flips the sense (e.g. efflux substrate: True is bad)."""
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return '<td style="text-align:center;color:#ccc">n/a</td>'
    ok = bool(v) if good_is_true else (not bool(v))
    return (f'<td style="text-align:center;color:{_GREEN if ok else _WEAK};font-weight:700">'
            f'{"✓" if ok else "✗"}</td>')


def _scorecard_table(df):
    head = ["compound", "target", "med IC50", "gap↑", "G1", "G2a", "G2b", "G5", "G4", "G6", "known payload?"]
    th = "".join(f"<th>{html.escape(h)}</th>" for h in head)
    rows = []
    for _, r in df.iterrows():
        gap = r.get("potency_gap_log10")
        pay = str(r.get("adc_payload_status") or "")
        known = pay and "not a payload" not in pay.lower()
        tds = [
            f'<td>{html.escape(str(r["compound"]))}</td>',
            f'<td style="font-family:monospace;font-size:11px;color:#666">{html.escape(str(r.get("target"))[:24])}</td>',
            f'<td>{_num(r.get("median_ic50_nM"))} nM</td>',
            f'<td style="text-align:center">{"-" if pd.isna(gap) else f"{gap:.1f}"}</td>',
            _flag_cell(r.get("g1_potency_pass")),
            _flag_cell(r.get("g2a_substrate"), good_is_true=False),
            _flag_cell(r.get("g2b_independent")),
            _flag_cell(r.get("g5_handle")),
            _flag_cell(r.get("g4_bystander")),
            _flag_cell(r.get("g6_window_ok")),
            f'<td style="font-size:11px;color:{RRB if known else "#999"}">'
            f'{("● " if known else "") + html.escape(pay[:26]) if pay else "-"}</td>',
        ]
        rows.append("<tr>" + "".join(tds) + "</tr>")
    return (f'<table class="rank"><thead><tr>{th}</tr></thead>'
            f'<tbody>{"".join(rows)}</tbody></table>')


# ---- Cascade (gate reference) tab -----------------------------------------
# thresholds mirror config/gates.yaml; a "thr" value renders as a mono block, a
# "caveat" value as muted text, anything else as inline HTML. Symbols as entities.
_CASCADE_GATES = [
    {"id": "G1", "name": "Potency, done properly", "status": ("ok", "implemented + tested"),
     "asks": '<span class="q">Is the free drug potent enough to survive becoming an ADC?</span> An ADC '
             'delivers very little payload per cell (low DAR, inefficient internalisation and lysosomal '
             'escape), so a payload must be far more potent than an ordinary drug.',
     "rows": [("Method", "Refit every curve as a 4-parameter logistic with an <b>interval-censored (Tobit)</b> "
               "likelihood; keep IC50, Emax and Hill <b>separate</b>. Never use AUC. Assay floor/ceiling are "
               "censored, not imputed to the lowest dose."),
              ("Threshold", ("thr", "gate:  IC50 &le; <b>1.0 nM</b> in &ge; <b>20%</b> of lines\n"
                             "   AND  median Emax &lt; <b>0.15</b>  (complete kill)\n"
                             "refit: 4PL &middot; floor 0.02 &middot; ceiling 1.05 &middot; &ge;5 pts\n"
                             "       hill [0.3, 8.0] &middot; Emax [-0.2, 1.2]")),
              ("Data", "Any screen's per-line IC50 + Emax (GDSC / PRISM / CTRP&hellip;)."),
              ("Caveat", ("caveat", "Absolute IC50 is <b>not comparable across screens</b> &mdash; see below."))],
     "callout": ("Choosing a sensible potency cutoff",
                 'The <span class="mono">1 nM / 20%</span> bar is calibrated to <b>PRISM</b>, where the reference '
                 'payloads read true: exatecan <code>0.10 nM</code> (82% of lines sub-nM), SN-38 <code>0.93 nM</code>, '
                 'maytansinol <code>0.99 nM</code>. It must <b>not</b> be applied to absolute IC50 from a '
                 'differently-run screen &mdash; CTRP\'s 72&nbsp;h CellTiter-Glo reads a systematic <b>~4.7&times; '
                 'weaker</b> (SN-38 is <code>121 nM</code> there), so the same molecule fails a bar it clears '
                 'elsewhere. This is assay design, not biology. <b>The fix:</b> apply G1 <em>per screen, anchored '
                 'to a reference payload run in that same screen</em> &mdash; pass &asymp; within ~1 log of '
                 'exatecan/SN-38 in that screen. Distance-to-a-known-payload is portable; absolute nM is not.')},
    {"id": "G2", "name": "Does the mechanism survive delivery?", "status": None,
     "asks": "Three independent ways a potent free drug still fails once it is a conjugate: it gets pumped back "
             "out, it only works on fast-dividing cells, or its target offers no tumour-selective hook.",
     "subs": [
         {"name": "G2a &middot; Efflux dependence", "tag": "ABCB1 / ABCG2", "status": ("ok", "implemented + tested"),
          "rows": [("Asks", "Is it a P-gp / BCRP substrate &mdash; the shared failure mode of MMAE and DM1?"),
                   ("Method", "Regress the resistance axis (log10 IC50) on transporter expression across the panel. "
                    "A strong positive slope = IC50 rises with transporter level = substrate."),
                   ("Threshold", ("thr", "substrate if  std slope &ge; <b>0.25</b>  on either gene\n"
                                  "              AND  BH q &lt; <b>0.05</b>")),
                   ("Data", "DepMap expression (ABCB1, ABCG2) &times; screen IC50.")]},
         {"name": "G2b &middot; Proliferation independence", "tag": "the HCC lever", "status": ("ok", "implemented + report"),
          "rows": [("Asks", "Does potency depend on doubling time? A payload for a low-proliferative-index "
                    "indication (MASH-HCC) must kill cells that are not dividing."),
                   ("Method", "Regress log10 IC50 on doubling time. <b>Keep</b> compounds whose slope is not "
                    "distinguishable from zero. Auristatins steep-positive &rarr; fail; RNAPII / translation / "
                    "degraders flat &rarr; pass. Slow lines up-weighted (equal-width bins)."),
                   ("Threshold", ("thr", "pass if  BH q &gt; <b>0.05</b>  (slope not significant)\n"
                                  "std_slope_max <b>0.15</b> advisory &middot; &ge;25 lines &middot; balance, 5 bins")),
                   ("Data", "Screen IC50 &times; DepMap growth rate (doubling-time proxy).")]},
         {"name": "G2c &middot; Collateral-lethality scan", "tag": "generalises POLR2A", "status": ("ok", "implemented + report"),
          "rows": [("Asks", "Is there a payload <em>target</em> whose recurrent hemizygous loss in the tumour "
                    "opens a dependency the normal tissue does not share?"),
                   ("Method", "For every gene: does CN-loss in TCGA-COAD/READ or LIHC shift the DepMap Chronos "
                    "dependency leftward (Mann&ndash;Whitney, FDR)? Exclude pan-essentials. Anchor on POLR2A "
                    "(17p, with TP53); candidate ME2 (18q, with SMAD4)."),
                   ("Threshold", ("thr", "CN loss &lt; <b>-0.3</b> &middot; recurrence &ge; <b>0.20</b> of tumours\n"
                                  "pan-essential cut Chronos &lt; <b>-0.85</b> &middot; BH &lt; 0.05\n"
                                  "positive control <b>POLR2A</b> must recover first")),
                   ("Data", "DepMap CRISPR (Chronos) + CN &times; TCGA GISTIC2.")]}]},
    {"id": "G3", "name": "MOA novelty &amp; orthogonality", "status": None,
     "asks": "Two questions: is the mechanism genuinely new relative to known payload classes, and &mdash; for "
             "the dual-payload program &mdash; is it orthogonal to exatecan's resistance?",
     "subs": [
         {"name": "G3a &middot; Embed &amp; locate the holes", "status": ("defer", "deferred &mdash; contract"),
          "rows": [("Asks", "How far is the candidate from the nearest known-payload class in mechanism space?"),
                   ("Method", "Joint embedding of LINCS L1000 + Tahoe pseudobulk + JUMP morphology; annotate known "
                    "payload regions from ADCdb; score distance to nearest centroid <b>conditional on passing G1</b>."),
                   ("Threshold", ("thr", "min distance to known centroid &ge; <b>0.35</b>")),
                   ("Data", "LINCS L1000 &middot; Tahoe &middot; JUMP-CP &middot; ADCdb.")]},
         {"name": "G3b &middot; Orthogonal-resistance search", "tag": "the exatecan partner", "status": ("ok", "implemented + report"),
          "rows": [("Asks", "Which compound is most potent on exactly the lines exatecan cannot handle?"),
                   ("Method", "Regress out the Top1i component (SLFN11 axis), rank by residual potency on the "
                    "Top1i-resistant (SLFN11-low) lines, stratified by ABCB1 so it is not just efflux escape."),
                   ("Threshold", ("thr", "SLFN11-low = lower <b>tertile</b> (0.333)\n"
                                  "partner's own SLFN11 slope &le; <b>0.1</b> &middot; top <b>50</b>")),
                   ("Data", "PRISM / CTRP sensitivity &times; DepMap SLFN11, ABCB1.")]}]},
    {"id": "G4", "name": "Bystander competence", "status": ("partial", "physchem window (partial)"),
     "asks": '<span class="q">Will the released payload diffuse into neighbouring antigen-negative cells?</span> '
             'A <b>design choice</b>, not a universal good &mdash; essential for a heterogeneous CRC tumour, a '
             'liability for a narrow-window target. <b>Tag, do not filter.</b>',
     "rows": [("Method", "Do not rebuild the predictor: use the Guo <em>et&nbsp;al.</em> 2024 graph-attention "
               "B-score. Layer on charge state at lysosomal pH 4.8 vs cytosol 7.2 (the MMAE/MMAF split) and "
               "cLogD(7.4), TPSA, MW via RDKit."),
              ("Threshold", ("thr", "B-score &ge; <b>1.5</b>  (Guo et al.)\n"
                             "cLogD <b>[-1, 3]</b> &middot; TPSA &le; <b>120</b> &middot; MW &le; <b>1000</b>")),
              ("Data", "RDKit over the G1-passing set &cap; COCONUT / NPAtlas / ChEMBL."),
              ("Caveat", ("caveat", "Currently the physchem window only, cLogP as a cLogD proxy &mdash; not the "
                          "trained B-score. Exatecan itself fails this window, so read G4 as <b>indicative</b>."))]},
    {"id": "G5", "name": "Conjugatability", "status": ("partial", "SMARTS done (partial)"),
     "asks": '<span class="q">Is there a conjugatable handle at a position SAR says tolerates substitution?</span> '
             'A handle in the pharmacophore is not a handle.',
     "rows": [("Method", "SMARTS match for a primary/secondary amine, hydroxyl, thiol or carboxylic acid, then "
               "cross-reference the scaffold's ChEMBL SAR series for a substitutable position."),
              ("Threshold", ("thr", "handle &isin; { 1&deg;/2&deg; amine, -OH, -SH, -COOH }\n"
                             "present AND at an SAR-tolerant position")),
              ("Data", "Structure (SMILES via PubChem) + ChEMBL SAR."),
              ("Caveat", ("caveat", "Handle detection is wired; the SAR-position check is pending, so a pass is "
                          "<b>necessary, not sufficient</b>."))]},
    {"id": "G6", "name": "Therapeutic window", "status": ("partial", "expression done (partial)"),
     "asks": '<span class="q">Is the payload\'s target enriched in the organs that recur as ADC dose-limiting '
             'toxicities?</span> Payload class sets the DLT organ far more than the antigen does.',
     "rows": [("Method", "(1) FAERS disproportionality (ROR / PRR, shrinkage) on ADC regimens, stratified by "
               "payload class &rarr; an empirical class&rarr;toxicity map. (2) Score the target's normal-tissue "
               "expression across the five recurring DLT compartments; low across all five is the window signal."),
              ("Threshold", ("thr", "target &le; <b>25th percentile</b> in ALL of:\n"
                             "  HSC/marrow &middot; GI crypt &middot; cornea\n"
                             "  alveolar type II &middot; peripheral nerve")),
              ("Data", "GTEx / HPA normal-tissue expression &middot; FAERS (ADC regimens)."),
              ("Caveat", ("caveat", "HPA expression implemented (retina / spinal-cord proxies for cornea / nerve); "
                          "scores <b>on-target risk only</b>. FAERS class&rarr;tox map still pending."))]},
]


def _cascade_rows(rows):
    cells = []
    for i, (label, val) in enumerate(rows):
        first = " first" if i == 0 else ""
        if isinstance(val, tuple) and val[0] == "thr":
            dd = f'<span class="thr">{val[1]}</span>'
        elif isinstance(val, tuple) and val[0] == "caveat":
            dd = f'<span class="caveat">{val[1]}</span>'
        else:
            dd = val
        cells.append(f'<dt class="l{first}">{label}</dt><dd class="l{first}">{dd}</dd>')
    return f'<dl class="grid">{"".join(cells)}</dl>'


def _cascade_card(g):
    h = ['<div class="gate"><div class="gate-head">',
         f'<span class="badge">{g["id"]}</span>',
         f'<h3 class="gate-name">{g["name"]}</h3>']
    if g.get("status"):
        h.append(f'<span class="pill {g["status"][0]}"><span class="dot"></span>{g["status"][1]}</span>')
    h.append(f'</div><p class="asks">{g["asks"]}</p>')
    if g.get("rows"):
        h.append(_cascade_rows(g["rows"]))
    for s in g.get("subs", []):
        h.append('<div class="g-sub"><div class="sub-name">' + s["name"])
        if s.get("tag"):
            h.append(f'<span class="sub-tag">{s["tag"]}</span>')
        if s.get("status"):
            h.append(f'<span class="pill {s["status"][0]}"><span class="dot"></span>{s["status"][1]}</span>')
        h.append('</div>' + _cascade_rows(s["rows"]) + '</div>')
    if g.get("callout"):
        h.append(f'<div class="callout"><div class="ct">&#9670; {g["callout"][0]}</div>{g["callout"][1]}</div>')
    h.append('</div>')
    return "".join(h)


def _cascade_div():
    rules = [("Conjunctive", "Gates are AND-ed. A fatal gate is fatal; nothing is averaged away."),
             ("Nothing dropped silently", 'Every gate returns <span class="mono">(passed, failed, reason)</span>.'),
             ("Free drug &ne; ADC", "Linker, DAR and internalisation move IC50 1&ndash;2 logs. Every output is a "
              "hypothesis <em>for conjugation</em>."),
             ("Identity is resolved", "Compound &rarr; InChIKey; cell line &rarr; DepMap ModelID, at ingest.")]
    parts = ['<div class="cascade">',
             '<p class="intro">An <b>explicit conjunctive filter</b> for prioritising novel ADC payload classes. '
             'Six gates, applied in order. A candidate that fails any gate is dead regardless of how it scores '
             'elsewhere &mdash; there is <b>no composite score</b>. Every threshold below is a stated, arguable '
             'number that lives in <span class="mono">config/gates.yaml</span>.</p>',
             '<div class="rules">']
    for k, v in rules:
        parts.append(f'<div class="rule"><div class="k">{k}</div><div class="v">{v}</div></div>')
    parts.append('</div><div class="legend">'
                 '<span class="pill ok"><span class="dot"></span>implemented + tested</span>'
                 '<span class="pill partial"><span class="dot"></span>partial / real-data pass</span>'
                 '<span class="pill defer"><span class="dot"></span>deferred &mdash; contract only</span></div>')
    for g in _CASCADE_GATES:
        parts.append(_cascade_card(g))
    parts.append('<p style="color:#8a8f98;font-size:12px;margin-top:22px">Source of truth: '
                 '<span class="mono">config/gates.yaml</span> (move a threshold there and re-run &mdash; '
                 'threshold sensitivity is a required deliverable). Full framing: '
                 '<span class="mono">WARHEAD.md</span>.</p></div>')
    return "".join(parts)


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
        elif sc["type"] == "conjugation":
            body.append('<h3>Conjugation-suitability ledger - PRISM top-20 (CRC)</h3>')
            body.append(
                '<p style="color:#444;font-size:13px;max-width:900px">Free-drug potency is not ADC '
                'potency. Each column is a separate, <b>conjunctive</b> question a payload must pass '
                '(no averaging): <b>G1</b> sub-nM potency + complete kill; <b>G2a</b> not an ABCB1/ABCG2 '
                'efflux substrate; <b>G2b</b> proliferation-independent; <b>G5</b> a conjugatable handle; '
                '<b>G4</b> bystander-permissive physchem; <b>G6</b> target not enriched in the DLT organs. '
                '<b>✓</b> = passes, <b>✗</b> = fails, n/a = not assessable.</p>')
            body.append(
                '<p style="color:#6E1426;font-size:13px;max-width:900px;font-weight:600">The compounds that '
                'clear the G1 potency bar are the validated payload chemotypes - exatecan-mesylate (0.1 nM, '
                'the deruxtecan warhead, sitting in the table as a reference anchor), maytansinol, '
                'dolastatin-10, plus gemcitabine and triptolide; most other hits are 1-1.5 logs short. All '
                'fail G6 - their targets are broadly-essential machinery expressed throughout the DLT organs. '
                'PRISM is used here (not CTRP) because its IC50 is calibrated so exatecan reads its true '
                'sub-nM potency; CTRP reads ~4.7x weaker and would fail every payload against an absolute bar.</p>')
            body.append(_scorecard_table(sc["scorecard"]))
            body.append(
                '<p style="color:#777;font-size:11.5px;max-width:900px">gap = log10 units the median IC50 '
                'sits above 1 nM (negative = below the bar). G2a: log10 IC50 vs ABCB1/ABCG2 expression '
                '(std slope ≥0.25 &amp; q&lt;0.05 = substrate). G2b: pan-panel Spearman vs DepMap growth '
                '(proxy; per-compound calls are noisy). G5 handle is necessary not SAR-verified; G4 is a '
                'cLogP-proxy window, not the Guo B-score (exatecan itself fails it, so read G4 as indicative); '
                'G6 uses HPA with retina/spinal-cord proxies for cornea/nerve, on-target risk only '
                '(FAERS class→tox pending).</p>')
        elif sc["type"] == "indication":
            ind = sc["indication"]; scrs = sc["screens"]
            body.append(f'<h3>{html.escape(ind)} potency &amp; selectivity - {" + ".join(scrs)}</h3>')
            if sc.get("caveat"):
                body.append(f'<p style="color:#444;font-size:13px;max-width:900px">{sc["caveat"]}</p>')
            body.append(_sel_pair_div([(f"{s} {ind}", sc["sel"][s]) for s in scrs], f"sel_{sid}"))
            for scr in scrs:
                body.append(f"<h3>Top 20 by EC90 - {scr} ({ind})</h3>")
                body.append(_table_cov(sc["rank"][scr].head(20), tested))
        elif sc["type"] == "cascade":
            body.append(_cascade_div())
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
    .searchbar{padding:0 26px 12px;display:flex;align-items:center;gap:10px}
    .searchbar input{width:min(460px,70vw);padding:8px 12px;font-size:13px;font-family:inherit;
      border:1px solid #d8c7cc;border-radius:8px;background:#fff;color:#17181d;outline:none}
    .searchbar input:focus{border-color:#6E1426;box-shadow:0 0 0 2px rgba(110,20,38,.12)}
    .searchbar .hits{font-size:12px;color:#888}
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
    /* --- Cascade (gate reference) tab --- */
    .cascade{max-width:940px}
    .cascade .intro{font-size:15px;color:#444;max-width:70ch;margin:2px 0 4px;line-height:1.55}
    .cascade .intro b{color:#17181d;font-weight:600}
    .cascade .rules{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:12px;margin:20px 0 8px}
    .cascade .rule{background:#fff;border:1px solid #eadfe2;border-radius:12px;padding:13px 15px}
    .cascade .rule .k{font-family:'Fraunces',serif;font-size:15px;font-weight:600;margin-bottom:3px;color:#17181d}
    .cascade .rule .v{font-size:12.5px;color:#6b6268;line-height:1.42}
    .cascade .legend{display:flex;flex-wrap:wrap;gap:16px;align-items:center;margin:22px 0 4px;font-size:12.5px;color:#6b6268}
    .cascade .pill{display:inline-flex;align-items:center;gap:6px;font-weight:600;font-size:11.5px;padding:3px 10px;border-radius:999px}
    .cascade .dot{width:8px;height:8px;border-radius:50%}
    .cascade .pill.ok{color:#2E7D6B;background:#e6f1ee} .cascade .pill.ok .dot{background:#2E7D6B}
    .cascade .pill.partial{color:#a06e1e;background:#f6ecd9} .cascade .pill.partial .dot{background:#B87C22}
    .cascade .pill.defer{color:#7d828b;background:#eef0f2} .cascade .pill.defer .dot{background:#8A8F98}
    .cascade .gate{background:#fff;border:1px solid #eadfe2;border-radius:16px;padding:22px 24px;margin:14px 0;box-shadow:0 1px 2px rgba(40,10,18,.04),0 6px 18px rgba(40,10,18,.045)}
    .cascade .gate-head{display:flex;align-items:baseline;gap:14px;flex-wrap:wrap}
    .cascade .badge{font-family:'Fraunces',serif;font-weight:600;font-size:14px;color:#fff;background:#6E1426;border-radius:8px;padding:5px 10px;line-height:1}
    .cascade .gate-name{font-family:'Fraunces',serif;font-weight:500;font-size:23px;margin:0;flex:1 1 auto;color:#17181d}
    .cascade .asks{font-size:15px;margin:12px 0 16px;color:#17181d;line-height:1.5}
    .cascade .asks .q{color:#6E1426;font-weight:600}
    .cascade dl.grid{display:grid;grid-template-columns:112px 1fr;gap:0 18px;font-size:14px;margin:0}
    .cascade dl.grid dt{color:#6b6268;font-size:11px;letter-spacing:.06em;text-transform:uppercase;font-weight:600;padding-top:10px;border-top:1px solid #f0e7ea}
    .cascade dl.grid dd{margin:0;padding-top:10px;border-top:1px solid #f0e7ea}
    .cascade dl.grid dt.first,.cascade dl.grid dd.first{border-top:0;padding-top:2px}
    .cascade .thr{font-family:'IBM Plex Mono',monospace;font-size:12.5px;background:#fbf6f7;border:1px solid #eadfe2;border-radius:7px;padding:8px 11px;display:block;line-height:1.7;white-space:pre-wrap;overflow-x:auto;margin:2px 0}
    .cascade .thr b{color:#6E1426;font-weight:500}
    .cascade .caveat{color:#6b6268;font-size:13px}
    .cascade .g-sub{border-left:2px solid #dcced2;padding-left:18px;margin:20px 0 4px}
    .cascade .g-sub + .g-sub{margin-top:24px}
    .cascade .sub-name{font-family:'Fraunces',serif;font-weight:600;font-size:17px;margin:0 0 6px;display:flex;align-items:center;gap:9px;flex-wrap:wrap;color:#17181d}
    .cascade .sub-tag{font-family:'IBM Plex Mono',monospace;font-size:11.5px;color:#6E1426;background:#f3e7ea;padding:2px 7px;border-radius:5px;font-weight:500}
    .cascade .callout{background:#FBEFC7;border:1px solid #E4C868;color:#5A4410;border-radius:12px;padding:14px 16px;margin:18px 0 2px;font-size:13.5px;line-height:1.55}
    .cascade .callout .ct{font-family:'Fraunces',serif;font-weight:600;font-size:14.5px;margin-bottom:5px}
    .cascade .callout code{font-family:'IBM Plex Mono',monospace;font-size:.92em;background:rgba(90,68,16,.12);padding:1px 5px;border-radius:4px}
    .cascade .mono{font-family:'IBM Plex Mono',monospace}
    """
    js = """
    function showTab(id){
      document.querySelectorAll('.tab-pane').forEach(p=>p.classList.remove('active'));
      document.querySelectorAll('.tab-btn').forEach(b=>b.classList.remove('active'));
      document.getElementById(id).classList.add('active');
      document.getElementById('btn_'+id).classList.add('active');
      document.querySelectorAll('#'+id+' .js-plotly-plot').forEach(d=>{try{Plotly.Plots.resize(d);}catch(e){}});
    }
    function highlightCompound(q){
      q=(q||'').trim().toLowerCase(); let hits=0;
      document.querySelectorAll('.js-plotly-plot').forEach(gd=>{
        if(!gd.data||!gd.layout) return;
        const anns=[];
        if(q.length>=2){
          gd.data.forEach(tr=>{
            if(!tr.customdata||!tr.x||!tr.y) return;
            for(let i=0;i<tr.customdata.length;i++){
              const cd=tr.customdata[i];
              if(cd&&String(cd[0]).toLowerCase().indexOf(q)>-1){
                hits++;
                anns.push({x:tr.x[i],y:tr.y[i],xref:tr.xaxis||'x',yref:tr.yaxis||'y',
                  text:'<b>'+cd[0]+'</b>',showarrow:true,arrowhead:2,arrowsize:1,arrowwidth:1.4,
                  arrowcolor:'#6E1426',ax:0,ay:-30,font:{color:'#6E1426',size:11},
                  bgcolor:'rgba(255,255,255,.9)',bordercolor:'#6E1426',borderwidth:1,borderpad:2});
              }
            }
          });
        }
        try{Plotly.relayout(gd,{annotations:anns});}catch(e){}
      });
      const h=document.getElementById('searchHits');
      if(h) h.textContent = q.length<2 ? '' : (hits? hits+' point'+(hits>1?'s':'')+' labelled' : 'no match');
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
           f"<link rel='preconnect' href='https://fonts.gstatic.com' crossorigin>"
           f"<link rel='stylesheet' href='https://fonts.googleapis.com/css2?"
           f"family=Fraunces:opsz,wght@9..144,500;9..144,600&family=IBM+Plex+Sans:wght@400;500;600&"
           f"family=IBM+Plex+Mono:wght@400;500&display=swap'>"
           f"<style>{css}</style><script>{get_plotlyjs()}</script></head><body>"
           f"<h1>WARHEAD - public drug-screen dashboard</h1>"
           f"<div class='sub'>EC90/IC50 potency, HCC/CRC selectivity and clinical/ADC context across five public screens. "
           f"Click a tab; hover points; click a table header to sort.</div>"
           f"<div class='searchbar'><input id='cmpdSearch' oninput='highlightCompound(this.value)' "
           f"placeholder='search a compound - label it across every selectivity plot (e.g. daporinad)' autocomplete='off'>"
           f"<span id='searchHits' class='hits'></span></div>"
           f"<div class='tabbar'>{''.join(tabs_btn)}</div>{''.join(tabs_div)}"
           f"<script>{js}\nshowTab('{first}');</script></body></html>")
    out_path.write_text(doc, encoding="utf-8")
    return out_path
