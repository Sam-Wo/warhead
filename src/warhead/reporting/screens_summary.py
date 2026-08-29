"""reports/screens_summary.pdf - cross-source view of the top CRC compounds.

Shows, for the union of the most potent compounds, the median EC90 in CRC in each
screen (so replication is visible), plus target, clinical status and whether the
chemotype is an ADC payload - the decision-relevant integration.
"""
from __future__ import annotations

import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

RRB_MAROON = "#6E1426"
_PHASE_COLOR = {"Launched": "#6E1426", "Approved": "#6E1426", "Phase 3": "#9E3A50",
                "Phase 2": "#C06A7C", "Phase 1": "#D9A6B0", "Preclinical": "#9AA0A6"}


def _norm(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(name).lower())


def build_summary(source_ranks: dict, clin_tox: pd.DataFrame, *, indication="CRC", top=15) -> pd.DataFrame:
    """source_ranks: {source_name: rank_potency_frame (has compound, target,
    median_ec90_nM, median_emax, clinical_phase)}. Returns one row per compound in
    the union of each source's top-N, with an EC90 column per source."""
    union = {}
    for src, r in source_ranks.items():
        for _, row in r.head(top).iterrows():
            union.setdefault(_norm(row["compound"]),
                             {"compound": row["compound"], "target": row.get("target")})
    # curated compounds -> their meaningful name tokens (>=5 chars) for matching
    def _tokens(name):
        return [t for t in re.split(r"[^a-z0-9]+", str(name).lower()) if len(t) >= 5]
    tox = clin_tox.assign(_tok=clin_tox["compound"].map(_tokens))
    rows = []
    for k, base in union.items():
        rec = dict(base)
        for src, r in source_ranks.items():
            m = r[r["compound"].map(_norm) == k]
            rec[f"ec90_{src}"] = float(m["median_ec90_nM"].iloc[0]) if len(m) else np.nan
            rec[f"ic50_{src}"] = float(m["median_ic50_nM"].iloc[0]) if len(m) else np.nan
            if len(m) and pd.isna(rec.get("clinical_phase")) and pd.notna(m["clinical_phase"].iloc[0]):
                rec["clinical_phase"] = m["clinical_phase"].iloc[0]
            if len(m) and (pd.isna(rec.get("emax")) if "emax" in rec else True) and pd.notna(m["median_emax"].iloc[0]):
                rec["emax"] = m["median_emax"].iloc[0]
        # curated clinical/payload: match if any curated token is inside the compound
        hit = tox[tox["_tok"].apply(lambda toks: any(t in k for t in toks))]
        rec["adc_payload"] = hit["adc_payload_status"].iloc[0] if len(hit) else ""
        rec["dlt"] = hit["patient_toxicity_DLT"].iloc[0] if len(hit) else ""
        if (not rec.get("clinical_phase") or pd.isna(rec.get("clinical_phase"))) and len(hit):
            rec["clinical_phase"] = hit["clinical_status"].iloc[0]
        rows.append(rec)
    df = pd.DataFrame(rows)
    src0 = f"ec90_{list(source_ranks)[0]}"
    sort_col = "ec90_PRISM Repurposing (secondary)" if "ec90_PRISM Repurposing (secondary)" in df else src0
    return df.sort_values(sort_col, na_position="last").reset_index(drop=True)


def _fmt(v):
    if pd.isna(v):
        return "-"
    if v >= 1000:
        return f"{v/1000:.1f}k" if v < 1e6 else ">1M"
    return f"{v:.0f}" if v >= 10 else f"{v:.1f}"


def render_summary_heatmap(df: pd.DataFrame, source_names: list, meta: pd.DataFrame | None = None,
                           *, out_path, indication="CRC", tested: dict | None = None) -> Path:
    """Heatmap table: per source an IC50 and EC90 cell (coloured by potency, value
    printed), plus target, Emax, clinical status and ADC-payload status. A screen
    metadata block sits on top.

    tested: optional {source_name: set(normalised compound names ever assayed)}. When
    given, an empty cell is disambiguated - a hatched grey "n/t" means the compound was
    never in that screen's library, versus a faint dot for tested-but-outside-the-gate.
    """
    from matplotlib.colors import LogNorm
    out_path = Path(out_path); out_path.parent.mkdir(parents=True, exist_ok=True)
    n = len(df)
    short = {s: s.split()[0] for s in source_names}
    meta_h = 0 if meta is None else len(meta) + 2

    fig_h = 0.34 * n + 0.24 * meta_h + 2.2
    fig = plt.figure(figsize=(15.5, fig_h))
    ax = fig.add_axes([0.008, 0.015, 0.984, 0.955]); ax.axis("off")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    fig.suptitle(f"WARHEAD - cross-source EC90 / IC50 of the top {indication} compounds",
                 color=RRB_MAROON, fontsize=15, fontweight="bold", y=.997)

    # ---- screen metadata block (top) ----
    top_y = 0.965
    if meta is not None:
        ax.text(0.0, top_y, "SCREEN METADATA", fontsize=8.5, fontweight="bold", color=RRB_MAROON)
        mcols = [("source", 0.0, 15), ("compounds", 0.16, 8), ("cell_lines", 0.225, 8),
                 ("CRC_lines", 0.29, 6), ("HCC_lines", 0.345, 6), ("dose_range", 0.40, 34),
                 ("n_doses", 0.60, 6), ("metrics", 0.66, 42)]
        hdr = {"source": "screen", "compounds": "cpds", "cell_lines": "lines", "CRC_lines": "CRC",
               "HCC_lines": "HCC", "dose_range": "dose range", "n_doses": "#dose", "metrics": "metrics"}
        for key, x, _ in mcols:
            ax.text(x, top_y - 0.018, hdr[key], fontsize=7, fontweight="bold", color="#444")
        for i, (_, r) in enumerate(meta.iterrows()):
            yy = top_y - 0.033 - i * 0.017
            for key, x, w in mcols:
                ax.text(x, yy, str(r[key])[:w], fontsize=6.6, color="#333",
                        family="monospace" if key in ("dose_range", "metrics", "n_doses") else "sans-serif")
        table_top = top_y - 0.033 - len(meta) * 0.017 - 0.03
    else:
        table_top = top_y
    ax.text(0.0, table_top + 0.012,
            f"COMPOUNDS  -  cell = median in {indication} (nM), coloured by potency (dark = potent); value printed.",
            fontsize=8, fontweight="bold", color=RRB_MAROON)

    # ---- column geometry for the compound table ----
    x_comp, x_tgt = 0.0, 0.15
    metric_x = {}
    x = 0.32
    for s in source_names:
        metric_x[(s, "ic50")] = x; metric_x[(s, "ec90")] = x + 0.052
        x += 0.125
    x_emax = x + 0.005; x_clin = x_emax + 0.055; x_pay = x_clin + 0.10
    cw = 0.024  # half cell width

    row_h = (table_top - 0.06) / max(n, 1)
    hy = table_top - 0.005
    ax.text(x_comp, hy, "compound", fontsize=8, fontweight="bold")
    ax.text(x_tgt, hy, "target", fontsize=8, fontweight="bold")
    for s in source_names:
        cx = (metric_x[(s, "ic50")] + metric_x[(s, "ec90")]) / 2
        ax.text(cx, hy + 0.008, short[s], fontsize=7.5, fontweight="bold", ha="center", color=RRB_MAROON)
        ax.text(metric_x[(s, "ic50")], hy - 0.006, "IC50", fontsize=6.3, ha="center", color="#666")
        ax.text(metric_x[(s, "ec90")], hy - 0.006, "EC90", fontsize=6.3, ha="center", color="#666")
    ax.text(x_emax, hy, "Emax", fontsize=7.5, fontweight="bold", ha="center")
    ax.text(x_clin, hy, "clinical", fontsize=7.5, fontweight="bold")
    ax.text(x_pay, hy, "ADC payload?", fontsize=7.5, fontweight="bold")

    norm = LogNorm(vmin=1, vmax=1e4); cmap = plt.cm.RdPu_r
    for i, row in df.reset_index(drop=True).iterrows():
        yc = table_top - 0.05 - (i + 0.5) * row_h
        ax.text(x_comp, yc, str(row["compound"])[:22], fontsize=7.6, va="center", fontweight="bold")
        ax.text(x_tgt, yc, str(row.get("target") or "")[:24], fontsize=6.3, va="center",
                color="#555", family="monospace")
        ckey = _norm(row["compound"])
        for s in source_names:
            not_tested = tested is not None and ckey not in tested.get(s, set())
            for metric in ("ic50", "ec90"):
                v = row.get(f"{metric}_{s}")
                mx = metric_x[(s, metric)]
                if pd.notna(v):
                    col = cmap(norm(max(v, 1)))
                    ax.add_patch(plt.Rectangle((mx - cw, yc - row_h * .42), 2 * cw, row_h * .84,
                                               facecolor=col, edgecolor="#eee", linewidth=.3, zorder=1))
                    lum = 0.299 * col[0] + 0.587 * col[1] + 0.114 * col[2]
                    ax.text(mx, yc, _fmt(v), fontsize=6.0, ha="center", va="center",
                            color="white" if lum < 0.5 else "#222", zorder=2)
                elif not_tested:
                    ax.add_patch(plt.Rectangle((mx - cw, yc - row_h * .42), 2 * cw, row_h * .84,
                                               facecolor="#EDEDED", edgecolor="#ddd", linewidth=.3,
                                               hatch="////", zorder=1))
                    ax.text(mx, yc, "n/t", fontsize=5.3, ha="center", va="center", color="#999", zorder=2)
                else:
                    ax.text(mx, yc, "·", fontsize=8, ha="center", va="center", color="#c7c7c7")
        em = row.get("emax")
        ax.text(x_emax, yc, f"{em:.2f}" if pd.notna(em) else "-", fontsize=6.8, ha="center", va="center")
        ph = str(row.get("clinical_phase") or "-")
        cc = "#6E1426" if ph.startswith(("Launched", "Approved", "FDA")) else ("#B05468" if ph.startswith("Phase") else "#888")
        ax.text(x_clin, yc, ph[:16], fontsize=6.8, va="center", color=cc, fontweight="bold")
        pay = str(row.get("adc_payload") or "")
        is_pay = "payload" in pay.lower() and "not a payload" not in pay.lower()
        ax.text(x_pay, yc, ("● " if is_pay else "○ ") + pay[:30], fontsize=6.2, va="center",
                color=RRB_MAROON if is_pay else "#999")

    legend = ("IC50 = fitted 50% viability; EC90 = 90% of max effect. GDSC EC90 is "
              "extrapolation-inflated (narrow windows, bottom=0); PRISM/CTRP EC90 use a real Emax.")
    if tested is not None:
        legend += ("   Cells:  coloured = potent hit (value in nM);  " + r"$\cdot$"
                   " = tested but outside the potency gate;  hatched n/t = not in that screen's library.")
    fig.text(0.008, 0.006, legend, fontsize=7, color="#777")
    fig.savefig(out_path, format=out_path.suffix.lstrip(".").lower() or "pdf", dpi=150)
    plt.close(fig)
    return out_path


def render_summary(df: pd.DataFrame, source_names: list, *, out_path, indication="CRC") -> Path:
    out_path = Path(out_path); out_path.parent.mkdir(parents=True, exist_ok=True)
    n = len(df)
    fig, ax = plt.subplots(figsize=(13, 0.42 * n + 2.2))
    ax.set_xlim(0, 12); ax.set_ylim(-0.5, n - 0.5); ax.invert_yaxis(); ax.axis("off")
    fig.suptitle(f"WARHEAD - cross-source EC90 of the top {indication} compounds",
                 color=RRB_MAROON, fontsize=14, fontweight="bold", y=.99)
    ax.text(0, -1.15, "per source: dot = median EC90 in " + indication + " (nM; darker = more potent), "
            "small grey = median IC50 (nM).  Emax = residual viability (PRISM/CTRP).",
            fontsize=8.5, color="#555")

    # column x positions
    x_comp, x_tgt = 0.0, 3.2
    src_x = {s: 6.0 + i * 0.9 for i, s in enumerate(source_names)}
    x_emax, x_clin, x_payload = 6.0 + len(source_names) * 0.9 + 0.2, None, None
    x_clin = x_emax + 0.7
    x_payload = x_clin + 1.4
    norm = LogNorm(vmin=1, vmax=3000)
    cmap = plt.cm.RdPu_r

    # headers
    ax.text(x_comp, -0.5, "compound", fontsize=8, fontweight="bold")
    ax.text(x_tgt, -0.5, "target", fontsize=8, fontweight="bold")
    for s in source_names:
        ax.text(src_x[s], -0.5, s.split()[0] + "\nEC90·IC50", fontsize=7, fontweight="bold", ha="center", va="top")
    ax.text(x_emax, -0.5, "Emax", fontsize=8, fontweight="bold", ha="center")
    ax.text(x_clin, -0.5, "clinical", fontsize=8, fontweight="bold")
    ax.text(x_payload, -0.5, "ADC payload?", fontsize=8, fontweight="bold")

    for i, row in df.iterrows():
        ax.text(x_comp, i, str(row["compound"])[:22], fontsize=8, va="center")
        ax.text(x_tgt, i, str(row.get("target") or "")[:26], fontsize=6.6, va="center",
                color="#555", family="monospace")
        for s in source_names:
            v = row.get(f"ec90_{s}"); ic = row.get(f"ic50_{s}")
            if pd.notna(v):
                ax.scatter(src_x[s], i - 0.12, s=150, c=[cmap(norm(v))], edgecolor="#333", linewidth=.4, zorder=3)
                ax.text(src_x[s], i - 0.12, f"{v:.0f}", fontsize=5.6, ha="center", va="center",
                        color="white" if v < 300 else "#222", zorder=4)
                ax.text(src_x[s], i + 0.30, f"{ic:.0f}" if pd.notna(ic) else "-", fontsize=5.2,
                        ha="center", va="center", color="#888")
            else:
                ax.text(src_x[s], i, "-", fontsize=8, ha="center", va="center", color="#ccc")
        em = row.get("emax")
        ax.text(x_emax, i, f"{em:.2f}" if pd.notna(em) else "-", fontsize=7, ha="center", va="center")
        ph = str(row.get("clinical_phase") or "-")
        ax.text(x_clin, i, ph[:16], fontsize=7, va="center",
                color=_PHASE_COLOR.get(ph.split("/")[0], "#666"), fontweight="bold")
        pay = str(row.get("adc_payload") or "")
        is_pay = "payload" in pay.lower() and "not a payload" not in pay.lower()
        ax.text(x_payload, i, ("● " if is_pay else "○ ") + pay[:34], fontsize=6.4, va="center",
                color=RRB_MAROON if is_pay else "#999")

    fig.tight_layout(rect=[0, 0, 1, .97])
    fig.savefig(out_path, format=out_path.suffix.lstrip(".").lower() or "pdf", dpi=150)
    plt.close(fig)
    return out_path
