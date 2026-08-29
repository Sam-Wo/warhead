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
            if len(m) and pd.isna(rec.get("clinical_phase")) and pd.notna(m["clinical_phase"].iloc[0]):
                rec["clinical_phase"] = m["clinical_phase"].iloc[0]
            if len(m) and "emax" not in rec:
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


def render_summary(df: pd.DataFrame, source_names: list, *, out_path, indication="CRC") -> Path:
    out_path = Path(out_path); out_path.parent.mkdir(parents=True, exist_ok=True)
    n = len(df)
    fig, ax = plt.subplots(figsize=(13, 0.42 * n + 2.2))
    ax.set_xlim(0, 12); ax.set_ylim(-0.5, n - 0.5); ax.invert_yaxis(); ax.axis("off")
    fig.suptitle(f"WARHEAD - cross-source EC90 of the top {indication} compounds",
                 color=RRB_MAROON, fontsize=14, fontweight="bold", y=.99)
    ax.text(0, -1.15, "dot = median EC90 in " + indication + " (nM; darker/left = more potent).  "
            "clin = clinical status.  Emax = residual viability (PRISM).", fontsize=8.5, color="#555")

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
        ax.text(src_x[s], -0.5, s.split()[0], fontsize=7.5, fontweight="bold", ha="center")
    ax.text(x_emax, -0.5, "Emax", fontsize=8, fontweight="bold", ha="center")
    ax.text(x_clin, -0.5, "clinical", fontsize=8, fontweight="bold")
    ax.text(x_payload, -0.5, "ADC payload?", fontsize=8, fontweight="bold")

    for i, row in df.iterrows():
        ax.text(x_comp, i, str(row["compound"])[:22], fontsize=8, va="center")
        ax.text(x_tgt, i, str(row.get("target") or "")[:26], fontsize=6.6, va="center",
                color="#555", family="monospace")
        for s in source_names:
            v = row.get(f"ec90_{s}")
            if pd.notna(v):
                ax.scatter(src_x[s], i, s=140, c=[cmap(norm(v))], edgecolor="#333", linewidth=.4, zorder=3)
                ax.text(src_x[s], i, f"{v:.0f}", fontsize=5.6, ha="center", va="center",
                        color="white" if v < 300 else "#222", zorder=4)
            else:
                ax.text(src_x[s], i, "-", fontsize=8, ha="center", va="center", color="#bbb")
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
