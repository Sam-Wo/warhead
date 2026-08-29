"""reports/clinical_toxicity.pdf - curated clinical validation + patient toxicity
+ ADC-payload status for the recurring top screen compounds (feeds G6)."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

RRB_MAROON = "#6E1426"


def _status_color(v: str) -> str:
    v = str(v)
    if v.startswith("Approved"):
        return "#6E1426"
    if v.startswith("Phase"):
        return "#B05468"
    return "#9AA0A6"


_COLS = [("compound", "compound", .00), ("target_class", "target / class", .17),
         ("clinical_status", "clinical status", .33), ("patient_toxicity_DLT", "characteristic patient toxicity (DLT)", .47),
         ("dlt_organ", "DLT organ (G6)", .74), ("adc_payload_status", "ADC payload?", .855)]
_MAXC = {"compound": 24, "target_class": 26, "clinical_status": 22,
         "patient_toxicity_DLT": 42, "dlt_organ": 22, "adc_payload_status": 24}


def render_clinical_tox(table, *, out_path):
    out_path = Path(out_path); out_path.parent.mkdir(parents=True, exist_ok=True)
    n = len(table)
    fig = plt.figure(figsize=(15.5, .40 * n + 1.6))
    ax = fig.add_axes([0.012, 0.03, 0.976, 0.9]); ax.axis("off")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    fig.suptitle("WARHEAD - clinical validation & patient toxicity of the recurring top screen compounds",
                 color=RRB_MAROON, fontsize=13, fontweight="bold", y=.99)

    for key, label, x in _COLS:
        ax.text(x, 1.0, label, fontsize=8.5, fontweight="bold", va="top")
    ax.axhline(0.975, color="#999", lw=.8)
    row_h = 0.95 / n
    for i, row in table.reset_index(drop=True).iterrows():
        yc = 0.955 - (i + 0.5) * row_h
        if i % 2 == 1:
            ax.axhspan(yc - row_h / 2, yc + row_h / 2, color="#F4F1F2", zorder=0)
        pay = str(row["adc_payload_status"])
        is_pay = "payload" in pay.lower() and "not a payload" not in pay.lower()
        for key, _, x in _COLS:
            val = str(row[key])[: _MAXC[key]]
            color, weight = "#222", "normal"
            mono = key in ("target_class", "dlt_organ")
            if key == "compound":
                weight = "bold"
            elif key == "clinical_status":
                color, weight = _status_color(val), "bold"
            elif key == "adc_payload_status":
                color = RRB_MAROON if is_pay else "#999"
                val = ("● " if is_pay else "○ ") + val
            ax.text(x, yc, val, fontsize=7.2, va="center", color=color, fontweight=weight,
                    family="monospace" if mono else "sans-serif")
    fig.text(0.012, 0.006, "curated from FDA labels / standard oncology references / ADC literature; "
             "characteristic clinical pattern, not a per-patient prediction.", fontsize=7, color="#777")
    fig.savefig(out_path, format=out_path.suffix.lstrip(".").lower() or "pdf", dpi=150)
    plt.close(fig)
    return out_path
