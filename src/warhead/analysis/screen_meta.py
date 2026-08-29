"""Per-screen metadata (assay design + coverage) for the summary and dashboard.

Counts (compounds / lines / CRC / HCC) are computed from the canonical frames
where available; the assay design fields (dose range, #doses, metrics) are curated
from each screen's documentation.
"""
from __future__ import annotations

import pandas as pd

# curated assay-design fields per source
_DESIGN = {
    "GDSC2": dict(dose_range="drug-specific ~1000-fold (median max 10 uM)", n_doses="5-9",
                  metrics="IC50, AUC (bottom fixed at 0 - no Emax)", assay="CellTiter-Glo, 72 h"),
    "PRISM Repurposing (secondary)": dict(dose_range="8-pt, up to ~10 uM", n_doses="8",
                  metrics="IC50, EC50, slope, Emax, clinical phase", assay="PRISM barcoded, 5 d"),
    "CTRP v2": dict(dose_range="16-pt, up to ~66 uM", n_doses="16",
                  metrics="IC50, EC50, Hill, Emax, FDA status", assay="CellTiter-Glo, 72 h"),
    "PDXE (Novartis)": dict(dose_range="in-vivo (dosed to tolerability)", n_doses="-",
                  metrics="BestAvgResponse, ORR (RECIST-like)", assay="1x1x1 PDX, tumour volume"),
    "NCI-60": dict(dose_range="5-pt, 10 nM - 100 uM", n_doses="5",
                  metrics="GI50/TGI/LC50 (public data z-scored)", assay="SRB, 48 h"),
}


def screen_metadata(counts_by_source: dict | None = None) -> pd.DataFrame:
    """counts_by_source: {source: {'compounds','lines','CRC','HCC'}}; missing sources
    fall back to the curated defaults below."""
    default_counts = {
        "GDSC2": dict(compounds=286, lines=969, CRC=46, HCC=15),
        "PRISM Repurposing (secondary)": dict(compounds=1448, lines=480, CRC=27, HCC=16),
        "CTRP v2": dict(compounds=544, lines=887, CRC=50, HCC=22),
        "PDXE (Novartis)": dict(compounds=63, lines=275, CRC=59, HCC=0),
        "NCI-60": dict(compounds=25515, lines=60, CRC=7, HCC=0),
    }
    counts_by_source = counts_by_source or {}
    rows = []
    for src, design in _DESIGN.items():
        c = {**default_counts[src], **counts_by_source.get(src, {})}
        rows.append({"source": src, "compounds": c["compounds"], "cell_lines": c["lines"],
                     "CRC_lines": c["CRC"], "HCC_lines": c["HCC"], **design})
    return pd.DataFrame(rows)
