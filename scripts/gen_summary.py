"""Build the cross-screen summary heatmaps + workbook for BOTH indications:
reports/screens_summary_{CRC,HCC}.pdf (numbers-in-coloured-cells heatmap with a
screen-metadata block, per-source IC50 + EC90, Emax, target, clinical/ADC status)
and reports/screens_summary.xlsx.

    PYTHONPATH=src py scripts/gen_summary.py
"""
from warhead.analysis.clinical_tox import clinical_tox_table
from warhead.analysis.screen_potency import rank_potency
from warhead.reporting.screen_data import load_dr_screens
from warhead.reporting.screen_overlap import render_overlap_venn, tested_sets
from warhead.reporting.screens_summary import build_summary, render_summary_heatmap
import pandas as pd

cans, counts, meta = load_dr_screens()
clin = clinical_tox_table()
tested = tested_sets(cans)                 # {source: set of normalised names ever assayed}

# compound-library overlap across the three dose-response screens
print("wrote", render_overlap_venn(cans, out_path="reports/screens_overlap.pdf"))

for ind in ("CRC", "HCC"):
    sr = {s: rank_potency(c, ind, emax_max=0.5) for s, c in cans.items()}
    summ = build_summary(sr, clin, indication=ind, top=20)
    out = render_summary_heatmap(summ, list(cans), meta, tested=tested,
                                 out_path=f"reports/screens_summary_{ind}.pdf", indication=ind)
    print("wrote", out)

with pd.ExcelWriter("reports/screens_summary.xlsx", engine="openpyxl") as xw:
    meta.to_excel(xw, sheet_name="screen_metadata", index=False)
    for ind in ("CRC", "HCC"):
        sr = {s: rank_potency(c, ind, emax_max=0.5) for s, c in cans.items()}
        build_summary(sr, clin, indication=ind, top=20).to_excel(xw, sheet_name=ind, index=False)
print("wrote reports/screens_summary.xlsx")
