"""Gastric (stomach-adenocarcinoma) potency + selectivity, same analysis as CRC/HCC,
on the three screens with gastric lines: CTRP v2 (30), GDSC2 (24 STAD), PRISM (17).
Gastric is a solid tumour, so - unlike AML - there is no global-sensitivity confound.

    PYTHONPATH=src py scripts/gen_gastric.py
-> reports/gastric_selectivity.pdf + reports/gastric_selectivity.xlsx
"""
import numpy as np
import pandas as pd

from warhead.analysis.indications import GASTRIC_CODES, label_indication, line_set
from warhead.analysis.screen_potency import rank_potency, selectivity
from warhead.reporting.indication_report import render_indication_report

IND = "Gastric"

ctrp = pd.read_pickle("data/interim/ctrp_canonical.pkl")
ctrp = label_indication(ctrp, line_set(ctrp["model_id"].unique(), codes=GASTRIC_CODES), IND)
prism = pd.read_pickle("data/interim/prism_canonical.pkl")
prism = label_indication(prism, line_set(prism["model_id"].unique(), codes=GASTRIC_CODES), IND)

g = pd.read_pickle("data/interim/gdsc2_ec90.pkl")
gmap = {"COREAD": "CRC", "LIHC": "HCC", "STAD": "Gastric"}
gdsc = pd.DataFrame({"source": "GDSC2", "compound": g.drug_name, "target": g.target, "moa": g.pathway,
                     "model_id": g.cell_line, "indication": g.tcga_desc.map(gmap).fillna("other"),
                     "ic50_nM": g.ic50_uM * 1e3, "ec90_nM": g.ec90_uM * 1e3, "emax": np.nan,
                     "ec90_extrapolated": g.ec90_range.eq("extrapolated"), "clinical_phase": pd.NA})

panels, sheets = [], {}
for label, frame in [("CTRP v2", ctrp), ("GDSC2", gdsc), ("PRISM", prism)]:
    rk = rank_potency(frame, IND, emax_max=0.5)
    sel = selectivity(frame, IND)
    n_lines = frame[frame.indication == IND].model_id.nunique()
    panels.append({"label": label, "rank": rk, "sel": sel,
                   "n_lines": n_lines, "n_total": frame.model_id.nunique()})
    sheets[f"{label} rank"] = rk
    sheets[f"{label} selectivity"] = sel
    n_sel = int(sel["selective"].sum()) if len(sel) else 0
    print(f"{label}: gastric lines {n_lines}  |  ranked {len(rk)}  |  "
          f"selective (q<0.1) {n_sel}  |  selective+potent {int(sel['selective_potent'].sum()) if len(sel) else 0}")

CAVEAT = ("Gastric is a solid tumour with no global-sensitivity confound; unlike CRC (MEK) or AML "
          "(Aurora/MCL1), it shows essentially NO potent+selective mechanism in these screens - it joins "
          "HCC as a lineage without a clean selective vulnerability here. Ranking (top) is still the "
          "most-potent-in-gastric compounds.")
out = render_indication_report(IND, panels, out_path="reports/gastric_selectivity.pdf", caveat=CAVEAT)
print("wrote", out)
with pd.ExcelWriter("reports/gastric_selectivity.xlsx", engine="openpyxl") as xw:
    for name, df in sheets.items():
        df.to_excel(xw, sheet_name=name[:31], index=False)
print("wrote reports/gastric_selectivity.xlsx")
