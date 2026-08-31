"""AML potency + selectivity, mirroring the CRC/HCC analysis, on the screens that
cover AML: CTRP v2 (30 AML lines) and GDSC2 (26 LAML lines). PRISM's secondary
subset has no blood-cancer lines, so it is skipped (as NCI-60 was for HCC).

    PYTHONPATH=src py scripts/gen_aml.py
-> reports/aml_selectivity.pdf + reports/aml_selectivity.xlsx
"""
import numpy as np
import pandas as pd

from warhead.analysis.indications import disease_line_set, label_indication
from warhead.analysis.screen_potency import rank_potency, selectivity
from warhead.reporting.indication_report import render_indication_report

IND = "AML"

# CTRP: relabel AML lines (cell-line names bridged to DepMap ModelIDs)
ctrp = pd.read_pickle("data/interim/ctrp_canonical.pkl")
ctrp = label_indication(ctrp, disease_line_set(ctrp["model_id"].unique(), "Acute Myeloid Leukemia"), IND)

# GDSC: per-line frame keyed by TCGA label; LAML = AML
g = pd.read_pickle("data/interim/gdsc2_ec90.pkl")
gmap = {"COREAD": "CRC", "LIHC": "HCC", "LAML": "AML"}
gdsc = pd.DataFrame({"source": "GDSC2", "compound": g.drug_name, "target": g.target, "moa": g.pathway,
                     "model_id": g.cell_line, "indication": g.tcga_desc.map(gmap).fillna("other"),
                     "ic50_nM": g.ic50_uM * 1e3, "ec90_nM": g.ec90_uM * 1e3, "emax": np.nan,
                     "ec90_extrapolated": g.ec90_range.eq("extrapolated"), "clinical_phase": pd.NA})

panels, sheets = [], {}
for label, frame in [("CTRP v2", ctrp), ("GDSC2", gdsc)]:
    rk = rank_potency(frame, IND, emax_max=0.5)
    sel = selectivity(frame, IND)
    n_lines = frame[frame.indication == IND].model_id.nunique()
    panels.append({"label": label, "rank": rk, "sel": sel,
                   "n_lines": n_lines, "n_total": frame.model_id.nunique()})
    sheets[f"{label} rank"] = rk
    sheets[f"{label} selectivity"] = sel
    n_sel = int(sel["selective"].sum()) if len(sel) else 0
    print(f"{label}: AML lines {n_lines}  |  ranked compounds {len(rk)}  |  "
          f"AML-selective (q<0.1) {n_sel}  |  selective+potent {int(sel['selective_potent'].sum()) if len(sel) else 0}")

CAVEAT = ("CAUTION: leukaemia lines are globally hypersensitive in vitro (median Delta-potency +0.2 "
          "to +0.5 logs), so the significance flag fires for ~90% of compounds - AML-vs-rest selectivity "
          "is confounded by a tissue-level sensitivity baseline. Read the TOP of each panel (largest "
          "Delta): those surface the real AML dependencies (Aurora B / barasertib, MCL1 / AZD5991, PLK1, "
          "BRD4, KIF11, cytarabine). A baseline-detrended version would give a cleaner selectivity call.")
out = render_indication_report(IND, panels, out_path="reports/aml_selectivity.pdf", caveat=CAVEAT)
print("wrote", out)
with pd.ExcelWriter("reports/aml_selectivity.xlsx", engine="openpyxl") as xw:
    for name, df in sheets.items():
        df.to_excel(xw, sheet_name=name[:31], index=False)
print("wrote reports/aml_selectivity.xlsx")
