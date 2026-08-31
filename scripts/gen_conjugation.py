"""Conjugation-suitability scorecard for the top hits. Writes
reports/conjugation_scorecard_{IND}.pdf + reports/conjugation_scorecard.xlsx.

    PYTHONPATH=src py scripts/gen_conjugation.py

Primary screen is PRISM Repurposing (secondary). PRISM is chosen over CTRP for the
G1 POTENCY bar because its absolute IC50 is calibrated so the reference ADC payloads
sit at their true potency (exatecan 0.10 nM / 82% of lines sub-nM, SN-38 0.93 nM,
maytansinol 0.99 nM) - CTRP reads a systematic ~4.7x weaker (SN-38 = 121 nM there),
which would fail every payload against an absolute sub-nM bar. PRISM also carries a
real Emax (lower_limit) and clinical phase. CTRP remains the source for the
wide-window measured curves elsewhere.
"""
import pandas as pd

from warhead.analysis.clinical_tox import clinical_tox_table
from warhead.analysis.conjugation import (add_chemistry, add_efflux, add_window,
                                          delivery_scorecard, growth_lookup)
from warhead.analysis.screen_potency import rank_potency
from warhead.io.pubchem import fetch_smiles
from warhead.reporting.conjugation_report import render_scorecard

SCREEN = "PRISM (secondary)"
prism = pd.read_pickle("data/interim/prism_canonical.pkl")
clin = clinical_tox_table()
growth = growth_lookup()

sheets = {}
for ind in ("CRC", "HCC"):
    rk = rank_potency(prism, ind, emax_max=0.5)
    sc = delivery_scorecard(prism, rk, clin, indication=ind, growth_fn=growth, top=20)
    smiles = fetch_smiles(sc["compound"].tolist())
    sc = add_chemistry(sc, smiles)            # G5 handle + G4 bystander physchem
    sc = add_efflux(sc, prism)                # G2a efflux (ABCB1/ABCG2)
    sc = add_window(sc)                        # G6 therapeutic window (HPA DLT organs)
    out = render_scorecard(sc, out_path=f"reports/conjugation_scorecard_{ind}.pdf",
                           indication=ind, screen=SCREEN)
    print("wrote", out, f"| G1 pass {sc.g1_potency_pass.sum()}/{len(sc)} | "
          f"G2b indep {sc.g2b_independent.sum()}/{len(sc)} | "
          f"G2a substrate {int(sc.g2a_substrate.fillna(False).sum())}/{len(sc)} | "
          f"G5 handle {int(sc.g5_handle.fillna(False).sum())}/{len(sc)} | "
          f"G4 bystd {int(sc.g4_bystander.fillna(False).sum())}/{len(sc)} | "
          f"G6 win {int((sc.g6_window_ok==True).sum())}/{len(sc)}")
    sheets[ind] = sc

with pd.ExcelWriter("reports/conjugation_scorecard.xlsx", engine="openpyxl") as xw:
    for ind, sc in sheets.items():
        sc.to_excel(xw, sheet_name=ind, index=False)
print("wrote reports/conjugation_scorecard.xlsx")
