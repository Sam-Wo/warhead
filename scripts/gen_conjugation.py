"""Phase-A conjugation-suitability scorecard (delivery survivability) for the top
hits, from data already on disk. Writes reports/conjugation_scorecard_{IND}.pdf +
reports/conjugation_scorecard.xlsx.

    PYTHONPATH=src py scripts/gen_conjugation.py

Uses CTRP v2 as the primary screen (widest window, real Emax, best line coverage).
G2a efflux and the chemical gates G4/G5 are Phase B (need DepMap expression /
per-compound SMILES + RDKit) and appear as pending columns.
"""
import pandas as pd

from warhead.analysis.clinical_tox import clinical_tox_table
from warhead.analysis.conjugation import (add_chemistry, add_efflux, add_window,
                                          delivery_scorecard, growth_lookup)
from warhead.analysis.screen_potency import rank_potency
from warhead.io.pubchem import fetch_smiles
from warhead.reporting.conjugation_report import render_scorecard

SCREEN = "CTRP v2"
ctrp = pd.read_pickle("data/interim/ctrp_canonical.pkl")
clin = clinical_tox_table()
growth = growth_lookup()

sheets = {}
for ind in ("CRC", "HCC"):
    rk = rank_potency(ctrp, ind, emax_max=0.5)
    sc = delivery_scorecard(ctrp, rk, clin, indication=ind, growth_fn=growth, top=20)
    # Phase B: fetch structures + run the chemical gates G4/G5, then G2a efflux
    smiles = fetch_smiles(sc["compound"].tolist())
    sc = add_chemistry(sc, smiles)
    sc = add_efflux(sc, ctrp)
    sc = add_window(sc)          # G6 therapeutic window (on-target, HPA DLT organs)
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
