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
from warhead.analysis.conjugation import delivery_scorecard, growth_lookup
from warhead.analysis.screen_potency import rank_potency
from warhead.reporting.conjugation_report import render_scorecard

SCREEN = "CTRP v2"
ctrp = pd.read_pickle("data/interim/ctrp_canonical.pkl")
clin = clinical_tox_table()
growth = growth_lookup()

sheets = {}
for ind in ("CRC", "HCC"):
    rk = rank_potency(ctrp, ind, emax_max=0.5)
    sc = delivery_scorecard(ctrp, rk, clin, indication=ind, growth_fn=growth, top=20)
    out = render_scorecard(sc, out_path=f"reports/conjugation_scorecard_{ind}.pdf",
                           indication=ind, screen=SCREEN)
    print("wrote", out, f"| G1 pass {sc.g1_potency_pass.sum()}/{len(sc)} | "
          f"G2b independent {sc.g2b_independent.sum()}/{len(sc)}")
    sheets[ind] = sc

with pd.ExcelWriter("reports/conjugation_scorecard.xlsx", engine="openpyxl") as xw:
    for ind, sc in sheets.items():
        sc.to_excel(xw, sheet_name=ind, index=False)
print("wrote reports/conjugation_scorecard.xlsx")
