"""Write data/interim/ctrp_top_drugs.csv = the current CTRP CRC top-20 single
agents (combinations dropped - they are not single-payload candidates). This is
the compound list that scripts/ctrp_export_curves.R pulls raw wells for, so run
this whenever the ranking changes, then re-run the R export and gen_curves.py.

    PYTHONPATH=src py scripts/gen_ctrp_top_drugs.py
"""
import pandas as pd

from warhead.analysis.screen_potency import rank_potency

ctrp = pd.read_pickle("data/interim/ctrp_canonical.pkl")
top = rank_potency(ctrp, "CRC", emax_max=0.5)
singles = top[~top["compound"].str.contains(":")].head(20)
singles[["compound"]].to_csv("data/interim/ctrp_top_drugs.csv", index=False)
print("wrote", len(singles), "compounds:", singles["compound"].tolist())
