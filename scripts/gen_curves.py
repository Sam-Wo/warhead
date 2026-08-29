"""Top-20 dose-response curves for the two wide-window screens:

  reports/ctrp_top20_curves.pdf   measured median+IQR, pooled from the 16-pt raw
                                  wells (data/interim/ctrp_curves.csv) with
                                  pool_measured_curves() to remove grid zig-zag;
  reports/prism_top20_curves.pdf  fitted 4PL with a free lower asymptote (Emax),
                                  params median-aggregated per compound so the
                                  drawn curve and its IC50/EC90 markers agree.

Both rank compounds by CRC potency (rank_potency, emax_max=0.5) for consistency
with the dashboard and summary. Run:

    PYTHONPATH=src py scripts/gen_curves.py
"""
import numpy as np
import pandas as pd

from warhead.analysis.screen_potency import rank_potency
from warhead.reporting.screen_curves import (pool_measured_curves,
                                             render_fitted_curves,
                                             render_measured_curves)

TOP = 20


# ---- CTRP: measured, pooled from raw wells ---------------------------------
ctrp = pd.read_pickle("data/interim/ctrp_canonical.pkl")
ctrp_top = rank_potency(ctrp, "CRC", emax_max=0.5).head(TOP).reset_index(drop=True)

raw = pd.read_csv("data/interim/ctrp_curves.csv")
raw = raw.rename(columns={"drug": "compound", "dose_uM": "conc_uM", "cellid": "model_id"})
raw["viability"] = raw["viability"] / 100.0          # CTRP raw viability is a PERCENT
raw = raw[raw["compound"].isin(ctrp_top["compound"])]

pooled = pool_measured_curves(raw)
maxc = raw.groupby("compound")["conc_uM"].max()
ctrp_summary = pd.DataFrame({
    "compound": ctrp_top["compound"], "target": ctrp_top["target"],
    "ic50_uM": ctrp_top["median_ic50_nM"] / 1e3, "ec90_uM": ctrp_top["median_ec90_nM"] / 1e3,
    "max_conc_uM": ctrp_top["compound"].map(maxc)})
ctrp_summary = ctrp_summary[ctrp_summary["compound"].isin(pooled["compound"])].reset_index(drop=True)
print("CTRP curves:", len(ctrp_summary), "compounds,", pooled["compound"].nunique(), "with pooled data")
print("wrote", render_measured_curves(pooled, ctrp_summary, source="CTRP v2",
                                       out_path="reports/ctrp_top20_curves.pdf"))


# ---- PRISM: fitted 4PL with Emax, params median-aggregated per compound -----
prism = pd.read_pickle("data/interim/prism_canonical.pkl")
prism_top = rank_potency(prism, "CRC", emax_max=0.5).head(TOP).reset_index(drop=True)

params = pd.read_pickle("data/interim/prism_params.pkl")
agg = (params.groupby("name")
       .agg(upper=("upper_limit", "median"), lower=("lower_limit", "median"),
            slope=("slope", "median"), ec50_uM=("ec50", "median"), ic50_uM=("ic50", "median"))
       .reset_index().rename(columns={"name": "compound"}))
prism_summary = prism_top.merge(agg, on="compound", how="inner")
# EC90 from the SAME fit that is drawn, so curve and markers stay consistent
prism_summary["ec90_uM"] = prism_summary["ec50_uM"] * 9.0 ** (1.0 / prism_summary["slope"].abs())
prism_summary = prism_summary[np.isfinite(prism_summary["ec90_uM"]) & (prism_summary["ec90_uM"] < 1e4)]
prism_summary = prism_summary[["compound", "target", "ic50_uM", "ec50_uM", "slope",
                               "upper", "lower", "ec90_uM"]].reset_index(drop=True)
print("PRISM curves:", len(prism_summary), "compounds")
print("wrote", render_fitted_curves(prism_summary, source="PRISM Repurposing (secondary)",
                                     out_path="reports/prism_top20_curves.pdf"))
