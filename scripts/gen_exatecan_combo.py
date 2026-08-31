"""Exatecan (Top1-inhibitor) combination-partner analysis for the deck.

Ranks LIKELY synergistic partners by three orthogonality axes (true synergy needs
the wetlab combination screen, which anchors this):
  1. complementary coverage across the panel  (PRISM + CTRP, reproducible)
  2. resistance orthogonality on the SLFN11 axis  (DepMap SLFN11 x PRISM)
  3. MOA orthogonality  (Tahoe-100M transcriptional distance from Top1i)

Writes reports/exatecan_partners_{map,slfn11,moa}.pdf + reports/exatecan_partners.xlsx.

    PYTHONPATH=src py scripts/gen_exatecan_combo.py
Needs: data/interim/depmap_genes.csv (scripts/fetch_depmap_genes.py) and
data/raw/tahoe/{gene_set_library_up_crisp,..._dn_crisp}.gmt.gz + drug_metadata.parquet.
"""
import re

import numpy as np
import pandas as pd

from warhead.analysis.exatecan_combo import TOP1I, combo_scores, slfn11_dependence
from warhead.analysis.tahoe_moa import moa_distance_table
from warhead.io.pubchem import fetch_smiles
from warhead.reporting.exatecan_combo_plot import (render_combo_consensus,
                                                   render_moa_orthogonality,
                                                   render_slfn11_orthogonality)


def _norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


prism = pd.read_pickle("data/interim/prism_canonical.pkl")
ctrp = pd.read_pickle("data/interim/ctrp_canonical.pkl")

# --- axis 1: complementary coverage, PRISM + CTRP consensus ---
P = combo_scores(prism); P["k"] = P["compound"].map(_norm)
C = combo_scores(ctrp); C["k"] = C["compound"].map(_norm)
m = P.merge(C, on="k", suffixes=("_p", "_c")).rename(
    columns={"compound_p": "compound", "target_p": "target", "moa_p": "moa"})
m["consensus"] = m["combo_score_p"] + m["combo_score_c"]
m = m.sort_values("consensus", ascending=False).reset_index(drop=True)
print("axis 1: complementary coverage -", len(m), "compounds scored in both screens")
print("  wrote", render_combo_consensus(m, out_path="reports/exatecan_partners_map.pdf"))

# --- axis 2: SLFN11 resistance orthogonality ---
ge = pd.read_csv("data/interim/depmap_genes.csv")
dep = slfn11_dependence(prism, ge)
partners = m.head(12)["compound"].tolist()
print("  wrote", render_slfn11_orthogonality(
    dep, anchor_norms=[_norm(x) for x in TOP1I], partner_names=partners,
    out_path="reports/exatecan_partners_slfn11.pdf"))

# --- axis 3: Tahoe MOA orthogonality (bridge combo <-> Tahoe via PubChem CID) ---
moa = moa_distance_table("data/raw/tahoe/gene_set_library_up_crisp.gmt.gz",
                         "data/raw/tahoe/gene_set_library_dn_crisp.gmt.gz")
anchor_dist = float(moa[moa["tahoe_drug"].map(_norm) == _norm("Topotecan (hydrochloride)")]["moa_distance"].iloc[0])
dm = pd.read_parquet("data/raw/tahoe/drug_metadata.parquet")
dm["cid"] = pd.to_numeric(dm["pubchem_cid"], errors="coerce"); dm["k"] = dm["drug"].map(_norm)
moa["k"] = moa["tahoe_drug"].map(_norm)
moa = moa.merge(dm[["k", "cid"]], on="k", how="left")
fetch_smiles(m["compound"].tolist())                         # populate the CID cache
cache = pd.read_csv("data/interim/pubchem_smiles.csv")
name2cid = {_norm(n): c for n, c in zip(cache["name"], cache["cid"]) if pd.notna(c)}
m["cid"] = m["compound"].map(lambda n: name2cid.get(_norm(n)))
matched = (m.dropna(subset=["cid"])
           .merge(moa.dropna(subset=["cid"])[["cid", "moa_distance", "tahoe_drug"]], on="cid", how="inner")
           .rename(columns={"consensus": "combo"}))
print(f"axis 3: MOA orthogonality - {len(matched)} partners bridged to Tahoe (anchor dist {anchor_dist:.2f})")
print("  wrote", render_moa_orthogonality(matched, anchor_dist=anchor_dist,
                                          out_path="reports/exatecan_partners_moa.pdf"))

# --- ranked workbook ---
out = m[["compound", "target", "combo_score_p", "combo_score_c", "consensus",
         "corr_p", "corr_c", "median_ic50_resistant_nM_p"]].copy()
out = out.merge(dep[["compound", "slfn11_corr"]], on="compound", how="left")
out = out.merge(matched[["compound", "moa_distance"]], on="compound", how="left")
with pd.ExcelWriter("reports/exatecan_partners.xlsx", engine="openpyxl") as xw:
    out.to_excel(xw, sheet_name="partners", index=False)
print("wrote reports/exatecan_partners.xlsx")
