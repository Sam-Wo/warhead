"""Conjugation-suitability scorecard (Phase A - delivery survivability).

Re-scores the top free-drug hits against the questions an ADC payload actually has
to answer, using only data already in hand (no chemistry required yet):

- G1 payload-potency bar (config g1.gate): an ADC delivers little payload per cell,
  so payloads must be far more potent than ordinary drugs - sub-nM IC50 in >=20% of
  lines AND median Emax < 0.15 (complete kill). `potency_gap_log10` = how many logs
  the median IC50 sits ABOVE the 1 nM bar (the size of the medchem gap).
- G2b proliferation independence (the HCC lever): Spearman of per-line log10(IC50)
  vs DepMap inferred growth rate. Strong NEGATIVE rho = more potent in fast-dividing
  lines = mitotic-dependent (an auristatin/taxane signature) = cannot carry a
  slow-growing indication. ~0 and non-significant = proliferation-independent = keep.
- Level-0 curated triage: is the chemotype already a known/failed/candidate ADC
  payload, and what DLT organ does the class carry (from analysis.clinical_tox).

G2a efflux (ABCB1/ABCG2) and the chemical gates G4/G5 need data not yet on disk
(DepMap expression; per-compound SMILES + RDKit) and are added in Phase B.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sstats

from ..stats import benjamini_hochberg


def _norm(s) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def growth_lookup(depmap_dir="data/raw/depmap"):
    """Return a fn model_id -> DepMap inferred growth rate. Accepts either a DepMap
    ModelID (ACH-######) or a cell-line name (bridged via Model.csv StrippedCellLineName)."""
    depmap_dir = Path(depmap_dir)
    gr = pd.read_csv(depmap_dir / "CRISPRInferredModelGrowthRate.csv")
    rate_cols = [c for c in gr.columns if c != "ModelID"]
    gr["growth"] = gr[rate_cols].mean(axis=1, skipna=True)
    by_ach = dict(zip(gr["ModelID"], gr["growth"]))
    mdl = pd.read_csv(depmap_dir / "Model.csv")
    name2ach = {_norm(v): k for k, v in zip(mdl["ModelID"], mdl["StrippedCellLineName"])
                if pd.notna(v)}

    def _g(model_id):
        if str(model_id).startswith("ACH-"):
            return by_ach.get(model_id, np.nan)
        ach = name2ach.get(_norm(model_id))
        return by_ach.get(ach, np.nan) if ach else np.nan

    return _g


def _clin_annotation(clin: pd.DataFrame):
    """Token matcher: curated-compound name tokens (>=5 chars) -> annotation row."""
    def _tokens(name):
        return [t for t in re.split(r"[^a-z0-9]+", str(name).lower()) if len(t) >= 5]
    tab = clin.assign(_tok=clin["compound"].map(_tokens))

    def _match(compound):
        k = _norm(compound)
        hit = tab[tab["_tok"].apply(lambda toks: any(t in k for t in toks))]
        if not len(hit):
            return {"adc_payload_status": "", "dlt_organ": "", "clinical_status": ""}
        r = hit.iloc[0]
        return {"adc_payload_status": r.get("adc_payload_status", ""),
                "dlt_organ": r.get("dlt_organ", ""), "clinical_status": r.get("clinical_status", "")}
    return _match


def delivery_scorecard(canonical: pd.DataFrame, ranking: pd.DataFrame, clin: pd.DataFrame,
                       *, indication="CRC", growth_fn=None, top=20,
                       sub_nM=1.0, emax_max=0.15, frac_min=0.20, alpha=0.05) -> pd.DataFrame:
    """canonical: the source's per-line long frame; ranking: rank_potency() output
    for `indication` (drives the top-N selection + medians). Returns one row per
    top-N compound with the delivery-survivability columns."""
    sub = canonical[canonical["indication"] == indication]
    growth_fn = growth_fn or (lambda m: np.nan)
    top_cmps = ranking.head(top)["compound"].tolist()
    match = _clin_annotation(clin)

    rows = []
    for cmp in top_cmps:
        g = sub[sub["compound"] == cmp]
        ic = g["ic50_nM"].to_numpy(dtype=float)
        ic = ic[np.isfinite(ic) & (ic > 0)]
        med_ic50 = float(np.median(ic)) if len(ic) else np.nan
        frac_sub = float(np.mean(ic <= sub_nM)) if len(ic) else np.nan
        med_emax = float(np.nanmedian(g["emax"])) if g["emax"].notna().any() else np.nan
        g1_pass = bool((frac_sub >= frac_min) and (not np.isnan(med_emax)) and (med_emax < emax_max))
        # proliferation dependence is a PAN-PANEL property (power comes from the full
        # doubling-time range), so regress across ALL lines, not just the indication's:
        gall = canonical[canonical["compound"] == cmp]
        gr = gall.assign(_g=gall["model_id"].map(growth_fn), _y=np.log10(gall["ic50_nM"]))
        gr = gr[np.isfinite(gr["_g"]) & np.isfinite(gr["_y"])]
        if len(gr) >= 8:
            rho, p = sstats.spearmanr(gr["_y"], gr["_g"])
        else:
            rho, p = np.nan, np.nan
        rec = {"compound": cmp, "target": g["target"].iloc[0] if len(g) else "",
               "moa": g["moa"].iloc[0] if len(g) and "moa" in g else "",
               "n_lines": int(len(ic)), "median_ic50_nM": med_ic50,
               "frac_subnM": frac_sub, "potency_gap_log10": float(np.log10(med_ic50)) if med_ic50 > 0 else np.nan,
               "median_emax": med_emax, "g1_potency_pass": g1_pass,
               "n_growth": int(len(gr)), "prolif_rho": float(rho) if np.isfinite(rho) else np.nan,
               "prolif_p": float(p) if np.isfinite(p) else np.nan}
        rec.update(match(cmp))
        rows.append(rec)
    df = pd.DataFrame(rows)
    # BH across the tested compounds; independent = NOT significantly proliferation-linked
    mask = df["prolif_p"].notna()
    df["prolif_q"] = np.nan
    if mask.any():
        df.loc[mask, "prolif_q"] = benjamini_hochberg(df.loc[mask, "prolif_p"].to_numpy())
    def _cls(r):
        if not np.isfinite(r["prolif_q"]):
            return "n/a"
        if r["prolif_q"] >= alpha:
            return "independent"
        return "mitotic-dependent" if r["prolif_rho"] < 0 else "prolif-linked (+)"
    df["prolif_class"] = df.apply(_cls, axis=1)
    df["g2b_independent"] = df["prolif_class"] == "independent"
    return df
