"""Source-agnostic potency ranking + HCC/CRC selectivity for any drug screen.

Every screen loader emits the same canonical long frame (one row per
compound x cell line):

    source, compound, target, moa, model_id, indication, ic50_nM, ec90_nM,
    emax, ec90_extrapolated, clinical_phase

``indication`` is 'CRC' / 'HCC' / 'other'. ``emax`` is the fitted lower asymptote
(residual viability; lower = more complete kill) where the screen provides it,
else NaN. This module then produces the identical outputs (EC90 ranking,
tissue-selectivity) across GDSC / PRISM / CTRP / NCI-60.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as sstats

from ..stats import benjamini_hochberg

CANON_COLS = ["source", "compound", "target", "moa", "model_id", "indication",
              "ic50_nM", "ec90_nM", "emax", "ec90_extrapolated", "clinical_phase"]


def rank_potency(df: pd.DataFrame, indication: str, *, min_lines: int = 5,
                 emax_max: float | None = None) -> pd.DataFrame:
    """Per-compound EC90 summary over the indication's lines, ranked by lowest
    (most potent) median EC90.

    ``emax_max`` (e.g. 0.5) drops compounds whose median residual viability is too
    high to count as a real kill - a filter GDSC's bottom=0 fit could not provide.
    """
    sub = df[df["indication"] == indication]
    cols = ["compound", "target", "moa", "clinical_phase", "n_lines", "median_ec90_nM",
            "median_ic50_nM", "median_emax", "frac_ec90_extrapolated"]
    if not len(sub):
        return pd.DataFrame(columns=cols)

    def _agg(g: pd.DataFrame) -> pd.Series:
        return pd.Series({
            "target": g["target"].iloc[0],
            "moa": g["moa"].iloc[0] if "moa" in g else np.nan,
            "clinical_phase": g["clinical_phase"].iloc[0] if "clinical_phase" in g else np.nan,
            "n_lines": len(g),
            "median_ec90_nM": float(np.nanmedian(g["ec90_nM"])),
            "median_ic50_nM": float(np.nanmedian(g["ic50_nM"])),
            "median_emax": float(np.nanmedian(g["emax"])) if g["emax"].notna().any() else np.nan,
            "frac_ec90_extrapolated": float(g["ec90_extrapolated"].mean()) if g["ec90_extrapolated"].notna().any() else np.nan,
        })

    agg = sub.groupby("compound").apply(_agg, include_groups=False).reset_index()
    agg = agg[agg["n_lines"] >= min_lines]
    if emax_max is not None:
        agg = agg[agg["median_emax"].isna() | (agg["median_emax"] <= emax_max)]
    return agg.sort_values("median_ec90_nM").reset_index(drop=True)


def selectivity(df: pd.DataFrame, indication: str, *, min_in: int = 5, min_out: int = 20,
                potent_ic50_max_nM: float = 1000.0) -> pd.DataFrame:
    """Per-compound potency in-indication vs the rest of the panel, on the fitted
    IC50 (nM). ``delta_potency`` > 0 = more potent in the indication; one-sided
    Mann-Whitney with BH-adjusted q; Cliff's delta effect size."""
    work = df.copy()
    work = work[np.isfinite(work["ic50_nM"]) & (work["ic50_nM"] > 0)]
    work["potency"] = -np.log10(work["ic50_nM"])
    rows = []
    for comp, g in work.groupby("compound"):
        pin = g.loc[g["indication"] == indication, "potency"].to_numpy()
        pout = g.loc[g["indication"] != indication, "potency"].to_numpy()
        if pin.size < min_in or pout.size < min_out:
            continue
        try:
            U, p = sstats.mannwhitneyu(pin, pout, alternative="greater")
            cliffs = 2.0 * U / (pin.size * pout.size) - 1.0
        except ValueError:
            p, cliffs = np.nan, np.nan
        med_in, med_out = float(np.median(pin)), float(np.median(pout))
        rows.append({
            "compound": comp, "target": g["target"].iloc[0],
            "moa": g["moa"].iloc[0] if "moa" in g else np.nan,
            "clinical_phase": g["clinical_phase"].iloc[0] if "clinical_phase" in g else np.nan,
            "n_in": int(pin.size), "n_out": int(pout.size),
            "potency_in": med_in, "potency_out": med_out,
            "delta_potency": med_in - med_out,
            "median_ic50_in_nM": float(10 ** (-med_in)),
            "cliffs_delta": float(cliffs), "p": float(p),
        })
    out = pd.DataFrame(rows)
    if len(out):
        out["q"] = benjamini_hochberg(out["p"].to_numpy())
        out["selective"] = (out["q"] < 0.1) & (out["delta_potency"] > 0)
        out["selective_potent"] = out["selective"] & (out["median_ic50_in_nM"] <= potent_ic50_max_nM)
        out = out.sort_values(["selective_potent", "delta_potency"],
                              ascending=[False, False]).reset_index(drop=True)
    return out
