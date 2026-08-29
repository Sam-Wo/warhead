"""GDSC EC90 potency ranking and HCC/CRC tissue-selectivity.

Two questions, per indication (HCC = LIHC, CRC = COREAD):
  1. Which compounds reach 90% effect at the lowest concentration in the
     indication's cell lines? -> ``indication_ranking``.
  2. Which compounds are *selectively* potent in the indication vs every other
     lineage - the signal for a mechanistic advantage? -> ``selectivity``.

Potency axis is ``-log10(EC90 [uM])`` (higher = more potent). Selectivity is the
median potency in-indication minus median potency in the rest, with a Mann-Whitney
test and BH-adjusted q. Caveats travel with the numbers: EC90 is extrapolated for
most GDSC curves (see io.gdsc), LIHC has only ~15 lines, and cell-line tissue
selectivity is hypothesis-generating - it can reflect a co-enriched dependency
rather than a tissue-intrinsic mechanism.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as sstats

from ..io.gdsc import INDICATION_TCGA
from ..stats import benjamini_hochberg


def _potency(ec90_uM: pd.Series) -> pd.Series:
    return -np.log10(ec90_uM)


def indication_ranking(df: pd.DataFrame, indication: str, *, min_lines: int = 5) -> pd.DataFrame:
    """Per-compound EC90 summary over the indication's cell lines, ranked by the
    lowest (most potent) median EC90."""
    tcga = INDICATION_TCGA[indication]
    sub = df[df["tcga_desc"] == tcga]

    def _agg(g: pd.DataFrame) -> pd.Series:
        n = len(g)
        return pd.Series({
            "target": g["target"].iloc[0],
            "pathway": g["pathway"].iloc[0],
            "n_lines": n,
            "median_ec90_uM": float(np.nanmedian(g["ec90_uM"])),
            "median_ic50_uM": float(np.nanmedian(g["ic50_uM"])),
            "median_auc": float(np.nanmedian(g["auc"])),
            "frac_ec90_within_range": float((g["ec90_range"] == "within").mean()),
            "frac_ic50_within_range": float((g["ic50_range"] == "within").mean()),
        })

    agg = sub.groupby("drug_name").apply(_agg, include_groups=False).reset_index()
    agg = agg[agg["n_lines"] >= min_lines]
    return agg.sort_values("median_ec90_uM").reset_index(drop=True)


def selectivity(df: pd.DataFrame, indication: str, *, min_in: int = 5, min_out: int = 20,
                potent_ic50_max_uM: float = 1.0) -> pd.DataFrame:
    """Per-compound potency in-indication vs the rest of the panel.

    Potency here is ``-log10(IC50 [uM])`` - the DIRECTLY fitted value, not the
    (mostly extrapolated) EC90 - so the selectivity call does not ride on model
    extrapolation. ``delta_potency`` > 0 means more potent in the indication; a
    one-sided Mann-Whitney tests whether in-indication potency exceeds the rest,
    with BH-adjusted q. (EC90 gives an almost identical ranking; it drives the
    potency panels A/B, IC50 drives selectivity C/D.)
    """
    tcga = INDICATION_TCGA[indication]
    work = df.copy()
    work["potency"] = _potency(work["ic50_uM"])
    work = work[np.isfinite(work["potency"])]

    rows = []
    for drug, g in work.groupby("drug_name"):
        pin = g.loc[g["tcga_desc"] == tcga, "potency"].to_numpy()
        pout = g.loc[g["tcga_desc"] != tcga, "potency"].to_numpy()
        if pin.size < min_in or pout.size < min_out:
            continue
        # one-sided: in-indication MORE potent than rest
        try:
            U, p = sstats.mannwhitneyu(pin, pout, alternative="greater")
            # Cliff's delta from U (scale-free effect size in [-1, 1]).
            cliffs = 2.0 * U / (pin.size * pout.size) - 1.0
        except ValueError:
            p, cliffs = np.nan, np.nan
        med_in, med_out = float(np.median(pin)), float(np.median(pout))
        rows.append({
            "drug_name": drug,
            "target": g["target"].iloc[0],
            "pathway": g["pathway"].iloc[0],
            "n_in": int(pin.size), "n_out": int(pout.size),
            "potency_in": med_in,            # -log10 IC50 (uM), in indication
            "potency_out": med_out,
            "delta_potency": med_in - med_out,   # >0 = more potent in indication
            "median_ic50_in_uM": float(10 ** (-med_in)),
            "cliffs_delta": float(cliffs),
            "p": float(p),
        })
    out = pd.DataFrame(rows)
    if len(out):
        out["q"] = benjamini_hochberg(out["p"].to_numpy())
        # Selective AND actually potent: a tissue difference on an inactive
        # compound (IC50 ~ mM) is noise, not a mechanistic advantage.
        out["selective"] = (out["q"] < 0.1) & (out["delta_potency"] > 0)
        out["selective_potent"] = out["selective"] & (out["median_ic50_in_uM"] <= potent_ic50_max_uM)
        out = out.sort_values(["selective_potent", "delta_potency"],
                              ascending=[False, False]).reset_index(drop=True)
    return out
