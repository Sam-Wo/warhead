"""Exatecan combination-partner scoring from single-agent sensitivity.

True synergy needs combination screening; this instead ranks LIKELY partners by
ORTHOGONALITY / complementary coverage across the cell-line panel - the assumption
behind a rational combo: a good partner is potent on the cells exatecan (a Top1i)
does NOT kill, and its response pattern is not just a copy of exatecan's.

Per candidate compound, over the lines it shares with the exatecan anchor:
- `corr`          Spearman of per-line potency vs exatecan. Low/negative = the two
                  hit different cells (orthogonal) = complementary.
- `pot_resistant` median potency (-log10 IC50) on the exatecan-RESISTANT lines
                  (bottom tertile of exatecan potency). High = covers the blind spot.
- `combo_score`   z(pot_resistant) - z(corr): high = potent-where-exatecan-fails AND
                  orthogonal.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd
from scipy import stats as sstats


def _norm(s) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


# Top1-inhibitor class: exatecan is the lead, the rest are surrogates so the anchor
# generalises across screens (CTRP has no exatecan) and is not tied to one molecule.
TOP1I = ["exatecan-mesylate", "SN-38", "irinotecan", "topotecan", "camptothecin",
         "10-hydroxycamptothecin"]


def top1i_axis(canonical: pd.DataFrame, anchors=TOP1I) -> pd.Series:
    """Consensus per-line Top1i sensitivity: each anchor drug's per-line potency is
    z-scored across lines (removing the large absolute-potency differences between,
    e.g., SN-38 and irinotecan), then averaged per line. Higher = more Top1i-sensitive."""
    anorm = {_norm(a) for a in anchors}
    p = canonical.copy()
    p = p[np.isfinite(p["ic50_nM"]) & (p["ic50_nM"] > 0)]
    p["pot"] = -np.log10(p["ic50_nM"])
    cls = p[p["compound"].map(_norm).isin(anorm)]
    zs = []
    for _, g in cls.groupby("compound"):
        s = g.groupby("model_id")["pot"].median()
        if len(s) >= 20 and s.std(ddof=0) > 0:
            zs.append((s - s.mean()) / s.std(ddof=0))
    if not zs:
        raise ValueError("no Top1i anchors found in this screen")
    return pd.concat(zs, axis=1).mean(axis=1)


def combo_scores(canonical: pd.DataFrame, *, anchors=TOP1I, min_shared=50,
                 resistant_q=1.0 / 3) -> pd.DataFrame:
    ax = top1i_axis(canonical, anchors)                  # consensus Top1i-sensitivity axis
    anorm = {_norm(a) for a in anchors}
    p = canonical.copy()
    p = p[np.isfinite(p["ic50_nM"]) & (p["ic50_nM"] > 0)]
    p["pot"] = -np.log10(p["ic50_nM"])
    res_lines = set(ax[ax <= ax.quantile(resistant_q)].index)   # Top1i-RESISTANT lines

    rows = []
    for cmp, g in p.groupby("compound"):
        if _norm(cmp) in anorm:
            continue
        gv = g.groupby("model_id")["pot"].median()
        shared = ax.index.intersection(gv.index)
        if len(shared) < min_shared:
            continue
        rho = sstats.spearmanr(ax[shared].to_numpy(), gv[shared].to_numpy()).correlation
        on_res = gv[gv.index.isin(res_lines)]
        rows.append({
            "compound": cmp, "target": g["target"].iloc[0],
            "moa": g["moa"].iloc[0] if "moa" in g else np.nan,
            "clinical_phase": g["clinical_phase"].iloc[0] if "clinical_phase" in g else np.nan,
            "n_shared": int(len(shared)), "corr": float(rho),
            "pot_resistant": float(on_res.median()) if len(on_res) else np.nan,
            "median_ic50_resistant_nM": float(10 ** (-on_res.median())) if len(on_res) else np.nan,
            "median_ic50_nM": float(10 ** (-gv.median())),
        })
    df = pd.DataFrame(rows)
    df = df[np.isfinite(df["corr"]) & np.isfinite(df["pot_resistant"])]

    def _z(s):
        return (s - s.mean()) / s.std(ddof=0)
    df["combo_score"] = _z(df["pot_resistant"]) - _z(df["corr"])
    return df.sort_values("combo_score", ascending=False).reset_index(drop=True)


# --- mechanistic potentiation vs antagonism, from target annotation --------------
# Path 2 (true synergy): inhibit the DDR / replication-stress checkpoint that RESOLVES
# Top1i damage -> forces damaged cells through replication into mitotic catastrophe.
DDR_TARGETS = {"ATR", "ATRIP", "CHEK1", "CHEK2", "WEE1", "ATM", "PRKDC", "TOPBP1", "CLSPN",
               "PARP1", "PARP2", "PARP3", "PARP", "RAD51", "BRCA1", "BRCA2", "MRE11", "MRE11A",
               "RAD50", "NBN", "POLQ", "RRM1", "RRM2", "PLK4"}
_DDR_KEYS = ("ATR", "CHK", "CHEK", "WEE1", "PARP", "DNA-PK", "DNAPK", "ATM ", "CHECKPOINT",
             "DNA DAMAGE", "DNA REPAIR", "HOMOLOGOUS RECOMB", "REPLICATION STRESS")
# Path 1 trap (antagonism risk): arrest the cycle -> pull cells OUT of S-phase ->
# protect from the S-phase-specific Top1i (a schedule collision).
ARREST_TARGETS = {"CDK1", "CDK2", "CDK4", "CDK6", "CDK7", "CDK9", "CCNB1", "KIF11", "KSP",
                  "AURKA", "AURKB", "AURKC", "PLK1", "BIRC5", "TUBB", "TUBB1", "TUBB3", "TUBB4A",
                  "TUBB4B", "TUBA1A", "TUBA1B", "TUBA1C", "TUBA3C", "TUBA4A", "TUBA8", "TUBD1"}
_ARREST_KEYS = ("TUBULIN", "MICROTUBULE", "MITOT", "KINESIN", "AURORA", "PLK", "CDK",
                "CYCLIN-DEPENDENT", "SURVIVIN")


def _hits(target, moa, gene_set, keys):
    ts = set(re.split(r"[;,/\s]+", str(target).upper()))
    n_genes = len(ts & gene_set)
    kw = 1 if any(k in (str(target) + " " + str(moa)).upper() for k in keys) else 0
    return n_genes, kw


def mechanism_axis(target, moa=""):
    """(ddr_score, arrest_score): DDR/checkpoint engagement (Path-2 synergy) and
    cell-cycle-arrest liability (Path-1 antagonism). Simple, interpretable target/MOA
    hits - a mechanistic prior, not a data fit."""
    dg, dk = _hits(target, moa, DDR_TARGETS, _DDR_KEYS)
    ag, ak = _hits(target, moa, ARREST_TARGETS, _ARREST_KEYS)
    ddr = min(dg, 3) + dk
    arrest = min(ag, 3) + ak
    return float(ddr), float(arrest)


def slfn11_dependence(canonical: pd.DataFrame, gene_expr: pd.DataFrame, *, gene="SLFN11",
                      low_q=1.0 / 3, min_shared=30) -> pd.DataFrame:
    """Per-compound dependence on SLFN11 - the axis that drives Top1i sensitivity.
    `slfn11_corr` = Spearman(potency, SLFN11 expression): strongly POSITIVE = SLFN11-
    dependent = shares the Top1i resistance mechanism; ~0 or negative = orthogonal
    (works regardless of SLFN11, i.e. on the Top1i-resistant, SLFN11-low cells).
    `ic50_slfn11low_nM` = potency on the SLFN11-low (Top1i-resistant) lines."""
    ex = gene_expr.set_index("ModelID")[gene].dropna()
    low_lines = set(ex[ex <= ex.quantile(low_q)].index)
    p = canonical.copy()
    p = p[np.isfinite(p["ic50_nM"]) & (p["ic50_nM"] > 0)]
    p["pot"] = -np.log10(p["ic50_nM"])
    rows = []
    for cmp, g in p.groupby("compound"):
        gv = g.groupby("model_id")["pot"].median()
        shared = ex.index.intersection(gv.index)
        if len(shared) < min_shared:
            continue
        rho = sstats.spearmanr(gv[shared].to_numpy(), ex[shared].to_numpy()).correlation
        on_low = gv[gv.index.isin(low_lines)]
        rows.append({"compound": cmp, "target": g["target"].iloc[0],
                     "moa": g["moa"].iloc[0] if "moa" in g else np.nan,
                     "n_shared": int(len(shared)), "slfn11_corr": float(rho),
                     "ic50_slfn11low_nM": float(10 ** (-on_low.median())) if len(on_low) else np.nan,
                     "median_ic50_nM": float(10 ** (-gv.median()))})
    return pd.DataFrame(rows)
