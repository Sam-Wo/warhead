"""Synthetic ADC-payload screening data with planted ground truth.

The point of this fixture is NOT to look real - it is to encode a known answer so
the pipeline can be validated end to end:

  * each compound has a true ``prolif_beta`` = how its log10(IC50) shifts per hour
    of doubling time. Antimitotic classes get a strong positive beta (lose potency
    in slow lines); non-mitotic classes get ~0. That dependence is baked into each
    line's IC50, a 4PL curve is simulated at that IC50, and the refit engine must
    recover it well enough for G2b to re-detect the planted beta.
  * pooled-PRISM under-representation of slow lines is mimicked by sampling
    doubling times skewed toward fast growth, so the G2b weighting has something
    to correct.
  * an efflux axis (ABCB1) and a POLR2A collateral-lethality signal are planted
    for G2a and G2c.

Determinism: everything is drawn from a seeded numpy Generator.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from ..curves.refit import four_pl
from ..stats import balance_weights


def _zscore(a: np.ndarray) -> np.ndarray:
    a = np.asarray(a, dtype=float)
    return (a - a.mean()) / (a.std() + 1e-9)


def _wremove(v: np.ndarray, u: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Remove the w-weighted projection of v onto u (weighted-centered)."""
    W = w.sum()
    vc = v - np.sum(w * v) / W
    uc = u - np.sum(w * u) / W
    coef = np.sum(w * vc * uc) / (np.sum(w * uc * uc) + 1e-12)
    return v - coef * uc


def _orthogonalise(cols: list[np.ndarray], anchor: np.ndarray, weights: np.ndarray) -> list[np.ndarray]:
    """Make each column w-orthogonal to ``anchor`` (doubling time) under the SAME
    balance weighting G2b uses, and pairwise-orthogonal to each other.

    The fixture plants a separate biological axis per gate (doubling time -> G2b,
    ABCB1 -> G2a, SLFN11 -> G3b). If SLFN11/ABCB1 have any WEIGHTED covariance with
    doubling time, a compound's SLFN11-driven potency leaks into G2b's weighted
    slope and flips its call. Removing it under the exact weighting G2b applies
    makes each gate isolable and seed-robust.
    """
    basis_u = [anchor]
    out = []
    for c in cols:
        v = _wremove(c, anchor, weights)          # kill weighted cov with doubling time
        for b in out:                              # decorrelate from earlier covariates
            v = v - (v @ b) / (b @ b + 1e-12) * b
        v = _wremove(v, anchor, weights)           # restore doubling-time orthogonality
        out.append(v)
    return [_zscore(v) for v in out]

# --- compound catalogue:
#   (id, class, base log10 IC50 @ ref dt, prolif_beta, emax, efflux_beta, slfn11_beta)
# prolif_beta : log10(IC50) shift per hour of doubling time (mitotic dependence).
# efflux_beta : log10(IC50) shift per z-unit of ABCB1 (efflux substrate strength).
# slfn11_beta : log10(IC50) shift per z-unit of SLFN11.
#   negative -> more potent when SLFN11 is HIGH (Top1i-like; exatecan needs SLFN11)
#   positive -> more potent when SLFN11 is LOW  (orthogonal partner; kills where
#              exatecan cannot -> the ATRi hypothesis)
_CATALOGUE = [
    # Strong mitotic dependence (antimitotics): steep positive slope -> fail G2b.
    # (DM1 beta is high because its large efflux variance dilutes the correlation.)
    ("MMAE_like",          "auristatin",            -9.5,  0.0200, 0.03, 0.35,  0.000),
    ("DM1_like",           "maytansinoid",          -9.3,  0.0160, 0.05, 0.45,  0.000),
    # Top1i: intermediate proliferation slope; potency tracks SLFN11 (slfn11_beta<0).
    ("exatecan_like",      "topoisomerase_I",       -9.2,  0.0022, 0.08, 0.05, -0.045),
    ("SN38_like",          "topoisomerase_I",       -9.1,  0.0022, 0.10, 0.05, -0.040),
    # Non-mitotic payload classes: flat (exact zero slope) -> pass G2b.
    ("amanitin_like",      "RNAPII_inhibitor",      -9.2,  0.0000, 0.06, 0.00,  0.000),
    ("PF846_like",         "translation_inhibitor", -9.4,  0.0000, 0.05, 0.00,  0.000),
    ("thailanstatin_like", "spliceosome",           -9.1,  0.0000, 0.07, 0.00,  0.000),
    ("degrader_like",      "protein_degrader",      -9.3,  0.0000, 0.05, 0.00,  0.000),
    # ATR inhibitor: the orthogonal-partner hypothesis. Potent on SLFN11-low
    # (Top1i-resistant) lines -> slfn11_beta > 0. Proliferation-flat -> passes G2b.
    ("ATRi_like",          "ATR_inhibitor",         -9.3,  0.0000, 0.06, 0.05,  0.050),
    # Efflux substrate: intermediate proliferation slope + strong ABCB1 signal
    # (fails G2a; here it also carries a modest proliferation slope).
    ("MDR_substrate_like", "efflux_control",        -9.0,  0.0025, 0.06, 0.60,  0.000),
    # Weak binders: not sub-nM -> filtered at G1, never reach G2b/G3b by default.
    ("weak_binder_1",      "non_payload",           -7.0,  0.0020, 0.35, 0.05,  0.000),
    ("weak_binder_2",      "non_payload",           -6.5,  0.0008, 0.45, 0.05,  0.000),
    ("flat_but_weak",      "non_payload",           -7.5,  0.0001, 0.20, 0.00,  0.000),
]

_LINEAGES = ["Colorectal", "Colorectal", "Liver", "Lung", "Breast", "Pancreas", "Ovary", "Skin"]
_ONCOTREE = {"Colorectal": "COAD", "Liver": "LIHC", "Lung": "LUAD", "Breast": "BRCA",
             "Pancreas": "PAAD", "Ovary": "OV", "Skin": "SKCM"}
_DT_REF = 50.0  # reference doubling time (hours)


@dataclass
class SyntheticData:
    model_meta: pd.DataFrame        # ModelID, CellLineName, OncotreeLineage, doubling_time_hours, ...
    dose_response: pd.DataFrame     # compound_id, ModelID, dose_M, viability  (tidy long)
    expression: pd.DataFrame        # ModelID, gene, expression  (ABCB1/ABCG2)
    chronos: pd.DataFrame           # ModelID, gene, chronos
    copy_number: pd.DataFrame       # ModelID, gene, cn_log2
    compound_truth: pd.DataFrame    # compound_id, class, prolif_beta, ...


def _make_models(rng: np.random.Generator, n: int) -> pd.DataFrame:
    # Doubling time: full 20-120 h coverage with a mild skew toward fast growth
    # (pooled PRISM under-represents slow lines, but not so extremely that a
    # single slow line dominates the balance-weighted regression).
    dt = 20.0 + 100.0 * rng.beta(1.6, 2.2, size=n)
    lineage = rng.choice(_LINEAGES, size=n)
    ids = [f"ACH-{i:06d}" for i in rng.choice(np.arange(1, 999999), size=n, replace=False)]
    names = [f"LINE{i:03d}" for i in range(n)]
    # Draw raw covariates, then orthogonalise each against doubling time (under
    # G2b's balance weighting) and each other so no gate's signal leaks elsewhere.
    w = balance_weights(dt, n_bins=5)
    abcb1_o, abcg2_o, slfn11_o = _orthogonalise(
        [rng.normal(size=n), rng.normal(size=n), rng.normal(size=n)], anchor=dt, weights=w
    )
    # Means kept high enough that clip(0) never triggers, so it cannot reintroduce
    # correlation with doubling time (which would leak into G2b).
    abcb1 = (6.5 + 1.5 * abcb1_o).clip(0)    # log2 TPM+1
    abcg2 = (5.5 + 1.3 * abcg2_o).clip(0)
    slfn11 = (5.5 + 1.5 * slfn11_o).clip(0)  # Top1i-sensitising axis (G3b)
    return pd.DataFrame(
        {
            "ModelID": ids,
            "CellLineName": names,
            "StrippedCellLineName": names,
            "OncotreeLineage": lineage,
            "OncotreeCode": [_ONCOTREE[l] for l in lineage],
            "doubling_time_hours": dt,
            "_ABCB1": abcb1,
            "_ABCG2": abcg2,
            "_SLFN11": slfn11,
        }
    )


def _dose_grid() -> np.ndarray:
    # 8-point, 0.1 nM -> 10 uM (matches PRISM secondary reach; floor sits above
    # the most potent lines so some curves are genuinely left-censored).
    return np.logspace(-10, -5, 8)


def make(seed: int = 20260827, n_lines: int = 80) -> SyntheticData:
    rng = np.random.default_rng(seed)
    models = _make_models(rng, n_lines)
    doses = _dose_grid()

    abcb1_z = _zscore(models["_ABCB1"].to_numpy())
    slfn11_z = _zscore(models["_SLFN11"].to_numpy())

    dr_rows = []
    truth_rows = []
    for cid, cls, base, beta, emax, efflux_beta, slfn11_beta in _CATALOGUE:
        truth_rows.append(
            {"compound_id": cid, "class": cls, "base_log10_ic50": base,
             "prolif_beta": beta, "emax_true": emax, "efflux_beta": efflux_beta,
             "slfn11_beta": slfn11_beta}
        )
        for i, m in models.reset_index(drop=True).iterrows():
            # Planted per-line log10 IC50:
            #   base + proliferation + efflux + SLFN11(Top1i axis) + noise.
            log_ic50 = (
                base
                + beta * (m["doubling_time_hours"] - _DT_REF)
                + efflux_beta * abcb1_z[i]
                + slfn11_beta * slfn11_z[i]
                + rng.normal(0, 0.05)
            )
            hill = rng.uniform(0.9, 1.4)
            v = four_pl(doses, top=1.0, emax=emax, log10_ic50=log_ic50, hill=hill)
            v = v + rng.normal(0, 0.03, size=doses.size)
            v = np.clip(v, 0.0, 1.08)
            for d, vi in zip(doses, v):
                dr_rows.append({"compound_id": cid, "ModelID": m["ModelID"],
                                "dose_M": float(d), "viability": float(vi)})

    dose_response = pd.DataFrame(dr_rows)

    expression = (
        models.melt(
            id_vars="ModelID", value_vars=["_ABCB1", "_ABCG2", "_SLFN11"],
            var_name="gene", value_name="expression",
        ).assign(gene=lambda d: d["gene"].str.lstrip("_"))
    )

    chronos, copy_number = _make_collateral(rng, models["ModelID"].tolist())

    model_meta = models.drop(columns=["_ABCB1", "_ABCG2", "_SLFN11"])
    compound_truth = pd.DataFrame(truth_rows)
    return SyntheticData(model_meta, dose_response, expression, chronos, copy_number, compound_truth)


def _make_collateral(rng: np.random.Generator, model_ids: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Plant a POLR2A positive control: lines with POLR2A CN loss are more
    dependent (lower Chronos). ME2 gets a weaker signal; other genes are null."""
    genes = ["POLR2A", "ME2", "TP53", "GAPDH", "ACTB", "KRAS"]
    # CN loss signal: strong for POLR2A, moderate for ME2, ~none otherwise.
    dep_effect = {"POLR2A": -0.55, "ME2": -0.30}
    ch_rows, cn_rows = [], []
    for gene in genes:
        loss = rng.random(len(model_ids)) < 0.35  # ~35% of lines hemizygous
        cn = np.where(loss, rng.normal(-0.55, 0.15, len(model_ids)),
                      rng.normal(0.02, 0.15, len(model_ids)))
        base_dep = {"POLR2A": -0.6, "ME2": -0.1, "TP53": 0.05,
                    "GAPDH": -1.1, "ACTB": -1.0, "KRAS": -0.2}[gene]
        eff = dep_effect.get(gene, 0.0)
        chronos = base_dep + eff * loss.astype(float) + rng.normal(0, 0.12, len(model_ids))
        for mid, c, ch in zip(model_ids, cn, chronos):
            cn_rows.append({"ModelID": mid, "gene": gene, "cn_log2": float(c)})
            ch_rows.append({"ModelID": mid, "gene": gene, "chronos": float(ch)})
    return pd.DataFrame(ch_rows), pd.DataFrame(cn_rows)
