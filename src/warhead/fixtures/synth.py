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

# --- compound catalogue: (id, class, base log10 IC50 @ ref dt, prolif_beta, emax, efflux_beta)
# prolif_beta in log10(IC50) units per hour of doubling time.
_CATALOGUE = [
    # Strong mitotic dependence (antimitotics): steep positive slope -> fail G2b.
    ("MMAE_like",          "auristatin",            -9.5,  0.0200, 0.03, 0.35),
    ("DM1_like",           "maytansinoid",          -9.3,  0.0120, 0.05, 0.45),
    # Intermediate (Top1i): partial dependence -> fail G2b, but less extreme.
    ("exatecan_like",      "topoisomerase_I",       -9.2,  0.0022, 0.08, 0.05),
    ("SN38_like",          "topoisomerase_I",       -9.0,  0.0020, 0.10, 0.05),
    # Non-mitotic payload classes: ~flat -> pass G2b (the ones we want).
    ("amanitin_like",      "RNAPII_inhibitor",      -9.2,  0.0002, 0.06, 0.00),
    ("PF846_like",         "translation_inhibitor", -9.4,  0.0000, 0.05, 0.00),
    ("thailanstatin_like", "spliceosome",           -9.1,  0.0003, 0.07, 0.00),
    ("degrader_like",      "protein_degrader",      -9.3,  0.0002, 0.05, 0.00),
    # Efflux substrate: intermediate proliferation slope + strong ABCB1 signal
    # (fails G2a; here it also carries a modest proliferation slope).
    ("MDR_substrate_like", "efflux_control",        -9.0,  0.0025, 0.06, 0.60),
    # Weak binders: not sub-nM -> filtered at G1, never reach G2b by default.
    ("weak_binder_1",      "non_payload",           -7.0,  0.0020, 0.35, 0.05),
    ("weak_binder_2",      "non_payload",           -6.5,  0.0008, 0.45, 0.05),
    ("flat_but_weak",      "non_payload",           -7.5,  0.0001, 0.20, 0.00),
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
    # Doubling time skewed toward fast growth (gamma), clipped to a plausible range.
    dt = np.clip(rng.gamma(shape=2.4, scale=16.0, size=n) + 16.0, 16.0, 130.0)
    lineage = rng.choice(_LINEAGES, size=n)
    ids = [f"ACH-{i:06d}" for i in rng.choice(np.arange(1, 999999), size=n, replace=False)]
    names = [f"LINE{i:03d}" for i in range(n)]
    abcb1 = rng.normal(4.0, 1.6, size=n).clip(0)   # log2 TPM+1
    abcg2 = rng.normal(3.0, 1.4, size=n).clip(0)
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
        }
    )


def _dose_grid() -> np.ndarray:
    # 8-point, 0.1 nM -> 10 uM (matches PRISM secondary reach; floor sits above
    # the most potent lines so some curves are genuinely left-censored).
    return np.logspace(-10, -5, 8)


def make(seed: int = 20260827, n_lines: int = 60) -> SyntheticData:
    rng = np.random.default_rng(seed)
    models = _make_models(rng, n_lines)
    doses = _dose_grid()

    dr_rows = []
    truth_rows = []
    for cid, cls, base, beta, emax, efflux_beta in _CATALOGUE:
        truth_rows.append(
            {"compound_id": cid, "class": cls, "base_log10_ic50": base,
             "prolif_beta": beta, "emax_true": emax, "efflux_beta": efflux_beta}
        )
        abcb1_c = models["_ABCB1"].to_numpy()
        abcb1_z = (abcb1_c - abcb1_c.mean()) / (abcb1_c.std() + 1e-9)
        for i, m in models.reset_index(drop=True).iterrows():
            # Planted per-line log10 IC50: base + proliferation term + efflux term + noise.
            log_ic50 = (
                base
                + beta * (m["doubling_time_hours"] - _DT_REF)
                + efflux_beta * abcb1_z[i]
                + rng.normal(0, 0.08)
            )
            hill = rng.uniform(0.9, 1.4)
            v = four_pl(doses, top=1.0, emax=emax, log10_ic50=log_ic50, hill=hill)
            v = v + rng.normal(0, 0.045, size=doses.size)
            v = np.clip(v, 0.0, 1.08)
            for d, vi in zip(doses, v):
                dr_rows.append({"compound_id": cid, "ModelID": m["ModelID"],
                                "dose_M": float(d), "viability": float(vi)})

    dose_response = pd.DataFrame(dr_rows)

    expression = (
        models.melt(
            id_vars="ModelID", value_vars=["_ABCB1", "_ABCG2"],
            var_name="gene", value_name="expression",
        ).assign(gene=lambda d: d["gene"].str.lstrip("_"))
    )

    chronos, copy_number = _make_collateral(rng, models["ModelID"].tolist())

    model_meta = models.drop(columns=["_ABCB1", "_ABCG2"])
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
