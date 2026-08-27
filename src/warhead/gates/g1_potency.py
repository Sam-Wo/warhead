"""G1 - Potency, done properly.

Consumes a per-(compound x line) refit frame (see curves.refit) and applies the
compound-level gate: sub-nM IC50 in a real fraction of lines AND a median Emax
low enough to represent complete kill. IC50 and Emax are kept separate; there is
no AUC anywhere in here by design.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from ..config import load_gates
from .base import GateResult


def aggregate_potency(
    fits: pd.DataFrame,
    *,
    compound_col: str = "compound_id",
    ic50_nM_max: float = 1.0,
    use_qc_pass: bool = True,
) -> pd.DataFrame:
    """Per-compound potency summary from a refit frame.

    ``frac_lines_sub_threshold`` counts lines whose fitted IC50 is below
    ``ic50_nM_max`` (left-censored curves, whose potency exceeds the assay floor,
    count as sub-threshold by construction). ``median_emax`` is over the same
    lines.
    """
    df = fits
    if use_qc_pass and "qc_pass" in df.columns:
        df = df[df["qc_pass"]]
    thr_M = ic50_nM_max * 1e-9

    def _agg(g: pd.DataFrame) -> pd.Series:
        n = len(g)
        sub = int((g["ic50_M"] < thr_M).sum())
        return pd.Series(
            {
                "n_lines": n,
                "n_lines_sub_threshold": sub,
                "frac_lines_sub_threshold": (sub / n) if n else 0.0,
                "median_ic50_nM": float(np.nanmedian(g["ic50_M"]) * 1e9),
                "median_emax": float(np.nanmedian(g["emax"])),
                "frac_left_censored": float((g["ic50_censoring"] == "left").mean()),
            }
        )

    return df.groupby(compound_col).apply(_agg, include_groups=False).reset_index()


def gate_g1(fits: pd.DataFrame, *, compound_col: str = "compound_id", config: dict | None = None) -> GateResult:
    cfg = (config or load_gates())["g1"]["gate"]
    agg = aggregate_potency(
        fits, compound_col=compound_col, ic50_nM_max=cfg["ic50_nM_max"]
    )
    potent = agg["frac_lines_sub_threshold"] >= cfg["ic50_frac_lines_min"]
    complete = agg["median_emax"] < cfg["emax_median_max"]
    agg = agg.assign(passed=(potent & complete))

    reasons = []
    for _, row in agg.iterrows():
        if row["passed"]:
            reasons.append("")
        elif not (row["frac_lines_sub_threshold"] >= cfg["ic50_frac_lines_min"]):
            reasons.append(
                f"sub-nM in only {row['frac_lines_sub_threshold']:.0%} of lines "
                f"(need {cfg['ic50_frac_lines_min']:.0%})"
            )
        else:
            reasons.append(f"median Emax {row['median_emax']:.2f} >= {cfg['emax_median_max']}")
    agg["g1_reason"] = reasons

    passed = agg[agg["passed"]].drop(columns=["passed"]).reset_index(drop=True)
    failed = agg[~agg["passed"]].drop(columns=["passed"]).reset_index(drop=True)
    return GateResult(
        gate="G1",
        passed=passed,
        failed=failed,
        reason_col="g1_reason",
        config=cfg,
        summary={"n_compounds": len(agg), "n_pass": int(passed.shape[0])},
    )
