"""G6 - Therapeutic window. (Deferred: build order step 8.)

Payload class determines DLT organ far more than antigen does (ocular=MMAF,
ILD=DXd, thrombocytopenia=maytansinoids, neuropathy=MMAE).

  1. FAERS disproportionality (ROR/PRR with shrinkage) restricted to ADC
     regimens, stratified by payload class -> empirical class->tox map.
  2. Score GTEx/HPA expression of each candidate payload TARGET across the five
     recurring DLT compartments (HSC/marrow, GI crypt, cornea, alveolar type II,
     peripheral nerve). Low across all five = window signal.
"""
from __future__ import annotations

import pandas as pd

from ..config import load_gates


def score_window(target_expression: pd.DataFrame, *, config: dict | None = None) -> pd.DataFrame:
    """``target_expression``: rows = payload targets, columns = the five DLT
    compartments (expression percentiles). Flags targets that are low across all
    five."""
    cfg = (config or load_gates())["g6"]
    comps = cfg["dlt_compartments"]
    missing = set(comps) - set(target_expression.columns)
    if missing:
        raise NotImplementedError(f"missing DLT compartment columns: {sorted(missing)}")
    df = target_expression.copy()
    thr = cfg["expression_percentile_max"]
    df["window_ok"] = (df[comps] <= thr).all(axis=1)
    return df
