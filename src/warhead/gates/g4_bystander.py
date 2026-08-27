"""G4 - Bystander competence. (Deferred: build order step 7.)

Do NOT rebuild the bystander predictor from scratch. Prior art: Guo et al.,
Adv. Sci. 2024 (10.1002/advs.202306309) - graph-attention model on membrane
permeability, validated against 80+ ADCdb payloads, B-score threshold 1.5.

Layer on top (RDKit over the G1-passing set intersected with COCONUT/NPAtlas/
ChEMBL):
  * charge state at lysosomal pH ~4.8 vs cytosol 7.2 (the MMAE/MMAF split is
    essentially this),
  * cLogD(7.4), TPSA, MW.

Bystander is a DESIGN CHOICE, not a universal good. Tag, do not filter - let the
antigen program decide.
"""
from __future__ import annotations

import pandas as pd

from ..config import load_gates


def tag_bystander(payload_props: pd.DataFrame, *, config: dict | None = None) -> pd.DataFrame:
    """Tag (not filter) each payload as bystander-competent using physchem
    windows from gates.yaml. Placeholder for the Guo B-score integration."""
    cfg = (config or load_gates())["g4"]
    lo, hi = cfg["clogd_range"]
    df = payload_props.copy()
    if {"clogd", "tpsa", "mw"}.issubset(df.columns):
        df["bystander_tag"] = (
            df["clogd"].between(lo, hi)
            & (df["tpsa"] <= cfg["tpsa_max"])
            & (df["mw"] <= cfg["mw_max"])
        )
    else:
        raise NotImplementedError(
            "Provide clogd/tpsa/mw (RDKit) or a Guo B-score column; see docstring."
        )
    return df
