"""Label extra indications onto a canonical frame by DepMap Oncotree grouping.

The loaders label CRC / HCC (the program indications); this adds others on demand
- e.g. AML - by resolving each line's DepMap OncotreePrimaryDisease. CTRP lines are
cell-line names (bridged to a ModelID via Model.csv); PRISM lines are already ACH
ModelIDs. Selectivity then compares the labelled lines against the whole rest of
the panel, exactly as it does for CRC / HCC.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


def _norm(s) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def disease_line_set(model_ids, disease_substr: str, depmap_dir="data/raw/depmap") -> set:
    """Return the subset of `model_ids` whose DepMap OncotreePrimaryDisease matches
    `disease_substr` (case-insensitive). Accepts ACH ModelIDs or cell-line names."""
    mdl = pd.read_csv(Path(depmap_dir) / "Model.csv")
    hit = mdl["OncotreePrimaryDisease"].astype(str).str.contains(disease_substr, case=False, na=False)
    members = set(mdl.loc[hit, "ModelID"])
    name2ach = {_norm(v): k for k, v in zip(mdl["ModelID"], mdl["StrippedCellLineName"])
                if pd.notna(v)}
    out = set()
    for m in model_ids:
        if str(m).startswith("ACH-"):
            if m in members:
                out.add(m)
        else:
            a = name2ach.get(_norm(m))
            if a in members:
                out.add(m)
    return out


def label_indication(canonical: pd.DataFrame, line_set: set, label: str) -> pd.DataFrame:
    """Set indication=`label` for rows whose model_id is in `line_set` (others kept)."""
    c = canonical.copy()
    c.loc[c["model_id"].isin(line_set), "indication"] = label
    return c
