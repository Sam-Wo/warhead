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


def line_set(model_ids, *, disease_substr: str | None = None, codes=None,
             depmap_dir="data/raw/depmap") -> set:
    """Return the subset of `model_ids` (ACH ModelIDs or cell-line names) that match a
    DepMap OncotreePrimaryDisease substring and/or a set of OncotreeCodes. Use a
    substring for a broad grouping (e.g. 'Acute Myeloid Leukemia') or codes for a
    precise subtype family (e.g. the STAD stomach-adenocarcinoma codes)."""
    mdl = pd.read_csv(Path(depmap_dir) / "Model.csv")
    mask = pd.Series(False, index=mdl.index)
    if disease_substr:
        mask |= mdl["OncotreePrimaryDisease"].astype(str).str.contains(disease_substr, case=False, na=False)
    if codes:
        mask |= mdl["OncotreeCode"].isin(set(codes))
    members = set(mdl.loc[mask, "ModelID"])
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


def disease_line_set(model_ids, disease_substr: str, depmap_dir="data/raw/depmap") -> set:
    """Back-compat wrapper: match by OncotreePrimaryDisease substring."""
    return line_set(model_ids, disease_substr=disease_substr, depmap_dir=depmap_dir)


# stomach-adenocarcinoma OncotreeCodes (gastric); excludes esophageal / GEJ
GASTRIC_CODES = ["STAD", "DSTAD", "ISTAD", "TSTAD", "MSTAD", "SSRCC", "STAS", "GRC"]


def label_indication(canonical: pd.DataFrame, line_set: set, label: str) -> pd.DataFrame:
    """Set indication=`label` for rows whose model_id is in `line_set` (others kept)."""
    c = canonical.copy()
    c.loc[c["model_id"].isin(line_set), "indication"] = label
    return c
