"""Compound and cell-line identity resolution.

Non-negotiable (WARHEAD.md sec 4 Conventions):
  * Compound identity -> InChIKey at ingest. Cross-screen name matching
    otherwise quietly loses 20-30% of overlaps.
  * Cell-line identity -> DepMap ModelID (ACH-######) at ingest.

RDKit is optional. When it is absent, ``to_inchikey`` returns None and warns;
downstream code must treat a null InChIKey as "unresolved", never as a match.
"""
from __future__ import annotations

import re
import warnings
from typing import Iterable

import pandas as pd

try:  # pragma: no cover - exercised only when rdkit is installed
    from rdkit import Chem
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")
    _HAVE_RDKIT = True
except Exception:  # pragma: no cover
    _HAVE_RDKIT = False


ACH_RE = re.compile(r"^ACH-\d{6}$")


def have_rdkit() -> bool:
    return _HAVE_RDKIT


def to_inchikey(smiles: str | None) -> str | None:
    """Canonical InChIKey from a SMILES string, or None if unresolvable.

    Requires RDKit. Without it, returns None and warns once per session so that
    callers never silently treat unresolved compounds as identical.
    """
    if smiles is None or (isinstance(smiles, float)):
        return None
    smiles = str(smiles).strip()
    if not smiles:
        return None
    if not _HAVE_RDKIT:
        warnings.warn(
            "RDKit not installed: InChIKey resolution unavailable "
            "(pip install 'warhead[chem]'). Compounds left unresolved.",
            RuntimeWarning,
            stacklevel=2,
        )
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToInchiKey(mol)


def resolve_compound_inchikeys(
    df: pd.DataFrame, smiles_col: str = "smiles", out_col: str = "inchikey"
) -> pd.DataFrame:
    """Add an ``inchikey`` column derived from a SMILES column (idempotent)."""
    out = df.copy()
    out[out_col] = out[smiles_col].map(to_inchikey)
    return out


def normalise_cellline_name(name: str) -> str:
    """Aggressive normalisation for cell-line name matching: uppercase,
    strip everything that is not alphanumeric. 'NCI-H23' -> 'NCIH23'."""
    return re.sub(r"[^A-Z0-9]", "", str(name).upper())


def build_model_alias_map(model_meta: pd.DataFrame) -> dict[str, str]:
    """Map normalised aliases -> ModelID from a DepMap Model metadata frame.

    Uses ``ModelID`` plus any of the common name columns present
    (``CellLineName``, ``StrippedCellLineName``, ``cell_line_name``, ``aliases``).
    """
    name_cols = [
        c
        for c in ("CellLineName", "StrippedCellLineName", "cell_line_name", "aliases")
        if c in model_meta.columns
    ]
    alias: dict[str, str] = {}
    for _, row in model_meta.iterrows():
        mid = row["ModelID"]
        for col in name_cols:
            val = row.get(col)
            if pd.isna(val):
                continue
            for piece in str(val).split(";"):
                key = normalise_cellline_name(piece)
                if key:
                    alias.setdefault(key, mid)
    return alias


def resolve_model_ids(
    df: pd.DataFrame,
    name_col: str,
    model_meta: pd.DataFrame,
    out_col: str = "ModelID",
) -> pd.DataFrame:
    """Resolve a free-text cell-line name column to DepMap ModelID.

    Names already in ACH-###### form pass through unchanged.
    """
    alias = build_model_alias_map(model_meta)
    out = df.copy()

    def _resolve(v: object) -> str | None:
        s = str(v).strip()
        if ACH_RE.match(s):
            return s
        return alias.get(normalise_cellline_name(s))

    out[out_col] = out[name_col].map(_resolve)
    return out


def resolution_report(resolved: pd.Series) -> dict[str, float]:
    """Summarise how much of a resolved id column actually resolved."""
    n = len(resolved)
    n_ok = int(resolved.notna().sum())
    return {
        "n": n,
        "resolved": n_ok,
        "unresolved": n - n_ok,
        "resolved_frac": (n_ok / n) if n else 0.0,
    }
