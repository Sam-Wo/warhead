"""G5 - Conjugatability. (Deferred: build order step 7.)

SMARTS match for a primary/secondary amine, hydroxyl, thiol or carboxylic acid
at a position that published SAR says tolerates substitution. A handle in the
pharmacophore is NOT a handle - cross-reference ChEMBL SAR for the scaffold.
"""
from __future__ import annotations

from ..config import load_gates


def find_handles(smiles: str, *, config: dict | None = None) -> dict[str, bool]:
    """Return which conjugation handles are present. Requires RDKit."""
    cfg = (config or load_gates())["g5"]["handles_smarts"]
    try:
        from rdkit import Chem
    except Exception as exc:  # pragma: no cover
        raise NotImplementedError("RDKit required for SMARTS handle matching") from exc
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {name: False for name in cfg}
    return {
        name: mol.HasSubstructMatch(Chem.MolFromSmarts(smarts))
        for name, smarts in cfg.items()
    }
