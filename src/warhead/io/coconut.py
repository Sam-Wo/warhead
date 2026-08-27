"""COCONUT (~400k natural products) loader. Access: coconut.naturalproducts.net.

Tidy output schema: coconut_id, inchikey, smiles

Deferred until the raw source is local; see WARHEAD.md sec 2 for priority.
"""
from __future__ import annotations

import pandas as pd


def load_structures(*args, **kwargs) -> pd.DataFrame:
    raise NotImplementedError(
        "COCONUT (~400k natural products) loader not yet wired. Source/access: coconut.naturalproducts.net. "
        "Target schema: coconut_id, inchikey, smiles"
    )
