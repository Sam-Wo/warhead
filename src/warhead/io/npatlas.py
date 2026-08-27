"""NPAtlas (microbial natural products) loader. Access: npatlas.org.

Tidy output schema: npaid, inchikey, smiles

Deferred until the raw source is local; see WARHEAD.md sec 2 for priority.
"""
from __future__ import annotations

import pandas as pd


def load_structures(*args, **kwargs) -> pd.DataFrame:
    raise NotImplementedError(
        "NPAtlas (microbial natural products) loader not yet wired. Source/access: npatlas.org. "
        "Target schema: npaid, inchikey, smiles"
    )
