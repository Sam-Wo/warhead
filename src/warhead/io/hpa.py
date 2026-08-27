"""Human Protein Atlas (protein/single-cell) loader. Access: proteinatlas.org.

Tidy output schema: gene, tissue, cell_type, level

Deferred until the raw source is local; see WARHEAD.md sec 2 for priority.
"""
from __future__ import annotations

import pandas as pd


def load_tissue_expression(*args, **kwargs) -> pd.DataFrame:
    raise NotImplementedError(
        "Human Protein Atlas (protein/single-cell) loader not yet wired. Source/access: proteinatlas.org. "
        "Target schema: gene, tissue, cell_type, level"
    )
