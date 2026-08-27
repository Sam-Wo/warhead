"""GTEx (bulk normal-tissue expression) loader. Access: gtexportal.org.

Tidy output schema: gene, tissue, median_tpm

Deferred until the raw source is local; see WARHEAD.md sec 2 for priority.
"""
from __future__ import annotations

import pandas as pd


def load_median_expression(*args, **kwargs) -> pd.DataFrame:
    raise NotImplementedError(
        "GTEx (bulk normal-tissue expression) loader not yet wired. Source/access: gtexportal.org. "
        "Target schema: gene, tissue, median_tpm"
    )
