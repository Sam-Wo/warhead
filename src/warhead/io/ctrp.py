"""CTRP v2 (dense dose-response, MOA-annotated) loader. Access: CTD2 portal; mirrored on DepMap.

Tidy output schema: compound_id, ModelID, dose_M, viability, moa

Deferred until the raw source is local; see WARHEAD.md sec 2 for priority.
"""
from __future__ import annotations

import pandas as pd


def load_doseresponse(*args, **kwargs) -> pd.DataFrame:
    raise NotImplementedError(
        "CTRP v2 (dense dose-response, MOA-annotated) loader not yet wired. Source/access: CTD2 portal; mirrored on DepMap. "
        "Target schema: compound_id, ModelID, dose_M, viability, moa"
    )
