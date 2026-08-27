"""GDSC1/GDSC2 (raw dose-response points) loader. Access: cancerrxgene.org (do not pool GDSC1+GDSC2 naively).

Tidy output schema: compound_id, ModelID, dose_M, viability, dataset

Deferred until the raw source is local; see WARHEAD.md sec 2 for priority.
"""
from __future__ import annotations

import pandas as pd


def load_doseresponse(*args, **kwargs) -> pd.DataFrame:
    raise NotImplementedError(
        "GDSC1/GDSC2 (raw dose-response points) loader not yet wired. Source/access: cancerrxgene.org (do not pool GDSC1+GDSC2 naively). "
        "Target schema: compound_id, ModelID, dose_M, viability, dataset"
    )
