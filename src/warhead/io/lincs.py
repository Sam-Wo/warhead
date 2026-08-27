"""LINCS L1000 (level-5 signatures) loader. Access: clue.io.

Tidy output schema: sig_id, compound_id, ModelID, gene(978 landmarks), zscore

Deferred until the raw source is local; see WARHEAD.md sec 2 for priority.
"""
from __future__ import annotations

import pandas as pd


def load_signatures(*args, **kwargs) -> pd.DataFrame:
    raise NotImplementedError(
        "LINCS L1000 (level-5 signatures) loader not yet wired. Source/access: clue.io. "
        "Target schema: sig_id, compound_id, ModelID, gene(978 landmarks), zscore"
    )
