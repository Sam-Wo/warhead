"""JUMP Cell Painting (morphology) loader. Access: cellpainting-gallery on AWS Open Data (normalise by dose).

Tidy output schema: compound_id, feature, value

Deferred until the raw source is local; see WARHEAD.md sec 2 for priority.
"""
from __future__ import annotations

import pandas as pd


def load_profiles(*args, **kwargs) -> pd.DataFrame:
    raise NotImplementedError(
        "JUMP Cell Painting (morphology) loader not yet wired. Source/access: cellpainting-gallery on AWS Open Data (normalise by dose). "
        "Target schema: compound_id, feature, value"
    )
