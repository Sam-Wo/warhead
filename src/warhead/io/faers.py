"""FAERS (adverse event reports) loader. Access: FDA quarterly extract files.

Tidy output schema: case_id, drug, reaction, regimen

Deferred until the raw source is local; see WARHEAD.md sec 2 for priority.
"""
from __future__ import annotations

import pandas as pd


def load_reports(*args, **kwargs) -> pd.DataFrame:
    raise NotImplementedError(
        "FAERS (adverse event reports) loader not yet wired. Source/access: FDA quarterly extract files. "
        "Target schema: case_id, drug, reaction, regimen"
    )
