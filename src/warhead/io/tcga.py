"""TCGA/GDC (COAD, READ, LIHC): GISTIC2 CN, expr loader. Access: portal.gdc.cancer.gov; cBioPortal.

Tidy output schema: sample_id, gene, gistic2_cn, project

Deferred until the raw source is local; see WARHEAD.md sec 2 for priority.
"""
from __future__ import annotations

import pandas as pd


def load_copy_number(*args, **kwargs) -> pd.DataFrame:
    raise NotImplementedError(
        "TCGA/GDC (COAD, READ, LIHC): GISTIC2 CN, expr loader not yet wired. Source/access: portal.gdc.cancer.gov; cBioPortal. "
        "Target schema: sample_id, gene, gistic2_cn, project"
    )
