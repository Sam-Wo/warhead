"""NCI-60 / DTP (GI50/TGI/LC50, sub-nM reach) loader. Access: dtp.cancer.gov; CellMiner.

Tidy output schema: compound_id, cell_line, gi50_M, tgi_M, lc50_M

Deferred until the raw source is local; see WARHEAD.md sec 2 for priority.
"""
from __future__ import annotations

import pandas as pd


def load_gi50(*args, **kwargs) -> pd.DataFrame:
    raise NotImplementedError(
        "NCI-60 / DTP (GI50/TGI/LC50, sub-nM reach) loader not yet wired. Source/access: dtp.cancer.gov; CellMiner. "
        "Target schema: compound_id, cell_line, gi50_M, tgi_M, lc50_M"
    )
