"""Tahoe-100M (single-cell perturbation, 24h) loader. Access: HuggingFace tahoebio/Tahoe-100M.

Tidy output schema: compound_id, ModelID, gene, pseudobulk_logfc (pseudobulk first; stream for single-cell)

Deferred until the raw source is local; see WARHEAD.md sec 2 for priority.
"""
from __future__ import annotations

import pandas as pd


def load_pseudobulk(*args, **kwargs) -> pd.DataFrame:
    raise NotImplementedError(
        "Tahoe-100M (single-cell perturbation, 24h) loader not yet wired. Source/access: HuggingFace tahoebio/Tahoe-100M. "
        "Target schema: compound_id, ModelID, gene, pseudobulk_logfc (pseudobulk first; stream for single-cell)"
    )
