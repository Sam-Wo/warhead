"""ChEMBL (bioactivity/SAR/structures) loader. Access: ebi.ac.uk/chembl.

Tidy output schema: molecule_chembl_id, inchikey, target, standard_type, standard_value_nM

Deferred until the raw source is local; see WARHEAD.md sec 2 for priority.
"""
from __future__ import annotations

import pandas as pd


def load_bioactivity(*args, **kwargs) -> pd.DataFrame:
    raise NotImplementedError(
        "ChEMBL (bioactivity/SAR/structures) loader not yet wired. Source/access: ebi.ac.uk/chembl. "
        "Target schema: molecule_chembl_id, inchikey, target, standard_type, standard_value_nM"
    )
