"""ADCdb (labelled positive set of ADC payloads) loader. Access: adcdb.idrblab.net.

Tidy output schema: adc_id, payload_name, payload_class, smiles, inchikey, status

Deferred until the raw source is local; see WARHEAD.md sec 2 for priority.
"""
from __future__ import annotations

import pandas as pd


def load_payloads(*args, **kwargs) -> pd.DataFrame:
    raise NotImplementedError(
        "ADCdb (labelled positive set of ADC payloads) loader not yet wired. Source/access: adcdb.idrblab.net. "
        "Target schema: adc_id, payload_name, payload_class, smiles, inchikey, status"
    )
