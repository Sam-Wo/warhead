"""Download the Human Protein Atlas RNA consensus tissue table (normal-tissue nTPM
per gene x 51 tissues) for the G6 therapeutic-window gate. ~5 MB zip.

    PYTHONPATH=src py scripts/fetch_hpa_consensus.py
-> data/interim/hpa_consensus.tsv  (Gene, Gene name, Tissue, nTPM)
"""
import io
import urllib.request
import zipfile
from pathlib import Path

import pandas as pd

URL = "https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip"
OUT = Path("data/interim/hpa_consensus.tsv")

with urllib.request.urlopen(URL, timeout=180) as r:
    z = zipfile.ZipFile(io.BytesIO(r.read()))
df = pd.read_csv(z.open(z.namelist()[0]), sep="\t")
OUT.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(OUT, sep="\t", index=False)
print(f"wrote {OUT} ({df.shape[0]} rows, {df['Tissue'].nunique()} tissues, "
      f"{df['Gene name'].nunique()} genes)")
