"""Download the DepMap 24Q2 hotspot-mutation matrix (models x genes, binary), used to
test whether lineage-selectivity findings are confounded by driver genotype
(e.g. is 'CRC-selective' MEK sensitivity really KRAS/BRAF/NRAS-mutant selectivity).
Small (~4 MB). Source: figshare mirror of DepMap 24Q2 (file id 46500379).

    PYTHONPATH=src py scripts/fetch_depmap_hotspot.py
-> data/raw/depmap/OmicsSomaticMutationsMatrixHotspot.csv
"""
import os
import urllib.request
from pathlib import Path

OUT = Path("data/raw/depmap/OmicsSomaticMutationsMatrixHotspot.csv")
OUT.parent.mkdir(parents=True, exist_ok=True)
urllib.request.urlretrieve("https://ndownloader.figshare.com/files/46500379", OUT)
print(f"wrote {OUT} ({os.path.getsize(OUT) / 1e6:.1f} MB)")
