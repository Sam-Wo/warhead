"""Stream the DepMap 24Q2 protein-coding expression matrix and extract a small set of
named gene columns (never storing the full 461 MB). For the exatecan-partner axes:
SLFN11 (Top1i sensitiser), TOP1/TOP2A (targets), and the efflux transporters.

    PYTHONPATH=src py scripts/fetch_depmap_genes.py
-> data/interim/depmap_genes.csv  (ModelID + one column per gene, TPM log2p1)
"""
import csv
import io
import urllib.request
from pathlib import Path

URL = "https://ndownloader.figshare.com/files/46490878"
GENES = ["SLFN11", "TOP1", "TOP2A", "ABCB1", "ABCG2", "ABCC1", "SLFN12", "SCHLAP1"]
OUT = Path("data/interim/depmap_genes.csv")

req = urllib.request.Request(URL, headers={"Accept-Encoding": "identity"})
with urllib.request.urlopen(req, timeout=120) as resp:
    reader = csv.reader(io.TextIOWrapper(resp, encoding="utf-8", newline=""))
    header = next(reader)
    idx = {}
    for i, col in enumerate(header):
        sym = col.split(" (")[0]
        if sym in GENES and sym not in idx:
            idx[sym] = i
    genes = [g for g in GENES if g in idx]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with OUT.open("w", newline="") as fo:
        w = csv.writer(fo)
        w.writerow(["ModelID"] + genes)
        for row in reader:
            if not row or not row[0]:
                continue
            w.writerow([row[0]] + [row[idx[g]] for g in genes])
            n += 1
print(f"wrote {OUT} ({n} models; genes {genes})")
