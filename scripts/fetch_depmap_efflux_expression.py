"""Stream the DepMap 24Q2 protein-coding expression matrix and extract ONLY the two
efflux-transporter columns (ABCB1, ABCG2), so we never store the full 461 MB.

Source: figshare mirror of DepMap 24Q2, OmicsExpressionProteinCodingGenesTPMLogp1.csv
(id 46490878). The DepMap portal is Cloudflare-gated; this mirror is the sanctioned
route. Output: data/interim/depmap_abcb1_abcg2.csv  (ModelID, ABCB1, ABCG2 - TPM log2p1).

    PYTHONPATH=src py scripts/fetch_depmap_efflux_expression.py
"""
import csv
import io
import urllib.request
from pathlib import Path

URL = "https://ndownloader.figshare.com/files/46490878"
GENES = {"ABCB1": None, "ABCG2": None}          # column index filled from the header
OUT = Path("data/interim/depmap_abcb1_abcg2.csv")

req = urllib.request.Request(URL, headers={"Accept-Encoding": "identity"})
with urllib.request.urlopen(req, timeout=120) as resp:
    reader = csv.reader(io.TextIOWrapper(resp, encoding="utf-8", newline=""))
    header = next(reader)
    for i, col in enumerate(header):
        sym = col.split(" (")[0]
        if sym in GENES:
            GENES[sym] = i
    if any(v is None for v in GENES.values()):
        raise SystemExit(f"columns not found: {GENES}")
    i_abcb1, i_abcg2 = GENES["ABCB1"], GENES["ABCG2"]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with OUT.open("w", newline="") as fo:
        w = csv.writer(fo)
        w.writerow(["ModelID", "ABCB1", "ABCG2"])
        for row in reader:
            if not row or not row[0]:
                continue
            w.writerow([row[0], row[i_abcb1], row[i_abcg2]])
            n += 1
print(f"wrote {OUT} ({n} models; ABCB1 col {i_abcb1}, ABCG2 col {i_abcg2})")
