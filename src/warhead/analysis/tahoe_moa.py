"""Mechanism-of-action distance from a Top1-inhibitor anchor, using Tahoe-100M
transcriptional signatures (via the Harmonizome up/down gene-set libraries).

Each Harmonizome term is `drug_cellline`; we build a per-drug CONSENSUS signature
(genes moved in the same direction across a fraction of that drug's cell lines) and
score every drug's directional concordance with the Top1i anchor (Topotecan +/-
Irinotecan). Low concordance = the drug engages a DIFFERENT transcriptional program
= orthogonal MOA = a mechanistically non-redundant exatecan partner.
"""
from __future__ import annotations

import gzip
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd


def _norm(s) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def _parse_gmt(path, min_frac=0.34) -> dict:
    """drug -> consensus gene set (genes present in >= min_frac of the drug's lines)."""
    lines = {}
    with gzip.open(path, "rt") as f:
        for row in f:
            parts = row.rstrip("\n").split("\t")
            drug = parts[0].rsplit("_", 1)[0]                 # strip the _cellline suffix
            lines.setdefault(drug, []).append(set(parts[2:]))
    cons = {}
    for drug, sets in lines.items():
        c = Counter()
        for s in sets:
            c.update(s)
        thr = max(2, min_frac * len(sets))
        cons[drug] = {g for g, ct in c.items() if ct >= thr}
    return cons


def load_signatures(up_gz, dn_gz, *, min_frac=0.34):
    return _parse_gmt(up_gz, min_frac), _parse_gmt(dn_gz, min_frac)


def _concordance(up_c, dn_c, up_a, dn_a) -> float:
    """Directional signature agreement in [-1, 1]: +1 same MOA, ~0 orthogonal,
    <0 opposing. |up&up|+|dn&dn| - |up&dn|-|dn&up|, normalised by the anchor size."""
    denom = len(up_a) + len(dn_a)
    if not denom or not (len(up_c) + len(dn_c)):
        return np.nan
    agree = len(up_c & up_a) + len(dn_c & dn_a)
    oppose = len(up_c & dn_a) + len(dn_c & up_a)
    return (agree - oppose) / denom


def moa_distance_table(up_gz, dn_gz, *, anchors=("Topotecan (hydrochloride)", "Irinotecan"),
                       min_frac=0.34) -> pd.DataFrame:
    """One row per Tahoe drug: concordance with the Top1i anchor and moa_distance
    (1 - concordance, clipped to >=0; higher = more orthogonal)."""
    up, dn = load_signatures(up_gz, dn_gz, min_frac=min_frac)
    anorm = {_norm(a) for a in anchors}
    up_a = set().union(*[up[d] for d in up if _norm(d) in anorm] or [set()])
    dn_a = set().union(*[dn[d] for d in dn if _norm(d) in anorm] or [set()])
    if not up_a and not dn_a:
        raise ValueError("no Top1i anchor found in the Tahoe signature libraries")
    rows = []
    for d in up.keys() | dn.keys():
        conc = _concordance(up.get(d, set()), dn.get(d, set()), up_a, dn_a)
        rows.append({"tahoe_drug": d, "concordance": conc,
                     "moa_distance": (1 - conc) if np.isfinite(conc) else np.nan,
                     "n_sig": len(up.get(d, set())) + len(dn.get(d, set()))})
    return pd.DataFrame(rows).sort_values("moa_distance").reset_index(drop=True)
