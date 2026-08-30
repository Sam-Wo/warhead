"""Minimal PubChem name -> SMILES lookup (stdlib only), cached to data/interim.

Used by the Phase-B chemical gates (G4/G5), which need a structure. Lookup is by
compound name via PUG-REST; research-code names (e.g. CTRP internal ids) often do
not resolve and come back as None - which is itself informative for conjugatability
(no structure = cannot assess a linker handle).
"""
from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

_URL = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/SMILES/JSON")


def _query(name: str) -> tuple[int | None, str | None]:
    url = _URL.format(urllib.parse.quote(name))
    try:
        with urllib.request.urlopen(url, timeout=20) as r:
            props = json.load(r)["PropertyTable"]["Properties"][0]
            return props.get("CID"), props.get("SMILES")
    except Exception:
        return None, None


def fetch_smiles(names, cache_path="data/interim/pubchem_smiles.csv", *, pause=0.25) -> dict:
    """Return {name: SMILES or None}. Results are cached (name, cid, smiles); only
    uncached names hit the network (rate-limited by `pause` seconds)."""
    cache_path = Path(cache_path); cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache = pd.read_csv(cache_path) if cache_path.exists() else pd.DataFrame(columns=["name", "cid", "smiles"])
    known = dict(zip(cache["name"], cache["smiles"]))
    new = []
    for nm in dict.fromkeys(names):                      # de-dup, keep order
        if nm in known:
            continue
        cid, smi = _query(nm)
        new.append({"name": nm, "cid": cid, "smiles": smi})
        known[nm] = smi
        time.sleep(pause)
    if new:
        add = pd.DataFrame(new)
        cache = add if cache.empty else pd.concat([cache, add], ignore_index=True)
        cache.to_csv(cache_path, index=False)
    return {nm: (known.get(nm) if pd.notna(known.get(nm)) else None) for nm in names}
