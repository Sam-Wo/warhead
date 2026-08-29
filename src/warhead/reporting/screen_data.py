"""Load the canonical per-screen frames from data/interim and attach coverage
metadata. Shared by the dashboard / summary / curve generation scripts so they
all see the same frames, counts, and top-N selection.

The dose-response screens (GDSC2 / PRISM / CTRP) share the source-agnostic
canonical schema; PDXE and NCI-60 are loaded on demand by their own callers.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from warhead.analysis.screen_meta import screen_metadata

# canonical short label -> full source name used in the metadata table
DR_SOURCES = {
    "GDSC2": "GDSC2",
    "PRISM": "PRISM Repurposing (secondary)",
    "CTRP v2": "CTRP v2",
}

_GDSC_TCGA = {"COREAD": "CRC", "LIHC": "HCC"}


def _gdsc_canonical(interim: Path) -> pd.DataFrame:
    g = pd.read_pickle(interim / "gdsc2_ec90.pkl")
    return pd.DataFrame({
        "source": "GDSC2", "compound": g.drug_name, "target": g.target, "moa": g.pathway,
        "model_id": g.cell_line, "indication": g.tcga_desc.map(_GDSC_TCGA).fillna("other"),
        "ic50_nM": g.ic50_uM * 1e3, "ec90_nM": g.ec90_uM * 1e3, "emax": np.nan,
        "ec90_extrapolated": g.ec90_range.eq("extrapolated"), "clinical_phase": pd.NA})


def load_dr_screens(interim: str | Path = "data/interim"):
    """Return (cans, counts, meta):

    - cans:   {full_source_name: canonical DataFrame} for GDSC2 / PRISM / CTRP
    - counts: {full_source_name: {compounds, lines, CRC, HCC}}
    - meta:   screen_metadata() DataFrame (all five screens, coverage + assay design)
    """
    interim = Path(interim)
    cans = {
        "GDSC2": _gdsc_canonical(interim),
        "PRISM Repurposing (secondary)": pd.read_pickle(interim / "prism_canonical.pkl"),
        "CTRP v2": pd.read_pickle(interim / "ctrp_canonical.pkl"),
    }
    counts = {s: {"compounds": int(c.compound.nunique()), "lines": int(c.model_id.nunique()),
                  "CRC": int(c[c.indication == "CRC"].model_id.nunique()),
                  "HCC": int(c[c.indication == "HCC"].model_id.nunique())}
              for s, c in cans.items()}
    return cans, counts, screen_metadata(counts)
