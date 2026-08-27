"""PRISM Repurposing loader (DepMap portal).

Only the SECONDARY screen is usable for G1: it carries 8-point dose-response,
whereas the primary screen is single-dose (WARHEAD.md G1 / data table). We load
the per-dose viability so the curves can be refit from raw, not read off
someone else's IC50.

Pooled barcoded format under-represents slow-growing lines; that bias is carried
downstream and corrected in G2b's weighting, not here.

Expected file: ``secondary-screen-dose-response-curve-info.csv`` plus the
treatment-level viability matrix, dropped into ``data/raw/prism/``. Output is the
tidy long schema every downstream module expects:

    compound_id, ModelID, dose_M, viability   (+ optional inchikey, replicate)
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from ..config import DATA_RAW
from .base import RawMissingError

PRISM_RAW = DATA_RAW / "prism"

TIDY_COLUMNS = ["compound_id", "ModelID", "dose_M", "viability"]


def load_secondary_doseresponse(raw_dir: Path | str = PRISM_RAW) -> pd.DataFrame:
    """Load PRISM secondary as a tidy long dose-response frame.

    Because PRISM's public layout has shifted across releases, this reads a
    normalised intermediate if present (``prism_secondary_long.parquet`` /
    ``.csv`` with the TIDY_COLUMNS), and otherwise raises with instructions. The
    normalisation from the raw matrix lives in ``notebooks/prism_ingest`` (to be
    added when the raw files are local).
    """
    raw_dir = Path(raw_dir)
    for name in ("prism_secondary_long.parquet", "prism_secondary_long.csv"):
        path = raw_dir / name
        if path.exists():
            df = pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path)
            missing = set(TIDY_COLUMNS) - set(df.columns)
            if missing:
                raise ValueError(f"{path} missing columns {sorted(missing)}")
            return df[TIDY_COLUMNS + [c for c in df.columns if c not in TIDY_COLUMNS]]
    raise RawMissingError(
        raw_dir / "prism_secondary_long.parquet",
        "PRISM secondary",
        "Download the secondary screen from depmap.org/portal/download and "
        "normalise to columns [compound_id, ModelID, dose_M, viability]",
    )
