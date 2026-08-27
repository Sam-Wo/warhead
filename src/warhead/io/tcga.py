"""TCGA / GDC copy-number loader for the collateral-lethality scan (G2c).

Access: portal.gdc.cancer.gov / cBioPortal (projects COAD, READ, LIHC). We need,
per gene, the hemizygous-loss frequency in each indication - the recurrence half
of the collateral-lethality argument (the dependency half comes from DepMap).

Feeds the tidy schema the cascade expects:
    gene, indication, loss_frequency, co_deleted
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from ..config import DATA_RAW
from .base import RawMissingError

TCGA_RAW = DATA_RAW / "tcga"

# GISTIC2 thresholded call for hemizygous loss.
_GISTIC_LOSS = -1

# TCGA project -> WARHEAD indication.
_PROJECT_INDICATION = {"COAD": "CRC", "READ": "CRC", "LIHC": "HCC"}


def load_loss_recurrence(raw_dir: Path | str = TCGA_RAW) -> pd.DataFrame:
    """Per-gene hemizygous-loss frequency by indication.

    Reads a normalised intermediate if present (``tcga_loss_recurrence.csv`` with
    columns [gene, indication, loss_frequency, co_deleted]). Otherwise computes it
    from a GISTIC2 thresholded matrix ``gistic2_all_thresholded.by_genes.csv``
    (genes x samples, integer calls) plus a ``sample_project.csv`` mapping sample
    -> TCGA project, counting calls <= -1 as hemizygous loss.
    """
    raw_dir = Path(raw_dir)

    norm = raw_dir / "tcga_loss_recurrence.csv"
    if norm.exists():
        df = pd.read_csv(norm)
        need = {"gene", "indication", "loss_frequency"}
        missing = need - set(df.columns)
        if missing:
            raise ValueError(f"{norm} missing columns {sorted(missing)}")
        if "co_deleted" not in df.columns:
            df["co_deleted"] = None
        return df

    gistic = raw_dir / "gistic2_all_thresholded.by_genes.csv"
    mapping = raw_dir / "sample_project.csv"
    if not (gistic.exists() and mapping.exists()):
        raise RawMissingError(
            norm, "TCGA COAD/READ/LIHC",
            "Provide 'tcga_loss_recurrence.csv' [gene, indication, loss_frequency, "
            "co_deleted], or GISTIC2 'gistic2_all_thresholded.by_genes.csv' + "
            "'sample_project.csv' from GDC/cBioPortal",
        )

    calls = pd.read_csv(gistic, index_col=0)
    smap = pd.read_csv(mapping).set_index("sample")["project"].map(_PROJECT_INDICATION)
    rows = []
    for indication in sorted(set(smap.dropna())):
        samples = smap[smap == indication].index.intersection(calls.columns)
        if not len(samples):
            continue
        sub = calls[samples]
        freq = (sub <= _GISTIC_LOSS).mean(axis=1)
        for gene, f in freq.items():
            rows.append({"gene": gene, "indication": indication,
                         "loss_frequency": float(f), "co_deleted": None})
    return pd.DataFrame(rows)
