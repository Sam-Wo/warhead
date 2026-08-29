"""Real-data G2b: GDSC potency vs DepMap growth rate.

Runs WARHEAD's proliferation-independence gate on real data by regressing GDSC2
log10(IC50) against a DepMap screen-inferred growth-rate proxy for doubling time,
across the ~600 cell lines shared by the two resources (joined on SangerModelID).

This is the real-data check of G2b's premise: agents that need active
proliferation (DNA-replication / mitotic) should lose potency in slow-growing
lines (positive slope); mechanism-targeted agents should stay flat.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import DATA_RAW, load_gates
from ..gates.g2_delivery import gate_g2b, proliferation_stats
from ..io.depmap import build_growth_model_meta
from ..io.gdsc import load_fitted, map_to_depmap


@dataclass
class RealG2B:
    stats: pd.DataFrame       # per-compound slope + target/pathway
    gate: object              # GateResult
    sensitivity: pd.DataFrame # compound_id, ModelID, sensitivity(=log10 IC50)
    model_meta: pd.DataFrame  # ModelID, doubling_time_hours (proxy)
    n_lines: int


def run_real_g2b(
    gdsc_dir: Path | str = DATA_RAW / "gdsc",
    depmap_dir: Path | str = DATA_RAW / "depmap",
    *,
    dataset: str = "GDSC2",
    config: dict | None = None,
) -> RealG2B:
    cfg = config or load_gates()
    meta = build_growth_model_meta(depmap_dir)
    gd = map_to_depmap(load_fitted(gdsc_dir, dataset), meta)
    gd["sensitivity"] = gd["ln_ic50"] / np.log(10)  # log10 IC50

    sensitivity = gd.rename(columns={"drug_name": "compound_id"})[
        ["compound_id", "ModelID", "sensitivity"]
    ]
    model_meta = meta[["ModelID", "doubling_time_hours", "OncotreeCode"]].drop_duplicates("ModelID")

    stats = proliferation_stats(sensitivity, model_meta, config=cfg)
    gate = gate_g2b(sensitivity, model_meta, config=cfg)

    # annotate with mechanism for interpretation
    ann = gd.groupby("drug_name")[["target", "pathway"]].first().reset_index()
    stats = stats.merge(ann, left_on="compound_id", right_on="drug_name", how="left").drop(columns="drug_name")
    return RealG2B(stats, gate, sensitivity, model_meta, int(sensitivity["ModelID"].nunique()))
