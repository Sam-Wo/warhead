"""Cascade orchestration + per-gate provenance.

Only the G1 -> G2b slice is wired end to end for now (build order steps 2-3).
The remaining gates have interfaces in ``warhead.gates`` and slot in here as the
data lands. Provenance is a running record so every shortlist entry can be traced
to the gate decisions and thresholds that produced it.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import pandas as pd

from .config import load_gates
from .curves.qc import qc_flags, summarise_qc
from .curves.refit import refit_frame
from .gates.g1_potency import gate_g1
from .gates.g2_delivery import gate_g2b, proliferation_stats, sensitivity_from_fits


@dataclass
class CascadeState:
    fits: pd.DataFrame | None = None
    g1: Any = None
    sensitivity: pd.DataFrame | None = None
    g2b_stats: pd.DataFrame | None = None
    g2b: Any = None
    provenance: list[dict] = field(default_factory=list)

    def add(self, step: str, info: dict) -> None:
        self.provenance.append({"step": step, **info})


def run_g2b_slice(
    dose_response: pd.DataFrame,
    model_meta: pd.DataFrame,
    *,
    config: dict | None = None,
    apply_g1_filter: bool = True,
) -> CascadeState:
    """Refit -> QC -> G1 -> G2b on a tidy dose-response frame.

    ``apply_g1_filter`` restricts G2b to compounds that clear G1 (potency is a
    precondition for the proliferation argument to matter). Set False to inspect
    G2b behaviour across all compounds, including the planted controls that fail
    G1 on purpose.
    """
    cfg = config or load_gates()
    state = CascadeState()

    # --- G1a: refit every curve, then QC.
    fits = refit_frame(dose_response, **cfg["g1"]["refit"])
    fits = qc_flags(fits)
    state.fits = fits
    state.add("refit", summarise_qc(fits))

    # --- G1 gate (potency).
    g1 = gate_g1(fits, config=cfg)
    state.g1 = g1
    state.add("G1", g1.provenance())

    # --- Build the G2b sensitivity axis.
    metric = cfg["g2"]["proliferation"]["sensitivity_metric"]
    sens = sensitivity_from_fits(fits, metric=metric)
    if apply_g1_filter:
        keep = set(g1.passed["compound_id"])
        sens = sens[sens["compound_id"].isin(keep)]
    state.sensitivity = sens

    # --- G2b.
    stats = proliferation_stats(sens, model_meta, config=cfg)
    g2b = gate_g2b(sens, model_meta, config=cfg)
    state.g2b_stats = stats
    state.g2b = g2b
    state.add("G2b", g2b.provenance())
    return state
