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
from .gates.g2_delivery import (
    collateral_lethality_scan,
    gate_g2b,
    proliferation_stats,
    rank_collateral_targets,
    sensitivity_from_fits,
)
from .gates.g3_moa import run_partner_search_from_expression

# CRC = COAD + READ; HCC = LIHC (WARHEAD.md scope).
_INDICATION_ONCOTREE = {"CRC": {"COAD", "READ"}, "HCC": {"LIHC"}}


@dataclass
class CascadeState:
    fits: pd.DataFrame | None = None
    g1: Any = None
    sensitivity: pd.DataFrame | None = None
    g2b_stats: pd.DataFrame | None = None
    g2b: Any = None
    g3b: pd.DataFrame | None = None
    g2c: pd.DataFrame | None = None
    provenance: list[dict] = field(default_factory=list)

    def add(self, step: str, info: dict) -> None:
        self.provenance.append({"step": step, **info})


def _refit_g1(dose_response: pd.DataFrame, cfg: dict, state: CascadeState):
    """Shared front half: refit every curve, QC, then the G1 potency gate."""
    fits = qc_flags(refit_frame(dose_response, **cfg["g1"]["refit"]))
    state.fits = fits
    state.add("refit", summarise_qc(fits))
    g1 = gate_g1(fits, config=cfg)
    state.g1 = g1
    state.add("G1", g1.provenance())
    return fits, g1


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

    # --- G1a refit + QC + G1 potency gate.
    fits, g1 = _refit_g1(dose_response, cfg, state)

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


def run_g3b_slice(
    dose_response: pd.DataFrame,
    expression: pd.DataFrame,
    *,
    config: dict | None = None,
    apply_g1_filter: bool = True,
) -> CascadeState:
    """Refit -> QC -> G1 -> G3b orthogonal-resistance (exatecan partner) search.

    Uses -log10(IC50) as the potency axis (higher = more potent) and pulls SLFN11
    and ABCB1 from the long expression frame. Restricts to G1 passers by default:
    a partner must be potent, not merely orthogonal.
    """
    cfg = config or load_gates()
    state = CascadeState()

    fits, g1 = _refit_g1(dose_response, cfg, state)

    # Potency axis for orthogonality search (higher = more potent).
    sens = sensitivity_from_fits(fits, metric="neg_log10_ic50")
    if apply_g1_filter:
        keep = set(g1.passed["compound_id"])
        sens = sens[sens["compound_id"].isin(keep)]
    state.sensitivity = sens

    g3b = run_partner_search_from_expression(sens, expression, config=cfg)
    state.g3b = g3b
    state.add(
        "G3b",
        {
            "n_candidates": len(g3b),
            "n_partner_candidates": int(g3b["is_partner_candidate"].sum()) if len(g3b) else 0,
            "top_partner": g3b.iloc[0]["compound_id"] if len(g3b) else None,
        },
    )
    return state


def run_g2c_slice(
    chronos: pd.DataFrame,
    copy_number: pd.DataFrame,
    tcga_recurrence: pd.DataFrame,
    model_meta: pd.DataFrame,
    *,
    indication: str = "CRC",
    config: dict | None = None,
) -> CascadeState:
    """G2c collateral-lethality scan for one indication (CRC or HCC).

    Restricts DepMap lines to the indication's lineages, runs the CN-loss vs
    CN-neutral Chronos scan, then joins TCGA loss recurrence to produce the ranked
    payload-target list. Does not use curve refit - it is dependency/CN based.
    """
    cfg = config or load_gates()
    state = CascadeState()

    codes = _INDICATION_ONCOTREE.get(indication.upper())
    if codes:
        ind_models = model_meta[model_meta["OncotreeCode"].isin(codes)]["ModelID"].tolist()
    else:
        ind_models = None  # whole panel

    scan = collateral_lethality_scan(
        chronos, copy_number, indication_models=ind_models, config=cfg
    )
    targets = rank_collateral_targets(
        scan, tcga_recurrence, indication=indication.upper(), config=cfg
    )
    state.g2c = targets

    pc = cfg["g2"]["collateral"]["positive_control"]
    pc_row = targets[targets["gene"] == pc]
    state.add(
        "G2c",
        {
            "indication": indication.upper(),
            "n_genes_tested": len(targets),
            "n_candidate_targets": int(targets["candidate"].sum()),
            "positive_control": pc,
            "positive_control_recovered": bool(pc_row["candidate"].iloc[0]) if len(pc_row) else False,
        },
    )
    return state
