"""G2b end to end on the synthetic fixture: the pipeline must re-detect the
planted proliferation dependence and separate mitotic from non-mitotic classes.
"""
import numpy as np
import pandas as pd

from warhead.cascade import run_g2b_slice
from warhead.fixtures import make
from warhead.gates.g2_delivery import (
    collateral_lethality_scan,
    gate_g2a,
    sensitivity_from_fits,
)
from warhead.curves.qc import qc_flags
from warhead.curves.refit import refit_frame


def test_recovers_planted_prolif_beta():
    data = make(seed=7, n_lines=60)
    state = run_g2b_slice(data.dose_response, data.model_meta, apply_g1_filter=False)
    st = state.g2b_stats.set_index("compound_id")
    truth = data.compound_truth.set_index("compound_id")["prolif_beta"]
    common = st.index.intersection(truth.index)
    corr = pd.Series({c: st.loc[c, "slope"] for c in common}).corr(truth.loc[common])
    assert corr > 0.9, f"slope recovery corr too low: {corr:.3f}"


def test_mitotic_fail_nonmitotic_pass():
    data = make(seed=7, n_lines=60)
    state = run_g2b_slice(data.dose_response, data.model_meta, apply_g1_filter=False)
    passed = set(state.g2b.passed["compound_id"])
    st = state.g2b_stats.set_index("compound_id")

    # Antimitotics: strong positive slope, excluded from the keep set.
    assert st.loc["MMAE_like", "std_slope"] > 0.5
    for mitotic in ("MMAE_like", "DM1_like", "exatecan_like", "SN38_like"):
        assert mitotic not in passed, f"{mitotic} should fail G2b"

    # Non-mitotic payload classes are proliferation-independent -> kept
    # (allow one to wobble to borderline across seeds).
    nonmitotic = {"amanitin_like", "PF846_like", "degrader_like", "thailanstatin_like"}
    assert len(nonmitotic & passed) >= 3

    # Slope ordering matches the spec: antimitotic > Top1i > non-mitotic.
    assert st.loc["MMAE_like", "std_slope"] > st.loc["exatecan_like", "std_slope"]
    assert st.loc["exatecan_like", "std_slope"] > st.loc["amanitin_like", "std_slope"]


def test_g2a_flags_efflux_substrates():
    data = make(seed=7, n_lines=60)
    fits = qc_flags(refit_frame(data.dose_response))
    # log10_ic50 (resistance axis): an efflux substrate's IC50 RISES with ABCB1,
    # so a positive slope flags it (WARHEAD.md G2a).
    sens = sensitivity_from_fits(fits, metric="log10_ic50")
    res = gate_g2a(sens, data.expression)
    failed = set(res.failed["compound_id"])
    # The planted efflux substrate must be caught by G2a...
    assert "MDR_substrate_like" in failed
    # ...and a clean non-mitotic (no ABCB1 dependence) must not be.
    assert "amanitin_like" in set(res.passed["compound_id"])


def test_collateral_recovers_polr2a():
    data = make(seed=7, n_lines=60)
    scan = collateral_lethality_scan(data.chronos, data.copy_number)
    row = scan.set_index("gene").loc["POLR2A"]
    assert row["delta"] < 0            # CN loss -> more dependent
    assert row["significant"]           # positive control recovers
    # A null housekeeping gene should not be flagged.
    assert not scan.set_index("gene").loc["GAPDH", "significant"]
