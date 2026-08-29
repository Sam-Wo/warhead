"""Generic screen potency/selectivity + PRISM EC90 recovery."""
import numpy as np
import pandas as pd

from warhead.analysis.screen_potency import rank_potency, selectivity
from warhead.io.prism import to_canonical


def _canon(rows):
    cols = ["source", "compound", "target", "moa", "model_id", "indication",
            "ic50_nM", "ec90_nM", "emax", "ec90_extrapolated", "clinical_phase"]
    return pd.DataFrame(rows, columns=cols)


def test_rank_emax_filter_drops_non_killers():
    rows = []
    for i in range(8):
        rows.append(["S", "killer", "T", "", f"M{i}", "CRC", 10.0, 90.0, 0.05, False, "Launched"])
        rows.append(["S", "nonkiller", "T2", "", f"M{i}", "CRC", 5.0, 45.0, 0.95, False, "Phase 1"])
    df = _canon(rows)
    r = rank_potency(df, "CRC", min_lines=3, emax_max=0.5)
    assert "killer" in set(r["compound"])
    assert "nonkiller" not in set(r["compound"])   # emax 0.95 -> not a real kill


def test_selectivity_direction():
    rng = np.random.default_rng(0)
    rows = []
    for i in range(10):
        rows.append(["S", "Sel", "MEK", "", f"C{i}", "CRC", 20 * np.exp(rng.normal(0, .2)), 200, 0.1, False, "Launched"])
    for i in range(40):
        rows.append(["S", "Sel", "MEK", "", f"O{i}", "other", 2000 * np.exp(rng.normal(0, .2)), 20000, 0.1, False, "Launched"])
    df = _canon(rows)
    s = selectivity(df, "CRC").set_index("compound")
    assert s.loc["Sel", "delta_potency"] > 1.0
    assert bool(s.loc["Sel", "selective_potent"])


def test_prism_ec90_from_params():
    # viability(d) = lower + (upper-lower)/(1+(d/ec50)^slope); EC90 = ec50*9^(1/slope)
    params = pd.DataFrame([{
        "depmap_id": "ACH-1", "ccle_name": "X", "screen_id": "HTS002",
        "upper_limit": 1.0, "lower_limit": 0.05, "slope": 1.0, "r2": .95, "auc": .3,
        "ec50": 0.01, "ic50": 0.012, "name": "drugA", "moa": "m", "target": "T", "phase": "Launched",
    }])
    can = to_canonical(params, {"ACH-1": "HCC"})
    assert len(can) == 1
    row = can.iloc[0]
    assert abs(row["ec90_nM"] - 0.01 * 9 ** (1 / 1.0) * 1e3) < 1e-6   # ec50*9 -> nM
    assert row["indication"] == "HCC"
    assert abs(row["emax"] - 0.05) < 1e-9
