"""GDSC EC90 recovery and tissue-selectivity mechanics (small synthetic frames)."""
import numpy as np
import pandas as pd

from warhead.io.gdsc import add_ec90, map_to_depmap, _mean_viability
from warhead.analysis.gdsc_ec90 import selectivity


def _synth_curve(ic50_uM, scal, lo=1e-3, hi=10.0):
    """One GDSC-shaped row whose AUC is the exact mean viability of the model."""
    xmid = np.log(ic50_uM)
    auc = _mean_viability(scal, xmid, np.log(lo), np.log(hi))
    return {"min_conc_uM": lo, "max_conc_uM": hi, "ln_ic50": xmid,
            "auc": auc, "ic50_uM": ic50_uM}


def test_ec90_recovers_scal_and_ec90():
    # IC50 off-centre (near the low edge) so scal is identifiable from AUC.
    scal_true = 0.8
    df = pd.DataFrame([_synth_curve(0.005, scal_true)])
    out = add_ec90(df)
    assert abs(out["scal"].iloc[0] - scal_true) < 1e-3
    assert abs(out["ec90_uM"].iloc[0] - 0.005 * 9 ** scal_true) < 1e-3
    assert out["ec90_confidence"].iloc[0] == "ok"


def test_ec90_flags_centred_ic50_as_low_confidence():
    # IC50 at the centre of the range -> AUC ~ 0.5 for any slope -> unidentifiable.
    df = pd.DataFrame([_synth_curve(0.1, 0.8, lo=1e-3, hi=10.0)])  # ln(0.1) centred
    out = add_ec90(df)
    assert out["ec90_confidence"].iloc[0] == "low"


def test_ec90_greater_than_ic50():
    df = pd.DataFrame([_synth_curve(0.05, 1.2)])
    out = add_ec90(df)
    assert out["ec90_uM"].iloc[0] > out["ic50_uM"].iloc[0]


def test_selectivity_flags_indication_potent_compound():
    rng = np.random.default_rng(0)
    rows = []
    # "SelDrug": 100x more potent (lower IC50) in CRC (COREAD) than elsewhere.
    for i in range(12):
        rows.append({"drug_name": "SelDrug", "target": "X", "pathway": "P",
                     "tcga_desc": "COREAD", "ic50_uM": 0.01 * np.exp(rng.normal(0, .2)),
                     "ec90_uM": 0.1})
    for i in range(60):
        rows.append({"drug_name": "SelDrug", "target": "X", "pathway": "P",
                     "tcga_desc": "OTHER", "ic50_uM": 1.0 * np.exp(rng.normal(0, .2)),
                     "ec90_uM": 10.0})
    # "FlatDrug": same potency everywhere.
    for t, n in (("COREAD", 12), ("OTHER", 60)):
        for i in range(n):
            rows.append({"drug_name": "FlatDrug", "target": "Y", "pathway": "Q",
                         "tcga_desc": t, "ic50_uM": 0.5 * np.exp(rng.normal(0, .2)),
                         "ec90_uM": 5.0})
    df = pd.DataFrame(rows)
    sel = selectivity(df, "CRC").set_index("drug_name")
    assert sel.loc["SelDrug", "delta_potency"] > 1.0     # ~2 logs more potent
    assert bool(sel.loc["SelDrug", "selective_potent"])
    assert not bool(sel.loc["FlatDrug", "selective"])


def test_map_to_depmap_joins_on_sanger_id():
    fitted = pd.DataFrame({
        "drug_name": ["D1", "D1"], "sanger_model_id": ["SIDM1", "SIDM2"],
        "ln_ic50": [0.0, 1.0], "auc": [.5, .6],
    })
    meta = pd.DataFrame({
        "ModelID": ["ACH-1", "ACH-2"], "SangerModelID": ["SIDM1", "SIDM2"],
        "doubling_time_hours": [24.0, 48.0], "OncotreeCode": ["LIHC", "COAD"],
    })
    out = map_to_depmap(fitted, meta)
    assert set(out["ModelID"]) == {"ACH-1", "ACH-2"}
    assert "doubling_time_hours" in out.columns
    # a GDSC line with no DepMap match is dropped (inner join), not fabricated
    fitted2 = pd.concat([fitted, pd.DataFrame([{
        "drug_name": "D1", "sanger_model_id": "SIDM_UNKNOWN", "ln_ic50": 2.0, "auc": .7}])])
    assert len(map_to_depmap(fitted2, meta)) == 2
