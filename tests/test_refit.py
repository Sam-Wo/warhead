"""4PL refit: parameter recovery and interval-censoring behaviour."""
import numpy as np

from warhead.curves.refit import four_pl, refit_curve, refit_frame
import pandas as pd


def _simulate(log10_ic50, emax, hill, doses, noise=0.02, seed=0):
    rng = np.random.default_rng(seed)
    v = four_pl(doses, top=1.0, emax=emax, log10_ic50=log10_ic50, hill=hill)
    return np.clip(v + rng.normal(0, noise, doses.size), 0, 1.05)


def test_recovers_known_params():
    doses = np.logspace(-10, -5, 10)
    v = _simulate(log10_ic50=-8.5, emax=0.08, hill=1.2, doses=doses, noise=0.02)
    fit = refit_curve(doses, v)
    assert fit.converged
    assert abs(fit.log10_ic50_M - (-8.5)) < 0.3     # within ~0.3 log
    assert abs(fit.emax - 0.08) < 0.06
    assert abs(fit.hill - 1.2) < 0.5
    assert fit.ic50_censoring == "none"


def test_right_censored_when_no_response():
    # IC50 two logs above the top tested dose: viability stays near 1 everywhere.
    doses = np.logspace(-11, -7, 8)
    v = _simulate(log10_ic50=-5.0, emax=0.05, hill=1.0, doses=doses, noise=0.01)
    fit = refit_curve(doses, v)
    assert fit.ic50_censoring == "right"


def test_left_censored_when_full_response_at_lowest_dose():
    # IC50 below the lowest tested dose: already killed at the assay floor.
    doses = np.logspace(-9, -5, 8)
    v = _simulate(log10_ic50=-10.5, emax=0.03, hill=1.1, doses=doses, noise=0.01)
    fit = refit_curve(doses, v)
    assert fit.ic50_censoring == "left"
    # And it must NOT be clamped up to the lowest tested dose (the classic error).
    assert fit.log10_ic50_M < -9.0


def test_refit_frame_shapes():
    doses = np.logspace(-10, -5, 8)
    rows = []
    for comp in ("A", "B"):
        for model in ("ACH-000001", "ACH-000002"):
            v = _simulate(-8.0, 0.1, 1.0, doses, seed=hash((comp, model)) % 100)
            for d, vi in zip(doses, v):
                rows.append({"compound_id": comp, "ModelID": model, "dose_M": d, "viability": vi})
    fits = refit_frame(pd.DataFrame(rows))
    assert len(fits) == 4
    assert {"ic50_M", "emax", "hill", "ic50_censoring"}.issubset(fits.columns)
