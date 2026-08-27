import numpy as np
from scipy import stats as sstats

from warhead.stats import balance_weights, benjamini_hochberg, weighted_linregress


def test_weighted_linregress_matches_scipy_unweighted():
    rng = np.random.default_rng(1)
    x = rng.normal(size=50)
    y = 2.0 * x + 1.0 + rng.normal(scale=0.3, size=50)
    f = weighted_linregress(x, y)
    sp = sstats.linregress(x, y)
    assert abs(f.slope - sp.slope) < 1e-6
    assert abs(f.p - sp.pvalue) < 1e-6


def test_bh_monotone_and_bounds():
    p = np.array([0.001, 0.01, 0.02, 0.5, np.nan])
    q = benjamini_hochberg(p)
    finite = np.isfinite(q)
    assert np.all(q[finite] >= p[finite] - 1e-12)  # q >= p
    assert np.all(q[finite] <= 1.0)
    assert np.isnan(q[-1])


def test_balance_upweights_sparse_slow_lines():
    # Fast-line-heavy sample (dense 20-60h) with a few rare slow lines (90-130h).
    fast = np.linspace(20.0, 60.0, 40)
    slow = np.array([95.0, 110.0, 130.0])
    v = np.concatenate([fast, slow])
    w = balance_weights(v, n_bins=5)
    # Equal-width binning up-weights the sparse slow bin above the dense fast one.
    assert w[-3:].mean() > w[:40].mean()
    assert abs(w.mean() - 1.0) < 1e-6
