"""Shared, session-scoped fixtures.

Refitting the synthetic dose-response is the expensive step, so we compute the
G2b and G3b cascade states once per test session and share them. Tests that need
the refit frame (e.g. G2a) reuse ``g2b_state.fits`` rather than refitting again.
"""
import pytest

from warhead.cascade import run_g2b_slice, run_g3b_slice
from warhead.fixtures import make

TEST_SEED = 7
TEST_N = 80


@pytest.fixture(scope="session")
def synth():
    return make(seed=TEST_SEED, n_lines=TEST_N)


@pytest.fixture(scope="session")
def g2b_state(synth):
    # apply_g1_filter=False so the planted controls that fail G1 on purpose are
    # still visible to the G2b assertions.
    return run_g2b_slice(synth.dose_response, synth.model_meta, apply_g1_filter=False)


@pytest.fixture(scope="session")
def g3b_state(synth):
    return run_g3b_slice(synth.dose_response, synth.expression)
