"""Synthetic, schema-faithful data for running and testing the cascade before
real sources are local. Frames match the tidy schemas the real loaders emit, so
the same downstream code runs on synthetic and real data unchanged."""
from .synth import SyntheticData, make  # noqa: F401
