"""One loader per source, all returning tidy frames.

Contract (WARHEAD.md sec 4): compound identity resolved to InChIKey and cell-line
identity resolved to DepMap ModelID at ingest. Loaders read from ``data/raw`` and
raise an actionable error naming the expected file if it is absent; they never
fabricate data. The synthetic fixture (warhead.fixtures) emits frames in these
exact schemas so the downstream cascade is source-agnostic.
"""
from .base import RawMissingError, download_file, sha256_of  # noqa: F401
