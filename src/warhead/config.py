"""Project paths and the single source of truth for gate thresholds.

All arguable numbers live in ``config/gates.yaml``. Nothing in the cascade
hard-codes a threshold; everything reads it from here so a reviewer can move a
number and re-run.
"""
from __future__ import annotations

import os
from functools import lru_cache
from pathlib import Path
from typing import Any

import yaml

# Repo root = three parents up from this file (src/warhead/config.py).
REPO_ROOT = Path(__file__).resolve().parents[2]

DATA_RAW = REPO_ROOT / "data" / "raw"
DATA_INTERIM = REPO_ROOT / "data" / "interim"
DATA_PROCESSED = REPO_ROOT / "data" / "processed"
REPORTS = REPO_ROOT / "reports"
CONFIG_DIR = REPO_ROOT / "config"
GATES_YAML = CONFIG_DIR / "gates.yaml"


def _resolve_gates_path(path: str | os.PathLike[str] | None) -> Path:
    if path is not None:
        return Path(path)
    env = os.environ.get("WARHEAD_GATES")
    return Path(env) if env else GATES_YAML


@lru_cache(maxsize=8)
def _load_cached(path_str: str) -> dict[str, Any]:
    with open(path_str, "r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def load_gates(path: str | os.PathLike[str] | None = None) -> dict[str, Any]:
    """Load ``gates.yaml`` (cached). Set ``WARHEAD_GATES`` to override for a
    threshold sensitivity sweep."""
    return _load_cached(str(_resolve_gates_path(path).resolve()))


def ensure_dirs() -> None:
    for d in (DATA_RAW, DATA_INTERIM, DATA_PROCESSED, REPORTS):
        d.mkdir(parents=True, exist_ok=True)
