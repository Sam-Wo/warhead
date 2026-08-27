"""Common gate contract.

WARHEAD.md sec 4 Conventions: "Every gate emits (passed_df, failed_df,
reason_column). Nothing is silently dropped." ``GateResult`` enforces that shape
and carries the config snapshot + summary for the per-gate audit trail.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import pandas as pd


@dataclass
class GateResult:
    gate: str                       # e.g. "G2b"
    passed: pd.DataFrame
    failed: pd.DataFrame
    reason_col: str                 # column in `failed` explaining the failure
    config: dict[str, Any] = field(default_factory=dict)
    summary: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.reason_col not in self.failed.columns and len(self.failed):
            raise ValueError(
                f"{self.gate}: failed frame missing reason column '{self.reason_col}'"
            )

    @property
    def n_in(self) -> int:
        return len(self.passed) + len(self.failed)

    @property
    def n_pass(self) -> int:
        return len(self.passed)

    def provenance(self) -> dict[str, Any]:
        return {
            "gate": self.gate,
            "n_in": self.n_in,
            "n_pass": self.n_pass,
            "n_fail": len(self.failed),
            "config": self.config,
            "summary": self.summary,
        }
