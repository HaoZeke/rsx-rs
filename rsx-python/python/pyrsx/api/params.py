"""Idiomatic parameter dataclasses for pyrsx high-level API."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal


@dataclass(frozen=True, kw_only=True)
class TriageParams:
    """
    Configuration for `MarkerTable.triage(...)`.

    Using frozen=True + kw_only=True is the modern, explicit, Pythonic way
    for configuration objects (prevents accidental positional args and
    makes the call site very readable).
    """

    min_depth: int = 10
    signif_threshold: float = 0.05
    posterior_threshold: float = 0.9
    bayes_factor_threshold: float = 10.0
    prior: float = 0.01
    linked_prob: float = 0.9
    correction: Literal["bonferroni", "fdr", "none"] = "bonferroni"
    test: Literal["chisq", "fisher", "gtest"] = "chisq"
    output_fasta: bool = False
    group1: str = "M"
    group2: str = "F"

    # Future: validation, __post_init__ checks, etc.
