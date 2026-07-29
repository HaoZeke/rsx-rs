"""Idiomatic parameter dataclasses for pyrsx high-level API."""

from __future__ import annotations

from dataclasses import dataclass
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
    null_prevalence: float = 0.5
    group1_linked_weight: float = 0.5
    bf_group1_alpha: float = 1.0
    bf_group1_beta: float = 1.0
    bf_group2_alpha: float = 1.0
    bf_group2_beta: float = 1.0
    bf_null_alpha: float = 1.0
    bf_null_beta: float = 1.0
    posterior_linked_family: Literal["fixed", "beta"] = "fixed"
    posterior_linked_alpha: float = 1.0
    posterior_linked_beta: float = 1.0
    posterior_null_family: Literal["fixed", "beta"] = "fixed"
    posterior_null_alpha: float = 1.0
    posterior_null_beta: float = 1.0
    correction: Literal["bonferroni", "fdr", "none"] = "bonferroni"
    test: Literal["chisq", "fisher", "gtest"] = "chisq"
    output_fasta: bool = False
    group1: str = "M"
    group2: str = "F"
