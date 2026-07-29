"""Strict, versioned configuration models for complete rsx invocations."""

from __future__ import annotations

from typing import Annotated, List, Literal, Optional, Union

from pydantic import BaseModel, ConfigDict, Field
from typing_extensions import TypeAliasType

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.9 and 3.10
    import tomli as tomllib


class StrictProfileModel(BaseModel):
    """Base for profile sections that reject coercion and unknown fields."""

    model_config = ConfigDict(extra="forbid", strict=True)


class ProcessProfile(StrictProfileModel):
    command: Literal["process"]
    input_dir: str
    output_file: str
    threads: int = Field(ge=1)
    min_depth: int = Field(ge=0, le=65535)
    kmer_dedup: Optional[int] = Field(default=None, ge=1)


class BetaPriorProfile(StrictProfileModel):
    """Positive shape parameters for a Beta distribution."""

    alpha: float = Field(default=1.0, gt=0.0)
    beta: float = Field(default=1.0, gt=0.0)


class BayesFactorProfile(StrictProfileModel):
    """Priors for separate group prevalences and the shared null prevalence."""

    alternative_group1: BetaPriorProfile = Field(default_factory=BetaPriorProfile)
    alternative_group2: BetaPriorProfile = Field(default_factory=BetaPriorProfile)
    null: BetaPriorProfile = Field(default_factory=BetaPriorProfile)


class DirectionalModelProfile(StrictProfileModel):
    linkage_prior: float = Field(gt=0.0, lt=1.0)
    linked_prevalence: float = Field(gt=0.0, lt=1.0)
    null_prevalence: float = Field(gt=0.0, lt=1.0)
    group1_linked_weight: float = Field(gt=0.0, lt=1.0)
    bayes_factor: BayesFactorProfile = Field(default_factory=BayesFactorProfile)


class DistribProfile(StrictProfileModel):
    command: Literal["distrib"]
    markers_table: str
    popmap: str
    output_file: str
    min_depth: int = Field(ge=0, le=65535)
    groups: Optional[List[str]] = None
    signif_threshold: float = Field(gt=0.0, lt=1.0)
    disable_correction: bool
    correction: str
    test_method: str
    output_bayes: bool
    bayes_model: DirectionalModelProfile


class SignifProfile(StrictProfileModel):
    command: Literal["signif"]
    markers_table: str
    popmap: str
    output_file: str
    min_depth: int = Field(ge=0, le=65535)
    groups: Optional[List[str]] = None
    signif_threshold: float = Field(gt=0.0, lt=1.0)
    correction: str
    test_method: str
    backend: str
    output_fasta: bool
    output_bayes: bool
    bayes_model: DirectionalModelProfile


class TriageProfile(StrictProfileModel):
    command: Literal["triage"]
    markers_table: str
    popmap: str
    output_file: str
    min_depth: int = Field(ge=0, le=65535)
    groups: Optional[List[str]] = None
    signif_threshold: float = Field(gt=0.0, lt=1.0)
    posterior_threshold: float = Field(gt=0.0, lt=1.0)
    bayes_factor_threshold: float = Field(gt=0.0)
    bayes_model: DirectionalModelProfile


class FreqProfile(StrictProfileModel):
    command: Literal["freq"]
    markers_table: str
    output_file: str
    min_depth: int = Field(ge=0, le=65535)


class DepthProfile(StrictProfileModel):
    command: Literal["depth"]
    markers_table: str
    popmap: str
    output_file: str
    min_frequency: float = Field(ge=0.0, le=1.0)


class MapProfile(StrictProfileModel):
    command: Literal["map"]
    markers_file: str
    output_file: str
    popmap: str
    genome_file: str
    min_depth: int = Field(ge=0, le=65535)
    groups: Optional[List[str]] = None
    min_quality: int = Field(ge=0)
    min_frequency: float = Field(ge=0.0, le=1.0)
    signif_threshold: float = Field(gt=0.0, lt=1.0)
    disable_correction: bool


class SubsetProfile(StrictProfileModel):
    command: Literal["subset"]
    markers_table: str
    popmap: str
    output_file: str
    min_depth: int = Field(ge=0, le=65535)
    groups: Optional[List[str]] = None
    signif_threshold: float = Field(gt=0.0, lt=1.0)
    disable_correction: bool
    output_fasta: bool
    min_group1: int = Field(ge=0)
    min_group2: int = Field(ge=0)
    max_group1: int = Field(ge=0)
    max_group2: int = Field(ge=0)
    min_individuals: int = Field(ge=0)
    max_individuals: int = Field(ge=0)


class MergeProfile(StrictProfileModel):
    command: Literal["merge"]
    input_files: List[str] = Field(min_length=2)
    output_file: str
    buffer_size: int = Field(ge=1)
    output_parquet: bool


class PcaProfile(StrictProfileModel):
    command: Literal["pca"]
    markers_table: str
    output_dir: str
    min_depth: int = Field(ge=0, le=65535)
    components: Optional[int] = Field(default=None, ge=1)


RunCommand = TypeAliasType(
    "RunCommand",
    Annotated[
        Union[
            ProcessProfile,
            DistribProfile,
            SignifProfile,
            TriageProfile,
            FreqProfile,
            DepthProfile,
            MapProfile,
            SubsetProfile,
            MergeProfile,
            PcaProfile,
        ],
        Field(discriminator="command"),
    ],
)


class RunProfile(StrictProfileModel):
    """Schema-versioned, fully specified rsx invocation."""

    schema_version: Literal[1]
    profile_name: str = Field(min_length=1)
    reproducibility_archive: Optional[str] = None
    write_hydrated_profile: Optional[str] = None
    run: RunCommand


def parse_run_profile_toml(text: str) -> RunProfile:
    """Parse TOML text and validate it against the versioned run schema."""

    return RunProfile.model_validate(tomllib.loads(text))


__all__ = [
    "BayesFactorProfile",
    "BetaPriorProfile",
    "DirectionalModelProfile",
    "RunCommand",
    "RunProfile",
    "parse_run_profile_toml",
]
