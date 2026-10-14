// GPL-3.0-or-later

//! Versioned Bayesian configuration with validation and parameter provenance.

use std::collections::BTreeMap;
use std::error::Error;
use std::fmt::{Display, Formatter};

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

use crate::stats::{
    BayesFactorModel, BetaPrior, DirectionalModel, PosteriorModel, PrevalencePrior,
};

/// Where a Bayesian profile field got its value.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum ParameterSource {
    /// The value was present in the input TOML.
    InputProfile,
    /// The value came from a named built-in profile.
    NamedProfile,
    /// The value was supplied on the command line.
    Cli,
    /// The value was filled in from another field.
    Derived,
}

/// Failure while reading, writing, or checking a Bayesian profile.
#[derive(Debug)]
pub enum ProfileError {
    /// The TOML document could not be decoded.
    Parse(toml::de::Error),
    /// The value could not be encoded as TOML.
    Encode(toml::ser::Error),
    /// A field is outside the domain the calculation accepts.
    InvalidValue {
        /// Profile field that failed validation.
        field: &'static str,
        /// Why the field was rejected.
        reason: String,
    },
}

impl Display for ProfileError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Parse(error) => write!(formatter, "invalid Bayesian profile TOML: {error}"),
            Self::Encode(error) => {
                write!(formatter, "cannot encode Bayesian profile TOML: {error}")
            }
            Self::InvalidValue { field, reason } => {
                write!(
                    formatter,
                    "invalid Bayesian profile field {field}: {reason}"
                )
            }
        }
    }
}

impl Error for ProfileError {}

impl From<toml::de::Error> for ProfileError {
    fn from(error: toml::de::Error) -> Self {
        Self::Parse(error)
    }
}

impl From<toml::ser::Error> for ProfileError {
    fn from(error: toml::ser::Error) -> Self {
        Self::Encode(error)
    }
}

/// Groups and depth cutoff recorded with a Bayesian profile.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservationProfile {
    /// First group name. Empty means the caller has not chosen one yet.
    pub group1: String,
    /// Second group name. Empty means the caller has not chosen one yet.
    pub group2: String,
    /// Depth at or above which an individual counts as carrying a marker.
    pub min_depth: u32,
}

/// Priors for the directional sex-linkage model.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ModelProfile {
    /// Prior probability that a marker is sex-linked.
    pub linkage_prior: f64,
    /// Expected prevalence in the linked group.
    pub linked_prevalence: f64,
    /// Expected prevalence under the shared-prevalence null.
    pub null_prevalence: f64,
    /// Mixture weight of the group1-linked direction.
    pub group1_linked_weight: f64,
    /// Prevalence priors used to compute the posterior.
    #[serde(default)]
    pub posterior: Option<PosteriorProfile>,
    /// Beta priors used to compute the Bayes factor.
    #[serde(default = "BayesFactorProfile::uniform_v1")]
    pub bayes_factor: BayesFactorProfile,
}

/// Beta shape parameters stored in a profile.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BetaPriorProfile {
    /// Beta shape alpha. Must be finite and greater than zero.
    pub alpha: f64,
    /// Beta shape beta. Must be finite and greater than zero.
    pub beta: f64,
}

impl BetaPriorProfile {
    /// Fn.
    pub const fn uniform() -> Self {
        Self {
            alpha: 1.0,
            beta: 1.0,
        }
    }

    const fn to_runtime(self) -> BetaPrior {
        BetaPrior {
            alpha: self.alpha,
            beta: self.beta,
        }
    }
}

/// Beta priors for the two group-specific hypotheses and the shared null.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BayesFactorProfile {
    /// Beta prior for prevalence in group 1 under separate prevalences.
    pub alternative_group1: BetaPriorProfile,
    /// Beta prior for prevalence in group 2 under separate prevalences.
    pub alternative_group2: BetaPriorProfile,
    /// Prevalence prior for the shared null.
    pub null: BetaPriorProfile,
}

impl BayesFactorProfile {
    /// Fn.
    pub const fn uniform_v1() -> Self {
        Self {
            alternative_group1: BetaPriorProfile::uniform(),
            alternative_group2: BetaPriorProfile::uniform(),
            null: BetaPriorProfile::uniform(),
        }
    }

    const fn to_runtime(self) -> BayesFactorModel {
        BayesFactorModel {
            alternative_group1: self.alternative_group1.to_runtime(),
            alternative_group2: self.alternative_group2.to_runtime(),
            null: self.null.to_runtime(),
        }
    }
}

/// Prevalence prior stored in a profile: a fixed probability or a Beta distribution.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "family", rename_all = "kebab-case")]
pub enum PrevalencePriorProfile {
    /// Use `probability` as the prevalence.
    Fixed {
        /// Prevalence strictly between 0 and 1.
        probability: f64,
    },
    /// Use a Beta(`alpha`, `beta`) prevalence prior.
    Beta {
        /// Beta shape alpha. Must be finite and greater than zero.
        alpha: f64,
        /// Beta shape beta. Must be finite and greater than zero.
        beta: f64,
    },
}

impl PrevalencePriorProfile {
    const fn to_runtime(self) -> PrevalencePrior {
        match self {
            Self::Fixed { probability } => PrevalencePrior::Fixed { probability },
            Self::Beta { alpha, beta } => PrevalencePrior::Beta(BetaPrior { alpha, beta }),
        }
    }
}

/// Prevalence priors used for the posterior probability of sex linkage.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PosteriorProfile {
    /// Prevalence prior for the sex-linked hypothesis.
    pub linked: PrevalencePriorProfile,
    /// Prevalence prior for the shared null.
    pub null: PrevalencePriorProfile,
}

impl PosteriorProfile {
    /// Fn.
    pub const fn fixed(linked_probability: f64, null_probability: f64) -> Self {
        Self {
            linked: PrevalencePriorProfile::Fixed {
                probability: linked_probability,
            },
            null: PrevalencePriorProfile::Fixed {
                probability: null_probability,
            },
        }
    }

    const fn to_runtime(self) -> PosteriorModel {
        PosteriorModel {
            linked: self.linked.to_runtime(),
            null: self.null.to_runtime(),
        }
    }
}

impl ModelProfile {
    /// Build the runtime model, rejecting parameters outside their domains.
    pub fn to_runtime(&self) -> Result<DirectionalModel, ProfileError> {
        probability("model.linkage_prior", self.linkage_prior)?;
        probability("model.linked_prevalence", self.linked_prevalence)?;
        probability("model.null_prevalence", self.null_prevalence)?;
        probability("model.group1_linked_weight", self.group1_linked_weight)?;
        let posterior = self
            .posterior
            .unwrap_or_else(|| {
                PosteriorProfile::fixed(self.linked_prevalence, self.null_prevalence)
            })
            .to_runtime();
        posterior
            .validate()
            .map_err(|error| ProfileError::InvalidValue {
                field: "model.posterior",
                reason: error.to_string(),
            })?;
        let bayes_factor = self.bayes_factor.to_runtime();
        bayes_factor
            .validate()
            .map_err(|error| ProfileError::InvalidValue {
                field: "model.bayes_factor",
                reason: error.to_string(),
            })?;
        Ok(DirectionalModel {
            linkage_prior: self.linkage_prior,
            linked_prevalence: self.linked_prevalence,
            null_prevalence: self.null_prevalence,
            group1_linked_weight: self.group1_linked_weight,
            posterior,
            bayes_factor,
        })
    }
}

/// Settings for the empirical-Bayes prevalence estimate.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EstimationProfile {
    /// Prior probability that a marker is sex-linked.
    pub linkage_prior: String,
    /// Maximum EM iterations.
    pub max_iterations: u32,
    /// EM convergence tolerance.
    pub tolerance: f64,
    /// Lower bound applied to an estimated probability.
    pub minimum_probability: f64,
    /// Upper bound applied to an estimated probability.
    pub maximum_probability: f64,
}

/// Thresholds that turn evidence into a call.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DecisionProfile {
    /// Posterior probability at or above which a posterior call is true.
    pub posterior_threshold: f64,
    /// Bayes factor at or above which a Bayes-factor call is true.
    pub bayes_factor_threshold: f64,
    /// P-value cutoff before correction.
    pub significance_threshold: f64,
    /// Multiple-testing correction: `bonferroni`, `fdr`, or `none`.
    pub correction: String,
}

/// Test and backend named in a Bayesian profile.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExecutionProfile {
    /// Association test name.
    pub test: String,
    /// Backend name: `cpu` or `cuda`.
    pub pvalue_backend: String,
}

/// Fully specified Bayesian analysis profile.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BayesProfile {
    /// Profile schema. This build accepts `1`.
    pub schema_version: u32,
    /// Name recorded in the profile. It must not be empty.
    pub profile_name: String,
    /// Identifier of the Bayesian model.
    pub model_id: String,
    /// Version of the Bayesian model.
    pub model_version: u32,
    /// Groups and minimum depth for the analysis.
    pub observation: ObservationProfile,
    /// Prior parameters.
    pub model: ModelProfile,
    /// Empirical-Bayes estimation settings.
    pub estimation: EstimationProfile,
    /// Call thresholds.
    pub decision: DecisionProfile,
    /// Test and backend.
    pub execution: ExecutionProfile,
}

/// Bayesian profile as read from TOML, before overrides are applied.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BayesProfileInput {
    /// Profile schema. This build accepts `1`.
    pub schema_version: u32,
    /// Name recorded in the profile. It must not be empty.
    pub profile_name: String,
    /// Identifier of the Bayesian model.
    pub model_id: String,
    /// Version of the Bayesian model.
    pub model_version: u32,
    /// Groups and minimum depth for the analysis.
    pub observation: ObservationProfile,
    /// Prior parameters.
    pub model: ModelProfile,
    /// Empirical-Bayes estimation settings.
    pub estimation: EstimationProfile,
    /// Call thresholds.
    pub decision: DecisionProfile,
    /// Test and backend.
    pub execution: ExecutionProfile,
}

/// Optional replacements applied on top of a Bayesian profile.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ProfileOverrides {
    /// Prior probability that a marker is sex-linked.
    pub linkage_prior: Option<f64>,
    /// Expected prevalence in the linked group.
    pub linked_prevalence: Option<f64>,
    /// Expected prevalence under the shared-prevalence null.
    pub null_prevalence: Option<f64>,
    /// Mixture weight of the group1-linked direction.
    pub group1_linked_weight: Option<f64>,
    /// Posterior probability at or above which a posterior call is true.
    pub posterior_threshold: Option<f64>,
    /// Bayes factor at or above which a Bayes-factor call is true.
    pub bayes_factor_threshold: Option<f64>,
    /// P-value cutoff before correction.
    pub significance_threshold: Option<f64>,
    /// Association test name.
    pub test: Option<String>,
    /// Backend name: `cpu` or `cuda`.
    pub pvalue_backend: Option<String>,
}

/// Bayesian profile after overrides, with a record of each field's source.
#[derive(Clone, Debug, PartialEq)]
pub struct HydratedBayesProfile {
    /// Profile after overrides.
    pub profile: BayesProfile,
    sources: BTreeMap<String, ParameterSource>,
}

impl BayesProfileInput {
    /// Parse and validate a TOML document.
    pub fn parse_toml(input: &str) -> Result<Self, ProfileError> {
        Ok(toml::from_str(input)?)
    }

    /// Apply `overrides` and record the source of each parameter.
    pub fn hydrate(
        &self,
        overrides: &ProfileOverrides,
    ) -> Result<HydratedBayesProfile, ProfileError> {
        let mut profile = BayesProfile {
            schema_version: self.schema_version,
            profile_name: self.profile_name.clone(),
            model_id: self.model_id.clone(),
            model_version: self.model_version,
            observation: self.observation.clone(),
            model: self.model.clone(),
            estimation: self.estimation.clone(),
            decision: self.decision.clone(),
            execution: self.execution.clone(),
        };
        let mut sources = input_sources();

        apply_override(
            &mut profile.model.linkage_prior,
            overrides.linkage_prior,
            "model.linkage_prior",
            &mut sources,
        );
        apply_override(
            &mut profile.model.linked_prevalence,
            overrides.linked_prevalence,
            "model.linked_prevalence",
            &mut sources,
        );
        apply_override(
            &mut profile.model.null_prevalence,
            overrides.null_prevalence,
            "model.null_prevalence",
            &mut sources,
        );
        apply_override(
            &mut profile.model.group1_linked_weight,
            overrides.group1_linked_weight,
            "model.group1_linked_weight",
            &mut sources,
        );
        apply_override(
            &mut profile.decision.posterior_threshold,
            overrides.posterior_threshold,
            "decision.posterior_threshold",
            &mut sources,
        );
        apply_override(
            &mut profile.decision.bayes_factor_threshold,
            overrides.bayes_factor_threshold,
            "decision.bayes_factor_threshold",
            &mut sources,
        );
        apply_override(
            &mut profile.decision.significance_threshold,
            overrides.significance_threshold,
            "decision.significance_threshold",
            &mut sources,
        );
        apply_string_override(
            &mut profile.execution.test,
            overrides.test.as_deref(),
            "execution.test",
            &mut sources,
        );
        apply_string_override(
            &mut profile.execution.pvalue_backend,
            overrides.pvalue_backend.as_deref(),
            "execution.pvalue_backend",
            &mut sources,
        );

        validate(&profile)?;
        Ok(HydratedBayesProfile { profile, sources })
    }
}

impl HydratedBayesProfile {
    /// Recorded origin of `field`, if this profile stored one.
    pub fn source_for(&self, field: &str) -> Option<ParameterSource> {
        self.sources.get(field).copied()
    }

    /// Encode this value as TOML.
    pub fn to_toml(&self) -> Result<String, ProfileError> {
        Ok(toml::to_string_pretty(&self.profile)?)
    }

    /// SHA-256 hex digest of the hydrated profile.
    pub fn digest_sha256(&self) -> String {
        let encoded = toml::to_string(&self.profile)
            .expect("a validated Bayesian profile contains only serializable fields");
        let digest = Sha256::digest(encoded.as_bytes());
        format!("{digest:x}")
    }
}

fn apply_override(
    target: &mut f64,
    override_value: Option<f64>,
    path: &str,
    sources: &mut BTreeMap<String, ParameterSource>,
) {
    if let Some(value) = override_value {
        *target = value;
        sources.insert(path.to_owned(), ParameterSource::Cli);
    }
}

fn apply_string_override(
    target: &mut String,
    override_value: Option<&str>,
    path: &str,
    sources: &mut BTreeMap<String, ParameterSource>,
) {
    if let Some(value) = override_value {
        *target = value.to_owned();
        sources.insert(path.to_owned(), ParameterSource::Cli);
    }
}

fn input_sources() -> BTreeMap<String, ParameterSource> {
    [
        "schema_version",
        "profile_name",
        "model_id",
        "model_version",
        "observation.group1",
        "observation.group2",
        "observation.min_depth",
        "model.linkage_prior",
        "model.linked_prevalence",
        "model.null_prevalence",
        "model.group1_linked_weight",
        "model.posterior.linked.family",
        "model.posterior.linked.probability",
        "model.posterior.linked.alpha",
        "model.posterior.linked.beta",
        "model.posterior.null.family",
        "model.posterior.null.probability",
        "model.posterior.null.alpha",
        "model.posterior.null.beta",
        "model.bayes_factor.alternative_group1.alpha",
        "model.bayes_factor.alternative_group1.beta",
        "model.bayes_factor.alternative_group2.alpha",
        "model.bayes_factor.alternative_group2.beta",
        "model.bayes_factor.null.alpha",
        "model.bayes_factor.null.beta",
        "estimation.linkage_prior",
        "estimation.max_iterations",
        "estimation.tolerance",
        "estimation.minimum_probability",
        "estimation.maximum_probability",
        "decision.posterior_threshold",
        "decision.bayes_factor_threshold",
        "decision.significance_threshold",
        "decision.correction",
        "execution.test",
        "execution.pvalue_backend",
    ]
    .into_iter()
    .map(|path| (path.to_owned(), ParameterSource::InputProfile))
    .collect()
}

fn validate(profile: &BayesProfile) -> Result<(), ProfileError> {
    if profile.schema_version != 1 {
        return invalid("schema_version", "only schema version 1 is supported");
    }
    if profile.model_version == 0 {
        return invalid("model_version", "must be at least 1");
    }
    if profile.profile_name.trim().is_empty() {
        return invalid("profile_name", "must not be empty");
    }
    if profile.model_id.trim().is_empty() {
        return invalid("model_id", "must not be empty");
    }
    if profile.observation.group1 == profile.observation.group2 {
        return invalid("observation.group2", "must differ from observation.group1");
    }
    probability("model.linkage_prior", profile.model.linkage_prior)?;
    probability("model.linked_prevalence", profile.model.linked_prevalence)?;
    probability("model.null_prevalence", profile.model.null_prevalence)?;
    probability(
        "model.group1_linked_weight",
        profile.model.group1_linked_weight,
    )?;
    probability(
        "estimation.minimum_probability",
        profile.estimation.minimum_probability,
    )?;
    probability(
        "estimation.maximum_probability",
        profile.estimation.maximum_probability,
    )?;
    if profile.estimation.minimum_probability >= profile.estimation.maximum_probability {
        return invalid(
            "estimation.maximum_probability",
            "must exceed estimation.minimum_probability",
        );
    }
    if profile.estimation.max_iterations == 0 {
        return invalid("estimation.max_iterations", "must be at least 1");
    }
    if !profile.estimation.tolerance.is_finite() || profile.estimation.tolerance <= 0.0 {
        return invalid(
            "estimation.tolerance",
            "must be finite and greater than zero",
        );
    }
    probability(
        "decision.posterior_threshold",
        profile.decision.posterior_threshold,
    )?;
    probability(
        "decision.significance_threshold",
        profile.decision.significance_threshold,
    )?;
    if !profile.decision.bayes_factor_threshold.is_finite()
        || profile.decision.bayes_factor_threshold <= 0.0
    {
        return invalid(
            "decision.bayes_factor_threshold",
            "must be finite and greater than zero",
        );
    }
    Ok(())
}

fn probability(field: &'static str, value: f64) -> Result<(), ProfileError> {
    if value.is_finite() && value > 0.0 && value < 1.0 {
        Ok(())
    } else {
        invalid(field, "must be finite and strictly between zero and one")
    }
}

fn invalid<T>(field: &'static str, reason: impl Into<String>) -> Result<T, ProfileError> {
    Err(ProfileError::InvalidValue {
        field,
        reason: reason.into(),
    })
}
