// GPL-3.0-or-later

//! Strict, versioned configuration for complete rsx command invocations.

use std::error::Error;
use std::fmt::{Display, Formatter};

use serde::{Deserialize, Serialize};

use crate::bayes_profile::ModelProfile;

/// Failure while reading or writing a run profile.
#[derive(Debug)]
pub enum RunProfileError {
    /// The TOML document could not be decoded.
    Parse(toml::de::Error),
    /// The value could not be encoded as TOML.
    Encode(toml::ser::Error),
    /// The document's `schema_version` is not 1.
    UnsupportedSchema(u32),
    /// The profile name is empty or only whitespace.
    EmptyProfileName,
}

impl Display for RunProfileError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Parse(error) => write!(formatter, "invalid run-profile TOML: {error}"),
            Self::Encode(error) => write!(formatter, "cannot encode run-profile TOML: {error}"),
            Self::UnsupportedSchema(version) => {
                write!(
                    formatter,
                    "unsupported run-profile schema version {version}; this build reads schema_version = 1"
                )
            }
            Self::EmptyProfileName => {
                write!(formatter, "run-profile profile_name must not be empty")
            }
        }
    }
}

impl Error for RunProfileError {}

impl From<toml::de::Error> for RunProfileError {
    fn from(error: toml::de::Error) -> Self {
        Self::Parse(error)
    }
}

impl From<toml::ser::Error> for RunProfileError {
    fn from(error: toml::ser::Error) -> Self {
        Self::Encode(error)
    }
}

/// One complete rsx invocation, including the values the command receives.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RunProfile {
    /// Profile schema. This build accepts `1`.
    pub schema_version: u32,
    /// Name recorded in the profile. It must not be empty.
    pub profile_name: String,
    /// Zip path written before analysis, when requested.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reproducibility_archive: Option<String>,
    /// Path of the resolved TOML, when requested.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub write_hydrated_profile: Option<String>,
    /// Command and the arguments that command receives.
    pub run: CommandProfile,
}

/// The command stored in a run profile, tagged by `command` in TOML.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "command", rename_all = "kebab-case")]
pub enum CommandProfile {
    /// A `process` invocation.
    Process(ProcessProfile),
    /// A `distrib` invocation.
    Distrib(DistribProfile),
    /// A `signif` invocation.
    Signif(SignifProfile),
    /// A `triage` invocation.
    Triage(TriageProfile),
    /// A `freq` invocation.
    Freq(FreqProfile),
    /// A `depth` invocation.
    Depth(DepthProfile),
    /// A `map` invocation.
    Map(MapProfile),
    /// A `subset` invocation.
    Subset(SubsetProfile),
    /// A `merge` invocation.
    Merge(MergeProfile),
    /// A `pca` invocation.
    Pca(PcaProfile),
}

/// Resolved arguments for `process`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ProcessProfile {
    /// Directory of demultiplexed FASTA or FASTQ files.
    pub input_dir: String,
    /// File the command writes.
    pub output_file: String,
    /// Worker threads. Must be at least 1.
    pub threads: u32,
    /// Depth at or above which an individual counts as carrying a marker.
    pub min_depth: u16,
    /// Canonical k-mer length used to collapse near-duplicate sequences, or none.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub kmer_dedup: Option<usize>,
}

/// Resolved arguments for `distrib`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DistribProfile {
    /// Marker depth table produced by `process`.
    pub markers_table: String,
    /// Population map path (`individual<TAB>group` per line).
    pub popmap: String,
    /// File the command writes.
    pub output_file: String,
    /// Depth at or above which an individual counts as carrying a marker.
    pub min_depth: u16,
    /// Two group names to compare, or none to use the only two groups in the popmap.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub groups: Option<Vec<String>>,
    /// P-value cutoff before multiple-testing correction.
    pub signif_threshold: f32,
    /// When true, compare the raw p-value with the threshold.
    pub disable_correction: bool,
    /// Multiple-testing correction: `bonferroni`, `fdr`, or `none`.
    pub correction: String,
    /// Association test: `chisq`, `fisher`, or `gtest`.
    pub test_method: String,
    /// When true, append the Bayes factor and the posterior probability.
    pub output_bayes: bool,
    /// Priors used for the Bayes factor and the posterior.
    pub bayes_model: ModelProfile,
}

/// Resolved arguments for `signif`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SignifProfile {
    /// Marker depth table produced by `process`.
    pub markers_table: String,
    /// Population map path (`individual<TAB>group` per line).
    pub popmap: String,
    /// File the command writes.
    pub output_file: String,
    /// Depth at or above which an individual counts as carrying a marker.
    pub min_depth: u16,
    /// Two group names to compare, or none to use the only two groups in the popmap.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub groups: Option<Vec<String>>,
    /// P-value cutoff before multiple-testing correction.
    pub signif_threshold: f32,
    /// Multiple-testing correction: `bonferroni`, `fdr`, or `none`.
    pub correction: String,
    /// Association test: `chisq`, `fisher`, or `gtest`.
    pub test_method: String,
    /// Execution backend: `cpu` or `cuda`.
    pub backend: String,
    /// When true, write FASTA instead of a table.
    pub output_fasta: bool,
    /// When true, append the Bayes factor and the posterior probability.
    pub output_bayes: bool,
    /// Priors used for the Bayes factor and the posterior.
    pub bayes_model: ModelProfile,
}

/// Resolved arguments for `triage`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TriageProfile {
    /// Marker depth table produced by `process`.
    pub markers_table: String,
    /// Population map path (`individual<TAB>group` per line).
    pub popmap: String,
    /// File the command writes.
    pub output_file: String,
    /// Depth at or above which an individual counts as carrying a marker.
    pub min_depth: u16,
    /// Two group names to compare, or none to use the only two groups in the popmap.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub groups: Option<Vec<String>>,
    /// P-value cutoff before multiple-testing correction.
    pub signif_threshold: f32,
    /// Posterior probability at or above which a posterior call is true.
    pub posterior_threshold: f64,
    /// Bayes factor at or above which a Bayes-factor call is true.
    pub bayes_factor_threshold: f64,
    /// Priors used for the Bayes factor and the posterior.
    pub bayes_model: ModelProfile,
}

/// Resolved arguments for `freq`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FreqProfile {
    /// Marker depth table produced by `process`.
    pub markers_table: String,
    /// File the command writes.
    pub output_file: String,
    /// Depth at or above which an individual counts as carrying a marker.
    pub min_depth: u16,
}

const fn default_streaming_threshold_bytes() -> u64 {
    2_000_000_000
}

/// How `depth` chooses between an in-memory median and an external sort.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum StreamingMode {
    /// Use streaming when the markers file exceeds `streaming_threshold_bytes`.
    #[default]
    Auto,
    /// Keep retained depths in memory.
    Memory,
    /// Sort (individual, depth) pairs on disk.
    Streaming,
}

impl StreamingMode {
    /// Fn.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Memory => "memory",
            Self::Streaming => "streaming",
        }
    }
}

/// Resolved arguments for `depth`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DepthProfile {
    /// Marker depth table produced by `process`.
    pub markers_table: String,
    /// Population map path (`individual<TAB>group` per line).
    pub popmap: String,
    /// File the command writes.
    pub output_file: String,
    /// Fraction of individuals that must carry a marker before it is retained.
    pub min_frequency: f32,
    /// Whether the median is computed in memory or by an external sort.
    #[serde(default)]
    pub streaming_mode: StreamingMode,
    /// File size above which `auto` uses the external sort.
    #[serde(default = "default_streaming_threshold_bytes")]
    pub streaming_threshold_bytes: u64,
}

/// Resolved arguments for `map`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MapProfile {
    /// Marker depth table or FASTA file to align.
    pub markers_file: String,
    /// File the command writes.
    pub output_file: String,
    /// Population map path (`individual<TAB>group` per line).
    pub popmap: String,
    /// Reference genome in FASTA format.
    pub genome_file: String,
    /// Depth at or above which an individual counts as carrying a marker.
    pub min_depth: u16,
    /// Two group names to compare, or none to use the only two groups in the popmap.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub groups: Option<Vec<String>>,
    /// Minimum mapping quality retained.
    pub min_quality: u32,
    /// Fraction of individuals that must carry a marker before it is retained.
    pub min_frequency: f32,
    /// P-value cutoff before multiple-testing correction.
    pub signif_threshold: f32,
    /// When true, compare the raw p-value with the threshold.
    pub disable_correction: bool,
}

/// Resolved arguments for `subset`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SubsetProfile {
    /// Marker depth table produced by `process`.
    pub markers_table: String,
    /// Population map path (`individual<TAB>group` per line).
    pub popmap: String,
    /// File the command writes.
    pub output_file: String,
    /// Depth at or above which an individual counts as carrying a marker.
    pub min_depth: u16,
    /// Two group names to compare, or none to use the only two groups in the popmap.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub groups: Option<Vec<String>>,
    /// P-value cutoff before multiple-testing correction.
    pub signif_threshold: f32,
    /// When true, compare the raw p-value with the threshold.
    pub disable_correction: bool,
    /// When true, write FASTA instead of a table.
    pub output_fasta: bool,
    /// Minimum number of group1 individuals that must carry the marker.
    pub min_group1: u32,
    /// Minimum number of group2 individuals that must carry the marker.
    pub min_group2: u32,
    /// Maximum number of group1 individuals that may carry the marker.
    pub max_group1: u32,
    /// Maximum number of group2 individuals that may carry the marker.
    pub max_group2: u32,
    /// Minimum number of individuals, from either group, that must carry the marker.
    pub min_individuals: u32,
    /// Maximum number of individuals, from either group, that may carry the marker.
    pub max_individuals: u32,
}

/// Resolved arguments for `merge`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MergeProfile {
    /// Marker depth tables, in the order their individuals are concatenated.
    pub input_files: Vec<String>,
    /// File the command writes.
    pub output_file: String,
    /// Rows held in memory before a sorted chunk is written.
    pub buffer_size: usize,
    /// When true, write Parquet instead of a TSV depth table.
    pub output_parquet: bool,
}

/// Resolved arguments for `pca`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PcaProfile {
    /// Marker depth table produced by `process`.
    pub markers_table: String,
    /// Directory created for the PCA tables.
    pub output_dir: String,
    /// Depth at or above which an individual counts as carrying a marker.
    pub min_depth: u16,
    /// How many leading components to write, or all of them when omitted.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub components: Option<usize>,
}

impl RunProfile {
    /// Parse and validate a TOML document.
    pub fn parse_toml(input: &str) -> Result<Self, RunProfileError> {
        let profile: Self = toml::from_str(input)?;
        profile.validate()?;
        Ok(profile)
    }

    /// Encode this value as TOML.
    pub fn to_toml(&self) -> Result<String, RunProfileError> {
        self.validate()?;
        Ok(toml::to_string_pretty(self)?)
    }

    /// Subcommand name (`process`, `signif`, and the other `rsx` commands).
    pub fn command_name(&self) -> &'static str {
        match &self.run {
            CommandProfile::Process(_) => "process",
            CommandProfile::Distrib(_) => "distrib",
            CommandProfile::Signif(_) => "signif",
            CommandProfile::Triage(_) => "triage",
            CommandProfile::Freq(_) => "freq",
            CommandProfile::Depth(_) => "depth",
            CommandProfile::Map(_) => "map",
            CommandProfile::Subset(_) => "subset",
            CommandProfile::Merge(_) => "merge",
            CommandProfile::Pca(_) => "pca",
        }
    }

    fn validate(&self) -> Result<(), RunProfileError> {
        if self.schema_version != 1 {
            return Err(RunProfileError::UnsupportedSchema(self.schema_version));
        }
        if self.profile_name.trim().is_empty() {
            return Err(RunProfileError::EmptyProfileName);
        }
        Ok(())
    }
}
