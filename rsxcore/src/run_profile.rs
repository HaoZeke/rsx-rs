// GPL-3.0-or-later

//! Strict, versioned configuration for complete rsx command invocations.

use std::error::Error;
use std::fmt::{Display, Formatter};

use serde::{Deserialize, Serialize};

#[derive(Debug)]
pub enum RunProfileError {
    Parse(toml::de::Error),
    Encode(toml::ser::Error),
    UnsupportedSchema(u32),
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
                    "unsupported run-profile schema version {version}"
                )
            }
            Self::EmptyProfileName => write!(formatter, "run-profile name must not be empty"),
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

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RunProfile {
    pub schema_version: u32,
    pub profile_name: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reproducibility_archive: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub write_hydrated_profile: Option<String>,
    pub run: CommandProfile,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "command", rename_all = "kebab-case")]
pub enum CommandProfile {
    Process(ProcessProfile),
    Distrib(DistribProfile),
    Signif(SignifProfile),
    Triage(TriageProfile),
    Freq(FreqProfile),
    Depth(DepthProfile),
    Map(MapProfile),
    Subset(SubsetProfile),
    Merge(MergeProfile),
    Pca(PcaProfile),
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ProcessProfile {
    pub input_dir: String,
    pub output_file: String,
    pub threads: u32,
    pub min_depth: u16,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub kmer_dedup: Option<usize>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DistribProfile {
    pub markers_table: String,
    pub popmap: String,
    pub output_file: String,
    pub min_depth: u16,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub groups: Option<Vec<String>>,
    pub signif_threshold: f32,
    pub disable_correction: bool,
    pub correction: String,
    pub test_method: String,
    pub output_bayes: bool,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SignifProfile {
    pub markers_table: String,
    pub popmap: String,
    pub output_file: String,
    pub min_depth: u16,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub groups: Option<Vec<String>>,
    pub signif_threshold: f32,
    pub correction: String,
    pub test_method: String,
    pub backend: String,
    pub output_fasta: bool,
    pub output_bayes: bool,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TriageProfile {
    pub markers_table: String,
    pub popmap: String,
    pub output_file: String,
    pub min_depth: u16,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub groups: Option<Vec<String>>,
    pub signif_threshold: f32,
    pub posterior_threshold: f64,
    pub bayes_factor_threshold: f64,
    pub prior_probability: f64,
    pub linked_probability: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FreqProfile {
    pub markers_table: String,
    pub output_file: String,
    pub min_depth: u16,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DepthProfile {
    pub markers_table: String,
    pub popmap: String,
    pub output_file: String,
    pub min_frequency: f32,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MapProfile {
    pub markers_file: String,
    pub output_file: String,
    pub popmap: String,
    pub genome_file: String,
    pub min_depth: u16,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub groups: Option<Vec<String>>,
    pub min_quality: u32,
    pub min_frequency: f32,
    pub signif_threshold: f32,
    pub disable_correction: bool,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SubsetProfile {
    pub markers_table: String,
    pub popmap: String,
    pub output_file: String,
    pub min_depth: u16,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub groups: Option<Vec<String>>,
    pub signif_threshold: f32,
    pub disable_correction: bool,
    pub output_fasta: bool,
    pub min_group1: u32,
    pub min_group2: u32,
    pub max_group1: u32,
    pub max_group2: u32,
    pub min_individuals: u32,
    pub max_individuals: u32,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MergeProfile {
    pub input_files: Vec<String>,
    pub output_file: String,
    pub buffer_size: usize,
    pub output_parquet: bool,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PcaProfile {
    pub markers_table: String,
    pub output_dir: String,
    pub min_depth: u16,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub components: Option<usize>,
}

impl RunProfile {
    pub fn parse_toml(input: &str) -> Result<Self, RunProfileError> {
        let profile: Self = toml::from_str(input)?;
        profile.validate()?;
        Ok(profile)
    }

    pub fn to_toml(&self) -> Result<String, RunProfileError> {
        self.validate()?;
        Ok(toml::to_string_pretty(self)?)
    }

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
