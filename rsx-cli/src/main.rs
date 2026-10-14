// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Command-line interface for rsx, a RAD-seq sex-marker analysis tool.
//!
//! `rsx --help` lists the commands. `rsx <command> --help` lists that
//! command's arguments and the table it writes. Failures are printed to
//! stderr as `error: ...` and exit with status 1.

use std::ffi::OsString;
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

use clap::{Args, Parser, Subcommand, ValueEnum};
use rsx_core::commands;
use rsx_core::run_profile::{self, CommandProfile, RunProfile};

mod profile_args;
mod repro_archive;

#[derive(Parser)]
#[command(
    name = "rsx",
    version,
    about = "RAD-seq sex-marker analysis",
    long_about = "\
rsx builds a marker depth table from demultiplexed RAD-seq reads and tests \
markers for association with two groups named in a population map. The \
commands follow RADSex (process, distrib, signif, freq, depth, map, subset) \
and add merge, pca, and triage.\n\
\n\
--profile loads a schema_version = 1 TOML file. Flags typed on the command \
line override the profile. --write-hydrated-profile writes the values the \
calculation receives. --reproducibility-archive writes a zip of that \
configuration before the analysis starts; a later failure updates the zip \
instead of deleting it.",
    after_help = "\
Examples:
  rsx process -i reads/ -o markers.tsv -T 8 -d 5
  rsx signif -t markers.tsv -p popmap.tsv -o signif.tsv -G M,F

Run `rsx <command> --help` for arguments and output columns.
-V and --version print the version. -v is also accepted as a version flag.",
    args_override_self = true
)]
struct Cli {
    /// TOML run profile (schema_version = 1) that supplies one complete command
    #[arg(long, global = true, value_name = "FILE")]
    profile: Option<String>,
    /// Write the resolved command, including defaults, before analysis starts
    #[arg(long, global = true, value_name = "FILE")]
    write_hydrated_profile: Option<String>,
    /// Write a zip of the resolved configuration, build, and citation files
    #[arg(long, global = true, value_name = "ZIP")]
    reproducibility_archive: Option<String>,
    #[command(subcommand)]
    command: Commands,
}

#[derive(Clone, Copy, Debug, Default, ValueEnum)]
enum PosteriorFamily {
    /// Use one prevalence for every marker.
    #[default]
    Fixed,
    /// Integrate the prevalence against a Beta(alpha, beta) prior.
    Beta,
}

impl PosteriorFamily {
    fn to_profile(
        self,
        probability: f64,
        alpha: f64,
        beta: f64,
    ) -> rsx_core::bayes_profile::PrevalencePriorProfile {
        match self {
            Self::Fixed => rsx_core::bayes_profile::PrevalencePriorProfile::Fixed { probability },
            Self::Beta => rsx_core::bayes_profile::PrevalencePriorProfile::Beta { alpha, beta },
        }
    }
}

#[derive(Args, Clone, Debug)]
struct BayesModelArgs {
    /// Prior probability that a marker is sex-linked
    #[arg(long = "prior-probability", default_value = "0.01")]
    linkage_prior: f64,
    /// Expected marker prevalence in the linked group
    #[arg(long = "linked-probability", default_value = "0.9")]
    linked_prevalence: f64,
    /// Expected marker prevalence under the directional null model
    #[arg(long = "null-prevalence", default_value = "0.5")]
    null_prevalence: f64,
    /// Mixture weight assigned to the group-1-linked direction
    #[arg(long = "group1-linked-weight", default_value = "0.5")]
    group1_linked_weight: f64,
    /// Posterior prevalence family for the linked group
    #[arg(long = "posterior-linked-family", value_enum, default_value = "fixed")]
    posterior_linked_family: PosteriorFamily,
    /// Fixed linked-group prevalence used by the posterior model
    #[arg(long = "posterior-linked-probability")]
    posterior_linked_probability: Option<f64>,
    /// Alpha shape for a Beta linked-group posterior prevalence
    #[arg(long = "posterior-linked-alpha", default_value = "1.0")]
    posterior_linked_alpha: f64,
    /// Beta shape for a Beta linked-group posterior prevalence
    #[arg(long = "posterior-linked-beta", default_value = "1.0")]
    posterior_linked_beta: f64,
    /// Posterior prevalence family for the shared null
    #[arg(long = "posterior-null-family", value_enum, default_value = "fixed")]
    posterior_null_family: PosteriorFamily,
    /// Fixed prevalence used by the posterior null model
    #[arg(long = "posterior-null-probability")]
    posterior_null_probability: Option<f64>,
    /// Alpha shape for a Beta posterior null prevalence
    #[arg(long = "posterior-null-alpha", default_value = "1.0")]
    posterior_null_alpha: f64,
    /// Beta shape for a Beta posterior null prevalence
    #[arg(long = "posterior-null-beta", default_value = "1.0")]
    posterior_null_beta: f64,
    /// Alpha shape for group 1 under the separate-prevalence hypothesis
    #[arg(long = "bf-group1-alpha", default_value = "1.0")]
    bf_group1_alpha: f64,
    /// Beta shape for group 1 under the separate-prevalence hypothesis
    #[arg(long = "bf-group1-beta", default_value = "1.0")]
    bf_group1_beta: f64,
    /// Alpha shape for group 2 under the separate-prevalence hypothesis
    #[arg(long = "bf-group2-alpha", default_value = "1.0")]
    bf_group2_alpha: f64,
    /// Beta shape for group 2 under the separate-prevalence hypothesis
    #[arg(long = "bf-group2-beta", default_value = "1.0")]
    bf_group2_beta: f64,
    /// Alpha shape for the shared-prevalence null hypothesis
    #[arg(long = "bf-null-alpha", default_value = "1.0")]
    bf_null_alpha: f64,
    /// Beta shape for the shared-prevalence null hypothesis
    #[arg(long = "bf-null-beta", default_value = "1.0")]
    bf_null_beta: f64,
}

impl BayesModelArgs {
    fn to_profile(&self) -> rsx_core::bayes_profile::ModelProfile {
        use rsx_core::bayes_profile::{
            BayesFactorProfile, BetaPriorProfile, ModelProfile, PosteriorProfile,
        };

        ModelProfile {
            linkage_prior: self.linkage_prior,
            linked_prevalence: self.linked_prevalence,
            null_prevalence: self.null_prevalence,
            group1_linked_weight: self.group1_linked_weight,
            posterior: Some(PosteriorProfile {
                linked: self.posterior_linked_family.to_profile(
                    self.posterior_linked_probability
                        .unwrap_or(self.linked_prevalence),
                    self.posterior_linked_alpha,
                    self.posterior_linked_beta,
                ),
                null: self.posterior_null_family.to_profile(
                    self.posterior_null_probability
                        .unwrap_or(self.null_prevalence),
                    self.posterior_null_alpha,
                    self.posterior_null_beta,
                ),
            }),
            bayes_factor: BayesFactorProfile {
                alternative_group1: BetaPriorProfile {
                    alpha: self.bf_group1_alpha,
                    beta: self.bf_group1_beta,
                },
                alternative_group2: BetaPriorProfile {
                    alpha: self.bf_group2_alpha,
                    beta: self.bf_group2_beta,
                },
                null: BetaPriorProfile {
                    alpha: self.bf_null_alpha,
                    beta: self.bf_null_beta,
                },
            },
        }
    }

    fn to_runtime(
        &self,
    ) -> Result<rsx_core::stats::DirectionalModel, rsx_core::bayes_profile::ProfileError> {
        self.to_profile().to_runtime()
    }
}

#[derive(Clone, Copy, Debug, Default, ValueEnum)]
enum DepthStreamingMode {
    /// Stream when the markers file is larger than --streaming-threshold-bytes.
    #[default]
    Auto,
    /// Keep every retained depth in memory and compute the median there.
    Memory,
    /// Compute the median with an external sort of (individual, depth) pairs.
    Streaming,
}

impl DepthStreamingMode {
    fn resolve(self, file_size: u64, threshold_bytes: u64) -> bool {
        match self {
            Self::Auto => file_size > threshold_bytes,
            Self::Memory => false,
            Self::Streaming => true,
        }
    }

    const fn to_profile(self) -> run_profile::StreamingMode {
        match self {
            Self::Auto => run_profile::StreamingMode::Auto,
            Self::Memory => run_profile::StreamingMode::Memory,
            Self::Streaming => run_profile::StreamingMode::Streaming,
        }
    }
}

#[derive(Subcommand)]
enum Commands {
    /// Count marker depths from demultiplexed reads.
    ///
    /// Each file in the directory is one individual. The individual name is
    /// the filename with the sequence extension removed.
    #[command(after_help = "\
Writes a marker depth table:
  #Number of markers : N
  id, sequence, then one column per individual
  one row per retained marker, depths as integers

Accepted read extensions: .fq .fq.gz .fastq .fastq.gz .fa .fa.gz .fasta .fasta.gz .fna .fna.gz")]
    Process {
        /// Directory of demultiplexed FASTA or FASTQ files, optionally gzipped
        #[arg(short = 'i', long = "input-dir", value_name = "DIR")]
        input_dir: String,
        /// Marker depth table to write
        #[arg(short = 'o', long = "output-file", value_name = "FILE")]
        output_file: String,
        /// Worker threads. Must be at least 1
        #[arg(short = 'T', long = "threads", default_value = "1", value_name = "N")]
        threads: u32,
        /// Keep a marker when at least one individual has depth >= this value
        #[arg(short = 'd', long = "min-depth", default_value = "1", value_name = "N")]
        min_depth: u16,
        /// Optional canonical k-mer length used to collapse near-duplicate sequences. Must be >= 1
        #[arg(short = 'k', long = "kmer-dedup", value_name = "K")]
        kmer_dedup: Option<usize>,
    },

    /// Count markers by how many individuals in each group carry them.
    ///
    /// Each output row is one (group1 count, group2 count) cell, not one marker.
    #[command(after_help = "\
Writes a distribution table. The comment line records min_depth, the threshold,
the correction, and the marker count. Columns:
  <group1>  <group2>  Markers  P  CorrectedP  Signif  Bias
With --bayes, Bayes_Factor and Posterior_SexLinked are appended.
--groups takes exactly two names, for example --groups M,F.
If the popmap has exactly two groups, --groups may be omitted.")]
    Distrib {
        /// Marker depth table from `process`
        #[arg(short = 't', long = "markers-table", value_name = "FILE")]
        markers_table: String,
        /// Population map: one `individual<TAB>group` row per line
        #[arg(short = 'p', long = "popmap", value_name = "FILE")]
        popmap: String,
        /// Distribution table to write
        #[arg(short = 'o', long = "output-file", value_name = "FILE")]
        output_file: String,
        /// Depth at or above which an individual counts as carrying the marker
        #[arg(short = 'd', long = "min-depth", default_value = "1", value_name = "N")]
        min_depth: u16,
        /// Exactly two group names, comma-separated, for example M,F
        #[arg(
            short = 'G',
            long = "groups",
            value_delimiter = ',',
            value_name = "A,B"
        )]
        groups: Option<Vec<String>>,
        /// P-value cutoff before correction. Default 0.05. Bonferroni divides it by the marker count
        #[arg(
            short = 'S',
            long = "signif-threshold",
            default_value = "0.05",
            value_name = "P"
        )]
        signif_threshold: f32,
        /// Set the correction to none. Prefer --correction none
        #[arg(short = 'C', long = "disable-correction")]
        disable_correction: bool,
        /// Multiple-testing correction: bonferroni (default), fdr, or none
        #[arg(
            long = "correction",
            default_value = "bonferroni",
            value_name = "METHOD"
        )]
        correction: String,
        /// Association test: chisq (default), fisher, or gtest
        #[arg(long = "test", default_value = "chisq", value_name = "TEST")]
        test_method: String,
        /// Append Bayes_Factor and Posterior_SexLinked to each cell
        #[arg(long = "bayes")]
        output_bayes: bool,
        #[command(flatten)]
        bayes_model: BayesModelArgs,
    },

    /// Write markers whose corrected p-value passes the threshold.
    ///
    /// Table output has the same columns as the input depth table.
    /// `--output-fasta` is rejected together with `--correction fdr`.
    #[command(after_help = "\
Table output starts with a #source:rsx-signif comment, then the input header.
With --bayes, Bayes_Factor and Posterior_SexLinked are appended.
FASTA records look like:
  >id_<group1>:<count>_<group2>:<count>_p:<p>_pcorr:<pcorr>_mindepth:<d>
  SEQUENCE
--groups takes exactly two names, for example --groups M,F.")]
    Signif {
        /// Marker depth table from `process`
        #[arg(short = 't', long = "markers-table", value_name = "FILE")]
        markers_table: String,
        /// Population map: one `individual<TAB>group` row per line
        #[arg(short = 'p', long = "popmap", value_name = "FILE")]
        popmap: String,
        /// Significant-marker table or FASTA file to write
        #[arg(short = 'o', long = "output-file", value_name = "FILE")]
        output_file: String,
        /// Depth at or above which an individual counts as carrying the marker
        #[arg(short = 'd', long = "min-depth", default_value = "1", value_name = "N")]
        min_depth: u16,
        /// Exactly two group names, comma-separated, for example M,F
        #[arg(
            short = 'G',
            long = "groups",
            value_delimiter = ',',
            value_name = "A,B"
        )]
        groups: Option<Vec<String>>,
        /// P-value cutoff before correction. Default 0.05
        #[arg(
            short = 'S',
            long = "signif-threshold",
            default_value = "0.05",
            value_name = "P"
        )]
        signif_threshold: f32,
        /// Multiple-testing correction: bonferroni (default), fdr, or none
        #[arg(
            long = "correction",
            default_value = "bonferroni",
            value_name = "METHOD"
        )]
        correction: String,
        /// Association test: chisq (default), fisher, or gtest
        #[arg(long = "test", default_value = "chisq", value_name = "TEST")]
        test_method: String,
        /// Where the p-values are computed: cpu (default) or cuda
        #[arg(long = "backend", default_value = "cpu", value_name = "BACKEND")]
        backend: String,
        /// Write FASTA instead of a depth table. Not valid with --correction fdr
        #[arg(short = 'a', long = "output-fasta")]
        output_fasta: bool,
        /// Append Bayes_Factor and Posterior_SexLinked to each table row
        #[arg(long = "bayes")]
        output_bayes: bool,
        #[command(flatten)]
        bayes_model: BayesModelArgs,
    },

    /// Rank markers by a strict test and by Bayesian sex-linkage evidence.
    ///
    /// Every marker with at least one individual above `--min-depth` is written.
    #[command(after_help = "\
Writes one row per retained marker. Columns:
  id sequence Group1 Group1_Present Group1_Total Group1_Penetrance
  Group2 Group2_Present Group2_Total Group2_Penetrance
  Bias_Direction Bias P CorrectedP Bayes_Factor Posterior_SexLinked
  Strict_Call Posterior_Call Bayes_Factor_Call Candidate_Class
Candidate_Class is one of strict+posterior, strict_only, posterior_only,
bayes_factor_only, exploratory.
Strict calls use a Bonferroni threshold (signif_threshold / n_markers).")]
    Triage {
        /// Marker depth table from `process`
        #[arg(short = 't', long = "markers-table", value_name = "FILE")]
        markers_table: String,
        /// Population map: one `individual<TAB>group` row per line
        #[arg(short = 'p', long = "popmap", value_name = "FILE")]
        popmap: String,
        /// Ranked-candidate table to write
        #[arg(short = 'o', long = "output-file", value_name = "FILE")]
        output_file: String,
        /// Depth at or above which an individual counts as carrying the marker
        #[arg(short = 'd', long = "min-depth", default_value = "1", value_name = "N")]
        min_depth: u16,
        /// Exactly two group names, comma-separated, for example M,F
        #[arg(
            short = 'G',
            long = "groups",
            value_delimiter = ',',
            value_name = "A,B"
        )]
        groups: Option<Vec<String>>,
        /// Raw p-value cutoff. The strict call uses this divided by the marker count
        #[arg(
            short = 'S',
            long = "signif-threshold",
            default_value = "0.05",
            value_name = "P"
        )]
        signif_threshold: f32,
        /// Posterior P(sex-linked) at or above which Posterior_Call is true
        #[arg(long = "posterior-threshold", default_value = "0.9", value_name = "P")]
        posterior_threshold: f64,
        /// Where the p-value and Bayesian evidence are computed: cpu (default) or cuda
        #[arg(long = "backend", default_value = "cpu", value_name = "BACKEND")]
        backend: String,
        /// Bayes factor at or above which Bayes_Factor_Call is true
        #[arg(
            long = "bayes-factor-threshold",
            default_value = "10.0",
            value_name = "BF"
        )]
        bayes_factor_threshold: f64,
        #[command(flatten)]
        bayes_model: BayesModelArgs,
    },

    /// Histogram of how many individuals carry each marker.
    #[command(after_help = "\
Writes:
  #source:rsx-freq;min_depth:N
  Frequency  Count
Frequency is the number of individuals with depth >= min-depth.
The Frequency 0 row is omitted.")]
    Freq {
        /// Marker depth table from `process`
        #[arg(short = 't', long = "markers-table", value_name = "FILE")]
        markers_table: String,
        /// Frequency table to write
        #[arg(short = 'o', long = "output-file", value_name = "FILE")]
        output_file: String,
        /// Depth at or above which an individual counts as carrying the marker
        #[arg(short = 'd', long = "min-depth", default_value = "1", value_name = "N")]
        min_depth: u16,
    },

    /// Summarize retained reads and marker depths for each individual.
    #[command(after_help = "\
Writes one row per individual. Columns:
  Sample Group Reads Markers Retained Min_depth Max_depth Median_depth Average_depth
Retained counts markers present in at least --min-frequency of individuals.
The median is exact in both memory and streaming mode.
A missing popmap group is left blank; the individual is still reported.")]
    Depth {
        /// Marker depth table from `process`
        #[arg(short = 't', long = "markers-table", value_name = "FILE")]
        markers_table: String,
        /// Population map: one `individual<TAB>group` row per line
        #[arg(short = 'p', long = "popmap", value_name = "FILE")]
        popmap: String,
        /// Per-individual depth summary to write
        #[arg(short = 'o', long = "output-file", value_name = "FILE")]
        output_file: String,
        /// Fraction of individuals that must carry a marker before it enters the median
        #[arg(
            short = 'f',
            long = "min-frequency",
            default_value = "0.75",
            value_name = "F"
        )]
        min_frequency: f32,
        /// auto, memory, or streaming. auto streams above --streaming-threshold-bytes
        #[arg(long = "streaming-mode", value_enum, default_value = "auto")]
        streaming_mode: DepthStreamingMode,
        /// File size in bytes above which --streaming-mode auto uses the external sort
        #[arg(
            long = "streaming-threshold-bytes",
            default_value = "2000000000",
            value_name = "BYTES"
        )]
        streaming_threshold_bytes: u64,
    },

    /// Align markers to a genome and report association at each hit.
    ///
    /// This command is omitted from binaries built without the `map` feature.
    #[command(after_help = "\
Writes one row per alignment that passes --min-quality. Columns:
  Contig Position Length Marker_id Bias P CorrectedP Signif
The markers file is a depth table or FASTA produced by rsx.
The genome file is FASTA. --groups takes exactly two names, for example M,F.")]
    #[cfg(feature = "map")]
    Map {
        /// Marker depth table or FASTA file to align
        #[arg(short = 't', long = "markers-file", value_name = "FILE")]
        markers_file: String,
        /// Alignment table to write
        #[arg(short = 'o', long = "output-file", value_name = "FILE")]
        output_file: String,
        /// Population map: one `individual<TAB>group` row per line
        #[arg(short = 'p', long = "popmap", value_name = "FILE")]
        popmap: String,
        /// Reference genome in FASTA format
        #[arg(short = 'g', long = "genome-file", value_name = "FASTA")]
        genome_file: String,
        /// Depth at or above which an individual counts as carrying the marker
        #[arg(short = 'd', long = "min-depth", default_value = "1", value_name = "N")]
        min_depth: u16,
        /// Exactly two group names, comma-separated, for example M,F
        #[arg(
            short = 'G',
            long = "groups",
            value_delimiter = ',',
            value_name = "A,B"
        )]
        groups: Option<Vec<String>>,
        /// Minimum mapping quality retained
        #[arg(
            short = 'q',
            long = "min-quality",
            default_value = "20",
            value_name = "Q"
        )]
        min_quality: u32,
        /// Minimum fraction of individuals that must carry a marker before it is aligned
        #[arg(
            short = 'Q',
            long = "min-frequency",
            default_value = "0.1",
            value_name = "F"
        )]
        min_frequency: f32,
        /// P-value cutoff before Bonferroni correction. Default 0.05
        #[arg(
            short = 'S',
            long = "signif-threshold",
            default_value = "0.05",
            value_name = "P"
        )]
        signif_threshold: f32,
        /// Compare the raw p-value with the threshold instead of the Bonferroni cutoff
        #[arg(short = 'C', long = "disable-correction")]
        disable_correction: bool,
    },

    /// Keep markers whose group counts fall inside the requested ranges.
    #[command(after_help = "\
Table output is a marker depth table with a #source:rsx-subset comment.
FASTA records use the same header shape as `signif --output-fasta`.
Counts are of individuals with depth >= --min-depth.
--groups takes exactly two names, for example --groups M,F.
Default maxima of 9999 mean those bounds are unset for typical cohorts.")]
    Subset {
        /// Marker depth table from `process`
        #[arg(short = 't', long = "markers-table", value_name = "FILE")]
        markers_table: String,
        /// Population map: one `individual<TAB>group` row per line
        #[arg(short = 'p', long = "popmap", value_name = "FILE")]
        popmap: String,
        /// Filtered depth table or FASTA file to write
        #[arg(short = 'o', long = "output-file", value_name = "FILE")]
        output_file: String,
        /// Depth at or above which an individual counts as carrying the marker
        #[arg(short = 'd', long = "min-depth", default_value = "1", value_name = "N")]
        min_depth: u16,
        /// Exactly two group names, comma-separated, for example M,F
        #[arg(
            short = 'G',
            long = "groups",
            value_delimiter = ',',
            value_name = "A,B"
        )]
        groups: Option<Vec<String>>,
        /// P-value cutoff before Bonferroni correction. Default 0.05
        #[arg(
            short = 'S',
            long = "signif-threshold",
            default_value = "0.05",
            value_name = "P"
        )]
        signif_threshold: f32,
        /// Compare the raw p-value with the threshold instead of the Bonferroni cutoff
        #[arg(short = 'C', long = "disable-correction")]
        disable_correction: bool,
        /// Write FASTA instead of a depth table
        #[arg(short = 'a', long = "output-fasta")]
        output_fasta: bool,
        /// Minimum number of group1 individuals that must carry the marker
        #[arg(
            short = 'm',
            long = "min-group1",
            default_value = "0",
            value_name = "N"
        )]
        min_group1: u32,
        /// Minimum number of group2 individuals that must carry the marker
        #[arg(
            short = 'n',
            long = "min-group2",
            default_value = "0",
            value_name = "N"
        )]
        min_group2: u32,
        /// Maximum number of group1 individuals that may carry the marker
        #[arg(
            short = 'M',
            long = "max-group1",
            default_value = "9999",
            value_name = "N"
        )]
        max_group1: u32,
        /// Maximum number of group2 individuals that may carry the marker
        #[arg(
            short = 'N',
            long = "max-group2",
            default_value = "9999",
            value_name = "N"
        )]
        max_group2: u32,
        /// Minimum number of individuals, from either group, that must carry the marker
        #[arg(
            short = 'i',
            long = "min-individuals",
            default_value = "0",
            value_name = "N"
        )]
        min_individuals: u32,
        /// Maximum number of individuals, from either group, that may carry the marker
        #[arg(
            short = 'I',
            long = "max-individuals",
            default_value = "9999",
            value_name = "N"
        )]
        max_individuals: u32,
    },

    /// Merge marker depth tables by sequence, concatenating individuals.
    #[command(after_help = "\
Input tables are marker depth tables. Individuals are appended in the order
the files are given. Equal sequences have their depths added with saturation
at 65535. The output is a marker depth table unless --output-parquet is set.
Progress messages are written to stderr and are not part of the table.")]
    Merge {
        /// Marker depth tables to merge, in individual-column order
        #[arg(required = true, num_args = 1.., value_name = "TABLE")]
        input_files: Vec<String>,
        /// Merged marker depth table to write
        #[arg(short = 'o', long = "output-file", value_name = "FILE")]
        output_file: String,
        /// Rows held in memory before a sorted chunk is written. Default 2000000
        #[arg(
            short = 'B',
            long = "buffer-size",
            default_value = "2000000",
            value_name = "N"
        )]
        buffer_size: usize,
        /// Write Parquet instead of TSV. Requires a build with the parquet-io feature
        #[arg(long = "output-parquet")]
        output_parquet: bool,
    },

    /// Sample-space PCA of the marker-by-individual depth matrix.
    #[command(after_help = "\
Creates OUTPUT_DIR and writes:
  eigenvalues.tsv     component, eigenvalue, variance_fraction, cumulative
  loadings.tsv        individual, PC1, PC2, ... (same contents as sample_scores.tsv)
  sample_scores.tsv   individual, PC1, PC2, ...
  summary.txt         marker count, individual count, and the leading fractions
--components limits how many components are written. The default is all of them.")]
    Pca {
        /// Marker depth table from `process`
        #[arg(short = 't', long = "markers-table", value_name = "FILE")]
        markers_table: String,
        /// Directory to create for eigenvalues.tsv, loadings.tsv, sample_scores.tsv, and summary.txt
        #[arg(short = 'o', long = "output-dir", value_name = "DIR")]
        output_dir: String,
        /// Depth at or above which an individual counts as carrying the marker
        #[arg(short = 'd', long = "min-depth", default_value = "1", value_name = "N")]
        min_depth: u16,
        /// How many leading components to write. Default: every component
        #[arg(short = 'r', long = "components", value_name = "N")]
        components: Option<usize>,
        /// Where the Gram matrix is accumulated: cpu (default) or cuda
        #[arg(long = "backend", default_value = "cpu", value_name = "BACKEND")]
        backend: String,
    },
}

fn extract_groups(
    groups: &Option<Vec<String>>,
) -> Result<(String, String), Box<dyn std::error::Error>> {
    match groups {
        None => Ok((String::new(), String::new())),
        Some(g) if g.len() == 2 && !g[0].is_empty() && !g[1].is_empty() => {
            Ok((g[0].clone(), g[1].clone()))
        }
        Some(g) => Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            format!(
                "--groups requires exactly two non-empty group names separated by a comma, for example --groups M,F; got {} ({})",
                g.len(),
                if g.is_empty() {
                    "none".to_string()
                } else {
                    g.join(", ")
                }
            ),
        )
        .into()),
    }
}

fn run() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();

    // Handle --version / -V early so it works even with required subcommands
    // (the derive version attr adds the flag; we short-circuit here for robustness).
    let raw_arguments: Vec<OsString> = std::env::args_os().collect();
    if raw_arguments
        .iter()
        .map(|argument| argument.to_string_lossy())
        .any(|argument| argument == "--version" || argument == "-V" || argument == "-v")
    {
        println!("rsx {}", env!("CARGO_PKG_VERSION"));
        std::process::exit(0);
    }

    let cli = match parse_cli_from(raw_arguments.clone()) {
        Ok(cli) => cli,
        Err(error) => {
            if let Err(archive_error) =
                capture_profile_resolution_failure(&raw_arguments, error.as_ref())
            {
                eprintln!("error: could not create configuration-failure archive: {archive_error}");
            }
            return Err(error);
        }
    };
    let hydrated = resolved_run_profile(&cli)?;
    write_hydrated_profile(&hydrated)?;
    let input_profile = match &cli.profile {
        Some(path) => fs::read_to_string(path)?,
        None => hydrated.to_toml()?,
    };
    let archive_path = hydrated.reproducibility_archive.clone();
    if let Some(path) = &archive_path {
        repro_archive::write_archive(
            Path::new(path),
            &input_profile,
            &hydrated,
            repro_archive::ArchiveStatus::Started,
        )?;
    }

    let result: Result<(), Box<dyn std::error::Error>> = match cli.command {
        Commands::Process {
            input_dir,
            output_file,
            threads,
            min_depth,
            kmer_dedup,
        } => {
            let params = commands::process::ProcessParams {
                input_dir_path: input_dir,
                output_file_path: output_file,
                n_threads: threads,
                min_depth,
                kmer_dedup,
            };
            // Use MPI if available (feature-gated), else single-node rayon
            #[cfg(feature = "mpi")]
            {
                commands::process_mpi::run_mpi(&params)
            }
            #[cfg(not(feature = "mpi"))]
            {
                commands::process::run(&params)
            }
        }

        Commands::Distrib {
            markers_table,
            popmap,
            output_file,
            min_depth,
            ref groups,
            signif_threshold,
            disable_correction,
            correction,
            test_method,
            output_bayes,
            bayes_model,
        } => {
            let (g1, g2) = extract_groups(groups)?;
            let mut corr = rsx_core::test_method::CorrectionMethod::parse_str(&correction)
                .map_err(|e| format!("distrib: {e}"))?;
            if disable_correction {
                corr = rsx_core::test_method::CorrectionMethod::None;
            }
            let test = rsx_core::test_method::TestMethod::parse_str(&test_method)
                .map_err(|e| format!("distrib: {e}"))?;
            commands::distrib::run(&commands::distrib::DistribParams {
                markers_table_path: markers_table,
                popmap_file_path: popmap,
                output_file_path: output_file,
                min_depth,
                signif_threshold,
                correction: corr,
                test_method: test,
                output_bayes,
                bayes_model: bayes_model.to_runtime()?,
                group1: g1,
                group2: g2,
            })
        }

        Commands::Signif {
            markers_table,
            popmap,
            output_file,
            min_depth,
            ref groups,
            signif_threshold,
            correction,
            test_method,
            backend,
            output_fasta,
            output_bayes,
            bayes_model,
        } => {
            let (g1, g2) = extract_groups(groups)?;
            let corr = rsx_core::test_method::CorrectionMethod::parse_str(&correction)
                .map_err(|e| format!("signif: {e}"))?;
            let test = rsx_core::test_method::TestMethod::parse_str(&test_method)
                .map_err(|e| format!("signif: {e}"))?;
            let backend = rsx_core::compute_backend::PValueBackend::parse_str(&backend)
                .map_err(|e| format!("signif: {e}"))?;
            commands::signif::run_with_backend(
                &commands::signif::SignifParams {
                    markers_table_path: markers_table,
                    popmap_file_path: popmap,
                    output_file_path: output_file,
                    min_depth,
                    signif_threshold,
                    correction: corr,
                    test_method: test,
                    output_fasta,
                    output_bayes,
                    bayes_model: bayes_model.to_runtime()?,
                    group1: g1,
                    group2: g2,
                },
                backend,
            )
        }

        Commands::Triage {
            markers_table,
            popmap,
            output_file,
            min_depth,
            ref groups,
            signif_threshold,
            posterior_threshold,
            bayes_factor_threshold,
            bayes_model,
            backend,
        } => {
            let (g1, g2) = extract_groups(groups)?;
            let backend = rsx_core::compute_backend::PValueBackend::parse_str(&backend)
                .map_err(|e| format!("triage: {e}"))?;
            commands::triage::run(&commands::triage::TriageParams {
                markers_table_path: markers_table,
                popmap_file_path: popmap,
                output_file_path: output_file,
                min_depth,
                signif_threshold,
                posterior_threshold,
                bayes_factor_threshold,
                bayes_model: bayes_model.to_runtime()?,
                group1: g1,
                group2: g2,
                backend,
            })
        }

        Commands::Freq {
            markers_table,
            output_file,
            min_depth,
        } => commands::freq::run(&commands::freq::FreqParams {
            markers_table_path: markers_table,
            output_file_path: output_file,
            min_depth,
        }),

        Commands::Depth {
            markers_table,
            popmap,
            output_file,
            min_frequency,
            streaming_mode,
            streaming_threshold_bytes,
        } => {
            let file_size = std::fs::metadata(&markers_table)
                .map(|m| m.len())
                .unwrap_or(0);
            let streaming = streaming_mode.resolve(file_size, streaming_threshold_bytes);
            if streaming {
                log::info!(
                    "using streaming depth mode (file bytes: {file_size}, auto threshold bytes: {streaming_threshold_bytes})"
                );
            }
            commands::depth::run(&commands::depth::DepthParams {
                markers_table_path: markers_table,
                popmap_file_path: popmap,
                output_file_path: output_file,
                min_frequency,
                streaming,
            })
        }

        #[cfg(feature = "map")]
        Commands::Map {
            markers_file,
            output_file,
            popmap,
            genome_file,
            min_depth,
            ref groups,
            min_quality,
            min_frequency,
            signif_threshold,
            disable_correction,
        } => {
            let (g1, g2) = extract_groups(groups)?;
            commands::map::run(&commands::map::MapParams {
                markers_table_path: markers_file,
                popmap_file_path: popmap,
                genome_file_path: genome_file,
                output_file_path: output_file,
                min_depth,
                min_quality,
                min_frequency,
                signif_threshold,
                correction: if disable_correction {
                    rsx_core::test_method::CorrectionMethod::None
                } else {
                    rsx_core::test_method::CorrectionMethod::Bonferroni
                },
                test_method: rsx_core::test_method::TestMethod::ChiSquared,
                output_bayes: false,
                group1: g1,
                group2: g2,
            })
        }

        Commands::Subset {
            markers_table,
            popmap,
            output_file,
            min_depth,
            ref groups,
            signif_threshold,
            disable_correction,
            output_fasta,
            min_group1,
            min_group2,
            max_group1,
            max_group2,
            min_individuals,
            max_individuals,
        } => {
            let (g1, g2) = extract_groups(groups)?;
            commands::subset::run(&commands::subset::SubsetParams {
                markers_table_path: markers_table,
                popmap_file_path: popmap,
                output_file_path: output_file,
                min_depth,
                signif_threshold,
                correction: if disable_correction {
                    rsx_core::test_method::CorrectionMethod::None
                } else {
                    rsx_core::test_method::CorrectionMethod::Bonferroni
                },
                test_method: rsx_core::test_method::TestMethod::ChiSquared,
                output_bayes: false,
                output_fasta,
                group1: g1,
                group2: g2,
                min_group1,
                min_group2,
                max_group1,
                max_group2,
                min_individuals,
                max_individuals,
            })
        }

        Commands::Merge {
            input_files,
            output_file,
            buffer_size,
            output_parquet,
        } => commands::merge::run(&commands::merge::MergeParams {
            input_files,
            buffer_size: Some(buffer_size),
            output_file_path: output_file,
            output_parquet,
        }),

        Commands::Pca {
            markers_table,
            output_dir,
            min_depth,
            components,
            backend,
        } => {
            let backend = rsx_core::compute_backend::PValueBackend::parse_str(&backend)
                .map_err(|e| format!("pca: {e}"))?;
            commands::pca::run(&commands::pca::PcaParams {
                markers_table_path: markers_table,
                output_dir,
                min_depth,
                n_components: components,
                backend,
            })
        }
    };

    if let Some(path) = &archive_path {
        let error_message = result.as_ref().err().map(ToString::to_string);
        let status = match error_message.as_deref() {
            Some(message) => repro_archive::ArchiveStatus::Failed(message),
            None => repro_archive::ArchiveStatus::Completed,
        };
        if let Err(archive_error) =
            repro_archive::write_archive(Path::new(path), &input_profile, &hydrated, status)
        {
            if result.is_ok() {
                return Err(archive_error);
            }
            eprintln!("error: could not enrich reproducibility archive: {archive_error}");
        }
    }
    result
}

fn capture_profile_resolution_failure(
    arguments: &[OsString],
    error: &dyn std::error::Error,
) -> Result<(), Box<dyn std::error::Error>> {
    let profile_path = option_value(arguments, "--profile").map(std::path::PathBuf::from);
    let input_profile = match &profile_path {
        Some(path) => fs::read_to_string(path).unwrap_or_default(),
        None => String::new(),
    };
    let archive_path = option_value(arguments, "--reproducibility-archive")
        .or_else(|| archive_path_from_loose_profile(&input_profile));
    if let Some(path) = archive_path {
        repro_archive::write_resolution_failure(
            Path::new(&path),
            &input_profile,
            &error.to_string(),
        )?;
    }
    Ok(())
}

fn option_value(arguments: &[OsString], name: &str) -> Option<String> {
    let prefix = format!("{name}=");
    let mut index = 1;
    while index < arguments.len() {
        let argument = arguments[index].to_string_lossy();
        if argument == name {
            return arguments
                .get(index + 1)
                .map(|value| value.to_string_lossy().into_owned());
        }
        if let Some(value) = argument.strip_prefix(&prefix) {
            if !value.is_empty() {
                return Some(value.to_owned());
            }
        }
        index += 1;
    }
    None
}

fn archive_path_from_loose_profile(input_profile: &str) -> Option<String> {
    toml::from_str::<toml::Value>(input_profile)
        .ok()?
        .get("reproducibility_archive")?
        .as_str()
        .map(str::to_owned)
}

/// Clap reports `--help` and `--version` as errors. They are requests for
/// output rather than failures, so they must not take the failure path, which
/// writes to the error log and would also archive a configuration failure.
fn is_display_request(kind: clap::error::ErrorKind) -> bool {
    matches!(
        kind,
        clap::error::ErrorKind::DisplayHelp
            | clap::error::ErrorKind::DisplayVersion
            | clap::error::ErrorKind::DisplayHelpOnMissingArgumentOrSubcommand
    )
}

fn parse_cli_from<I>(arguments: I) -> Result<Cli, Box<dyn std::error::Error>>
where
    I: IntoIterator<Item = OsString>,
{
    let expanded = profile_args::expand_profile_args(arguments)?;
    match Cli::try_parse_from(expanded) {
        Ok(cli) => Ok(cli),
        // `exit` writes help and version to stdout and exits 0.
        Err(error) if is_display_request(error.kind()) => error.exit(),
        Err(error) => Err(error.into()),
    }
}

fn resolved_run_profile(cli: &Cli) -> Result<RunProfile, Box<dyn std::error::Error>> {
    let profile_name = match &cli.profile {
        Some(path) => {
            let text = fs::read_to_string(path)
                .map_err(|error| format!("cannot read --profile {path}: {error}"))?;
            RunProfile::parse_toml(&text)
                .map_err(|error| {
                    format!("cannot parse --profile {path}: {error}. Expected schema_version = 1")
                })?
                .profile_name
        }
        None => "cli-v1".to_owned(),
    };
    let run = match &cli.command {
        Commands::Process {
            input_dir,
            output_file,
            threads,
            min_depth,
            kmer_dedup,
        } => CommandProfile::Process(run_profile::ProcessProfile {
            input_dir: input_dir.clone(),
            output_file: output_file.clone(),
            threads: *threads,
            min_depth: *min_depth,
            kmer_dedup: *kmer_dedup,
        }),
        Commands::Distrib {
            markers_table,
            popmap,
            output_file,
            min_depth,
            groups,
            signif_threshold,
            disable_correction,
            correction,
            test_method,
            output_bayes,
            bayes_model,
        } => CommandProfile::Distrib(run_profile::DistribProfile {
            markers_table: markers_table.clone(),
            popmap: popmap.clone(),
            output_file: output_file.clone(),
            min_depth: *min_depth,
            groups: groups.clone(),
            signif_threshold: *signif_threshold,
            disable_correction: *disable_correction,
            correction: correction.clone(),
            test_method: test_method.clone(),
            output_bayes: *output_bayes,
            bayes_model: bayes_model.to_profile(),
        }),
        Commands::Signif {
            markers_table,
            popmap,
            output_file,
            min_depth,
            groups,
            signif_threshold,
            correction,
            test_method,
            backend,
            output_fasta,
            output_bayes,
            bayes_model,
        } => CommandProfile::Signif(run_profile::SignifProfile {
            markers_table: markers_table.clone(),
            popmap: popmap.clone(),
            output_file: output_file.clone(),
            min_depth: *min_depth,
            groups: groups.clone(),
            signif_threshold: *signif_threshold,
            correction: correction.clone(),
            test_method: test_method.clone(),
            backend: backend.clone(),
            output_fasta: *output_fasta,
            output_bayes: *output_bayes,
            bayes_model: bayes_model.to_profile(),
        }),
        Commands::Triage {
            markers_table,
            popmap,
            output_file,
            min_depth,
            groups,
            signif_threshold,
            posterior_threshold,
            bayes_factor_threshold,
            bayes_model,
            backend: _,
        } => CommandProfile::Triage(run_profile::TriageProfile {
            markers_table: markers_table.clone(),
            popmap: popmap.clone(),
            output_file: output_file.clone(),
            min_depth: *min_depth,
            groups: groups.clone(),
            signif_threshold: *signif_threshold,
            posterior_threshold: *posterior_threshold,
            bayes_factor_threshold: *bayes_factor_threshold,
            bayes_model: bayes_model.to_profile(),
        }),
        Commands::Freq {
            markers_table,
            output_file,
            min_depth,
        } => CommandProfile::Freq(run_profile::FreqProfile {
            markers_table: markers_table.clone(),
            output_file: output_file.clone(),
            min_depth: *min_depth,
        }),
        Commands::Depth {
            markers_table,
            popmap,
            output_file,
            min_frequency,
            streaming_mode,
            streaming_threshold_bytes,
        } => CommandProfile::Depth(run_profile::DepthProfile {
            markers_table: markers_table.clone(),
            popmap: popmap.clone(),
            output_file: output_file.clone(),
            min_frequency: *min_frequency,
            streaming_mode: streaming_mode.to_profile(),
            streaming_threshold_bytes: *streaming_threshold_bytes,
        }),
        #[cfg(feature = "map")]
        Commands::Map {
            markers_file,
            output_file,
            popmap,
            genome_file,
            min_depth,
            groups,
            min_quality,
            min_frequency,
            signif_threshold,
            disable_correction,
        } => CommandProfile::Map(run_profile::MapProfile {
            markers_file: markers_file.clone(),
            output_file: output_file.clone(),
            popmap: popmap.clone(),
            genome_file: genome_file.clone(),
            min_depth: *min_depth,
            groups: groups.clone(),
            min_quality: *min_quality,
            min_frequency: *min_frequency,
            signif_threshold: *signif_threshold,
            disable_correction: *disable_correction,
        }),
        Commands::Subset {
            markers_table,
            popmap,
            output_file,
            min_depth,
            groups,
            signif_threshold,
            disable_correction,
            output_fasta,
            min_group1,
            min_group2,
            max_group1,
            max_group2,
            min_individuals,
            max_individuals,
        } => CommandProfile::Subset(run_profile::SubsetProfile {
            markers_table: markers_table.clone(),
            popmap: popmap.clone(),
            output_file: output_file.clone(),
            min_depth: *min_depth,
            groups: groups.clone(),
            signif_threshold: *signif_threshold,
            disable_correction: *disable_correction,
            output_fasta: *output_fasta,
            min_group1: *min_group1,
            min_group2: *min_group2,
            max_group1: *max_group1,
            max_group2: *max_group2,
            min_individuals: *min_individuals,
            max_individuals: *max_individuals,
        }),
        Commands::Merge {
            input_files,
            output_file,
            buffer_size,
            output_parquet,
        } => CommandProfile::Merge(run_profile::MergeProfile {
            input_files: input_files.clone(),
            output_file: output_file.clone(),
            buffer_size: *buffer_size,
            output_parquet: *output_parquet,
        }),
        Commands::Pca {
            markers_table,
            output_dir,
            min_depth,
            components,
            backend: _,
        } => CommandProfile::Pca(run_profile::PcaProfile {
            markers_table: markers_table.clone(),
            output_dir: output_dir.clone(),
            min_depth: *min_depth,
            components: *components,
        }),
    };
    Ok(RunProfile {
        schema_version: 1,
        profile_name,
        reproducibility_archive: cli.reproducibility_archive.clone(),
        write_hydrated_profile: cli.write_hydrated_profile.clone(),
        run,
    })
}

fn write_hydrated_profile(profile: &RunProfile) -> Result<(), Box<dyn std::error::Error>> {
    let Some(destination) = &profile.write_hydrated_profile else {
        return Ok(());
    };
    let destination = Path::new(destination);
    let parent = destination.parent().unwrap_or_else(|| Path::new("."));
    let file_name = destination
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("profile.hydrated.toml");
    static TEMP_SEQUENCE: AtomicU64 = AtomicU64::new(0);
    let sequence = TEMP_SEQUENCE.fetch_add(1, Ordering::Relaxed);
    let temporary = parent.join(format!(
        ".{file_name}.rsx-{}-{sequence}.tmp",
        std::process::id()
    ));
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&temporary)
        .map_err(|error| {
            format!(
                "cannot create hydrated profile {} (temporary file {}): {error}",
                destination.display(),
                temporary.display()
            )
        })?;
    file.write_all(profile.to_toml()?.as_bytes())
        .map_err(|error| {
            format!(
                "cannot write hydrated profile {}: {error}",
                destination.display()
            )
        })?;
    file.sync_all().map_err(|error| {
        format!(
            "cannot finish writing hydrated profile {}: {error}",
            destination.display()
        )
    })?;
    fs::rename(&temporary, destination).map_err(|error| {
        format!(
            "cannot replace hydrated profile {} with {}: {error}",
            destination.display(),
            temporary.display()
        )
    })?;
    Ok(())
}

fn report_error(error: &dyn std::error::Error) {
    let message = error.to_string();
    if message.starts_with("error:") {
        if message.ends_with('\n') {
            eprint!("{message}");
        } else {
            eprintln!("{message}");
        }
    } else {
        eprintln!("error: {message}");
    }
}

fn main() {
    if let Err(error) = run() {
        report_error(error.as_ref());
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use std::fs;

    use super::{Commands, extract_groups, is_display_request, parse_cli_from};
    use clap::Parser;

    #[test]
    fn help_and_version_are_display_requests_not_failures() {
        for arguments in [
            vec!["rsx", "--help"],
            vec!["rsx", "-h"],
            vec!["rsx", "freq", "--help"],
            vec!["rsx", "--version"],
        ] {
            let error = match super::Cli::try_parse_from(&arguments) {
                Ok(_) => panic!("{arguments:?} should be reported as a clap error"),
                Err(error) => error,
            };
            assert!(
                is_display_request(error.kind()),
                "{arguments:?} produced {:?}, which would exit non-zero",
                error.kind()
            );
        }
    }

    #[test]
    fn a_genuine_parse_error_is_not_a_display_request() {
        let error = match super::Cli::try_parse_from(["rsx", "no-such-command"]) {
            Ok(_) => panic!("an unknown subcommand must fail"),
            Err(error) => error,
        };
        assert!(!is_display_request(error.kind()));
    }

    #[test]
    fn missing_groups_uses_popmap_resolution() {
        let groups = None;
        let resolved = extract_groups(&groups).expect("missing groups are valid");
        assert_eq!(resolved, (String::new(), String::new()));
    }

    #[test]
    fn pair_groups_are_accepted() {
        let groups = Some(vec!["male".to_string(), "female".to_string()]);
        let resolved = extract_groups(&groups).expect("two groups are valid");
        assert_eq!(resolved, ("male".to_string(), "female".to_string()));
    }

    #[test]
    fn malformed_groups_are_rejected() {
        let groups = Some(vec!["male".to_string()]);
        let err = extract_groups(&groups).expect_err("single group must fail");
        assert!(
            err.to_string()
                .contains("exactly two non-empty group names")
        );
    }

    #[test]
    fn profile_supplies_required_arguments_and_cli_overrides_values() {
        let path = std::env::temp_dir().join(format!(
            "rsx-process-profile-{}-{}.toml",
            std::process::id(),
            std::thread::current().name().unwrap_or("test")
        ));
        fs::write(
            &path,
            r#"
schema_version = 1
profile_name = "process-v1"

[run]
command = "process"
input_dir = "reads"
output_file = "markers.tsv"
threads = 4
min_depth = 2
kmer_dedup = 31
"#,
        )
        .unwrap();

        let cli = parse_cli_from([
            "rsx".into(),
            "--profile".into(),
            path.clone().into_os_string(),
            "--threads".into(),
            "8".into(),
        ])
        .unwrap();

        match cli.command {
            Commands::Process {
                input_dir,
                output_file,
                threads,
                min_depth,
                kmer_dedup,
            } => {
                assert_eq!(input_dir, "reads");
                assert_eq!(output_file, "markers.tsv");
                assert_eq!(threads, 8);
                assert_eq!(min_depth, 2);
                assert_eq!(kmer_dedup, Some(31));
            }
            _ => panic!("profile selected the wrong command"),
        }

        fs::remove_file(path).unwrap();
    }

    #[test]
    fn bayesian_profile_fields_are_cli_overridable() {
        let path =
            std::env::temp_dir().join(format!("rsx-distrib-profile-{}.toml", std::process::id()));
        fs::write(
            &path,
            r#"
schema_version = 1
profile_name = "distrib-bayes-v1"

[run]
command = "distrib"
markers_table = "markers.tsv"
popmap = "popmap.tsv"
output_file = "distrib.tsv"
min_depth = 1
groups = ["M", "F"]
signif_threshold = 0.05
disable_correction = false
correction = "bonferroni"
test_method = "chisq"
output_bayes = true

[run.bayes_model]
linkage_prior = 0.01
linked_prevalence = 0.9
null_prevalence = 0.5
group1_linked_weight = 0.5
"#,
        )
        .unwrap();

        let cli = parse_cli_from([
            "rsx".into(),
            "--profile".into(),
            path.clone().into_os_string(),
            "--null-prevalence".into(),
            "0.4".into(),
        ])
        .unwrap();

        match cli.command {
            Commands::Distrib { bayes_model, .. } => {
                assert_eq!(bayes_model.linkage_prior, 0.01);
                assert_eq!(bayes_model.linked_prevalence, 0.9);
                assert_eq!(bayes_model.null_prevalence, 0.4);
                assert_eq!(bayes_model.group1_linked_weight, 0.5);
                assert_eq!(bayes_model.bf_group1_alpha, 1.0);
                assert_eq!(bayes_model.bf_null_beta, 1.0);
            }
            _ => panic!("profile selected the wrong command"),
        }

        fs::remove_file(path).unwrap();
    }
}
