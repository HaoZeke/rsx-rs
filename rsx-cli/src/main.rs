// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! CLI for RADSex: sex-determination analysis from RAD-Sequencing data.

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
    about = "rsx: sex-determination from RAD-seq data",
    version,
    args_override_self = true
)]
struct Cli {
    /// Versioned TOML profile supplying a complete command invocation
    #[arg(long, global = true)]
    profile: Option<String>,
    /// Write the fully resolved invocation before analysis
    #[arg(long, global = true)]
    write_hydrated_profile: Option<String>,
    /// Write a software reproducibility bundle before analysis
    #[arg(long, global = true)]
    reproducibility_archive: Option<String>,
    #[command(subcommand)]
    command: Commands,
}

#[derive(Clone, Copy, Debug, Default, ValueEnum)]
enum PosteriorFamily {
    #[default]
    Fixed,
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
    #[default]
    Auto,
    Memory,
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
    /// Compute a table of marker depths from demultiplexed reads files
    Process {
        /// Path to directory containing demultiplexed sequence files
        #[arg(short = 'i', long = "input-dir")]
        input_dir: String,
        /// Path to the output file
        #[arg(short = 'o', long = "output-file")]
        output_file: String,
        /// Number of threads to use
        #[arg(short = 'T', long = "threads", default_value = "1")]
        threads: u32,
        /// Minimum depth in at least one individual to retain a marker
        #[arg(short = 'd', long = "min-depth", default_value = "1")]
        min_depth: u16,
        /// Optional min-hash k-mer dedup (canonical k-mer LSH; heuristic)
        #[arg(short = 'k', long = "kmer-dedup")]
        kmer_dedup: Option<usize>,
    },

    /// Compute the distribution of markers between group1 and group2
    Distrib {
        /// Path to a marker depths table generated by "process"
        #[arg(short = 't', long = "markers-table")]
        markers_table: String,
        /// Path to a population map file
        #[arg(short = 'p', long = "popmap")]
        popmap: String,
        /// Path to the output file
        #[arg(short = 'o', long = "output-file")]
        output_file: String,
        /// Minimum depth to consider a marker present
        #[arg(short = 'd', long = "min-depth", default_value = "1")]
        min_depth: u16,
        /// Names of groups to compare (comma-separated)
        #[arg(short = 'G', long = "groups", value_delimiter = ',')]
        groups: Option<Vec<String>>,
        /// P-value significance threshold
        #[arg(short = 'S', long = "signif-threshold", default_value = "0.05")]
        signif_threshold: f32,
        /// Disable Bonferroni correction (legacy; prefer --correction none)
        #[arg(short = 'C', long = "disable-correction")]
        disable_correction: bool,
        /// Multiple testing correction: bonferroni (default), fdr, none
        #[arg(long = "correction", default_value = "bonferroni")]
        correction: String,
        /// Statistical test: chisq (default), fisher, gtest
        #[arg(long = "test", default_value = "chisq")]
        test_method: String,
        /// Include Bayes Factor and posterior P(sex-linked) per cell
        #[arg(long = "bayes")]
        output_bayes: bool,
        #[command(flatten)]
        bayes_model: BayesModelArgs,
    },

    /// Extract markers significantly associated with a group
    Signif {
        /// Path to a marker depths table
        #[arg(short = 't', long = "markers-table")]
        markers_table: String,
        /// Path to a population map file
        #[arg(short = 'p', long = "popmap")]
        popmap: String,
        /// Path to the output file
        #[arg(short = 'o', long = "output-file")]
        output_file: String,
        /// Minimum depth to consider a marker present
        #[arg(short = 'd', long = "min-depth", default_value = "1")]
        min_depth: u16,
        /// Names of groups to compare
        #[arg(short = 'G', long = "groups", value_delimiter = ',')]
        groups: Option<Vec<String>>,
        /// P-value significance threshold
        #[arg(short = 'S', long = "signif-threshold", default_value = "0.05")]
        signif_threshold: f32,
        /// Multiple testing correction: bonferroni (default), fdr, none
        #[arg(long = "correction", default_value = "bonferroni")]
        correction: String,
        /// Statistical test: chisq (default), fisher, gtest
        #[arg(long = "test", default_value = "chisq")]
        test_method: String,
        /// P-value execution backend: cpu (default) or cuda
        #[arg(long = "backend", default_value = "cpu")]
        backend: String,
        /// Output in FASTA format instead of table format
        #[arg(short = 'a', long = "output-fasta")]
        output_fasta: bool,
        /// Include Bayes Factor and posterior P(sex-linked) in output
        #[arg(long = "bayes")]
        output_bayes: bool,
        #[command(flatten)]
        bayes_model: BayesModelArgs,
    },

    /// Rank strict and Bayesian marker candidates for biological follow-up
    Triage {
        /// Path to a marker depths table
        #[arg(short = 't', long = "markers-table")]
        markers_table: String,
        /// Path to a population map file
        #[arg(short = 'p', long = "popmap")]
        popmap: String,
        /// Path to the output file
        #[arg(short = 'o', long = "output-file")]
        output_file: String,
        /// Minimum depth to consider a marker present
        #[arg(short = 'd', long = "min-depth", default_value = "1")]
        min_depth: u16,
        /// Names of groups to compare
        #[arg(short = 'G', long = "groups", value_delimiter = ',')]
        groups: Option<Vec<String>>,
        /// P-value significance threshold for strict calls
        #[arg(short = 'S', long = "signif-threshold", default_value = "0.05")]
        signif_threshold: f32,
        /// Posterior P(sex-linked) threshold
        #[arg(long = "posterior-threshold", default_value = "0.9")]
        posterior_threshold: f64,
        /// Bayes factor threshold
        #[arg(long = "bayes-factor-threshold", default_value = "10.0")]
        bayes_factor_threshold: f64,
        #[command(flatten)]
        bayes_model: BayesModelArgs,
    },

    /// Compute marker frequencies in all individuals
    Freq {
        /// Path to a marker depths table
        #[arg(short = 't', long = "markers-table")]
        markers_table: String,
        /// Path to the output file
        #[arg(short = 'o', long = "output-file")]
        output_file: String,
        /// Minimum depth to consider a marker present
        #[arg(short = 'd', long = "min-depth", default_value = "1")]
        min_depth: u16,
    },

    /// Compute number of retained reads for each individual
    Depth {
        /// Path to a marker depths table
        #[arg(short = 't', long = "markers-table")]
        markers_table: String,
        /// Path to a population map file
        #[arg(short = 'p', long = "popmap")]
        popmap: String,
        /// Path to the output file
        #[arg(short = 'o', long = "output-file")]
        output_file: String,
        /// Minimum frequency of a marker to retain it
        #[arg(short = 'f', long = "min-frequency", default_value = "0.75")]
        min_frequency: f32,
        /// Depth-table execution policy
        #[arg(long = "streaming-mode", value_enum, default_value = "auto")]
        streaming_mode: DepthStreamingMode,
        /// File-size boundary used by --streaming-mode auto
        #[arg(long = "streaming-threshold-bytes", default_value = "2000000000")]
        streaming_threshold_bytes: u64,
    },

    /// Align markers to a genome and compute metrics
    #[cfg(feature = "map")]
    Map {
        /// Path to a markers file (depth table or FASTA)
        #[arg(short = 't', long = "markers-file")]
        markers_file: String,
        /// Path to the output file
        #[arg(short = 'o', long = "output-file")]
        output_file: String,
        /// Path to a population map file
        #[arg(short = 'p', long = "popmap")]
        popmap: String,
        /// Path to the genome file in FASTA format
        #[arg(short = 'g', long = "genome-file")]
        genome_file: String,
        /// Minimum depth to consider a marker present
        #[arg(short = 'd', long = "min-depth", default_value = "1")]
        min_depth: u16,
        /// Names of groups to compare
        #[arg(short = 'G', long = "groups", value_delimiter = ',')]
        groups: Option<Vec<String>>,
        /// Minimum mapping quality
        #[arg(short = 'q', long = "min-quality", default_value = "20")]
        min_quality: u32,
        /// Minimum frequency of individuals to retain a marker
        #[arg(short = 'Q', long = "min-frequency", default_value = "0.1")]
        min_frequency: f32,
        /// P-value significance threshold
        #[arg(short = 'S', long = "signif-threshold", default_value = "0.05")]
        signif_threshold: f32,
        /// Disable Bonferroni correction
        #[arg(short = 'C', long = "disable-correction")]
        disable_correction: bool,
    },

    /// Extract a subset of a marker depths table
    Subset {
        /// Path to a marker depths table
        #[arg(short = 't', long = "markers-table")]
        markers_table: String,
        /// Path to a population map file
        #[arg(short = 'p', long = "popmap")]
        popmap: String,
        /// Path to the output file
        #[arg(short = 'o', long = "output-file")]
        output_file: String,
        /// Minimum depth to consider a marker present
        #[arg(short = 'd', long = "min-depth", default_value = "1")]
        min_depth: u16,
        /// Names of groups to compare
        #[arg(short = 'G', long = "groups", value_delimiter = ',')]
        groups: Option<Vec<String>>,
        /// P-value significance threshold
        #[arg(short = 'S', long = "signif-threshold", default_value = "0.05")]
        signif_threshold: f32,
        /// Disable Bonferroni correction
        #[arg(short = 'C', long = "disable-correction")]
        disable_correction: bool,
        /// Output in FASTA format
        #[arg(short = 'a', long = "output-fasta")]
        output_fasta: bool,
        /// Minimum individuals from group1
        #[arg(short = 'm', long = "min-group1", default_value = "0")]
        min_group1: u32,
        /// Minimum individuals from group2
        #[arg(short = 'n', long = "min-group2", default_value = "0")]
        min_group2: u32,
        /// Maximum individuals from group1
        #[arg(short = 'M', long = "max-group1", default_value = "9999")]
        max_group1: u32,
        /// Maximum individuals from group2
        #[arg(short = 'N', long = "max-group2", default_value = "9999")]
        max_group2: u32,
        /// Minimum total individuals
        #[arg(short = 'i', long = "min-individuals", default_value = "0")]
        min_individuals: u32,
        /// Maximum total individuals
        #[arg(short = 'I', long = "max-individuals", default_value = "9999")]
        max_individuals: u32,
    },

    /// Merge multiple marker depth tables by sequence identity
    Merge {
        /// Paths to input marker depth tables (positional, 2 or more)
        #[arg(required = true, num_args = 1..)]
        input_files: Vec<String>,
        /// Path to the output merged table
        #[arg(short = 'o', long = "output-file")]
        output_file: String,
        /// Number of entries to buffer before flushing to disk (default: 2M)
        #[arg(short = 'B', long = "buffer-size", default_value = "2000000")]
        buffer_size: usize,
        /// Output as Parquet instead of TSV (requires --features parquet-io)
        #[arg(long = "output-parquet")]
        output_parquet: bool,
    },

    /// Streaming PCA of the depth matrix (sample-space / Tucker mode-2 factors)
    Pca {
        /// Path to a marker depths table
        #[arg(short = 't', long = "markers-table")]
        markers_table: String,
        /// Output directory for eigenvalues, sample scores (loadings.tsv), summary
        #[arg(short = 'o', long = "output-dir")]
        output_dir: String,
        /// Minimum depth to consider a marker present
        #[arg(short = 'd', long = "min-depth", default_value = "1")]
        min_depth: u16,
        /// Number of principal components to output (default: all)
        #[arg(short = 'r', long = "components")]
        components: Option<usize>,
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
                "--groups requires exactly two non-empty group names separated by a comma; got {}",
                g.len()
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
                log::error!("could not create configuration-failure archive: {archive_error}");
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
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, e))?;
            if disable_correction {
                corr = rsx_core::test_method::CorrectionMethod::None;
            }
            let test = rsx_core::test_method::TestMethod::parse_str(&test_method)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, e))?;
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
                .unwrap_or_else(|e| {
                    log::error!("{e}");
                    std::process::exit(1);
                });
            let test =
                rsx_core::test_method::TestMethod::parse_str(&test_method).unwrap_or_else(|e| {
                    log::error!("{e}");
                    std::process::exit(1);
                });
            let backend = rsx_core::compute_backend::PValueBackend::parse_str(&backend)
                .unwrap_or_else(|e| {
                    log::error!("{e}");
                    std::process::exit(1);
                });
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
        } => {
            let (g1, g2) = extract_groups(groups)?;
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
        } => commands::pca::run(&commands::pca::PcaParams {
            markers_table_path: markers_table,
            output_dir,
            min_depth,
            n_components: components,
        }),
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
            log::error!("could not enrich reproducibility archive: {archive_error}");
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
        Some(path) => RunProfile::parse_toml(&fs::read_to_string(path)?)?.profile_name,
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
        .open(&temporary)?;
    file.write_all(profile.to_toml()?.as_bytes())?;
    file.sync_all()?;
    fs::rename(temporary, destination)?;
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        log::error!("{e}");
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
