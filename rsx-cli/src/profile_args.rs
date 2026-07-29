// GPL-3.0-or-later

use std::error::Error;
use std::ffi::{OsStr, OsString};
use std::fmt::Display;
use std::fs;
use std::io;
use std::path::PathBuf;

use rsx_core::run_profile::{CommandProfile, RunProfile};

const COMMANDS: &[&str] = &[
    "process", "distrib", "signif", "triage", "freq", "depth", "map", "subset", "merge", "pca",
];

pub fn expand_profile_args<I>(arguments: I) -> Result<Vec<OsString>, Box<dyn Error>>
where
    I: IntoIterator<Item = OsString>,
{
    let raw: Vec<OsString> = arguments.into_iter().collect();
    let Some(profile_path) = find_profile_path(&raw)? else {
        return Ok(raw);
    };
    let source = fs::read_to_string(&profile_path)?;
    let profile = RunProfile::parse_toml(&source)?;
    let command_name = profile.command_name();

    let executable = raw
        .first()
        .cloned()
        .unwrap_or_else(|| OsString::from("rsx"));
    let mut expanded = vec![executable, OsString::from("--profile")];
    expanded.push(profile_path.into_os_string());
    if let Some(path) = &profile.reproducibility_archive {
        option(&mut expanded, "--reproducibility-archive", path);
    }
    if let Some(path) = &profile.write_hydrated_profile {
        option(&mut expanded, "--write-hydrated-profile", path);
    }
    expanded.extend(command_args(&profile));
    expanded.extend(user_overrides(&raw, command_name)?);
    Ok(expanded)
}

fn find_profile_path(arguments: &[OsString]) -> Result<Option<PathBuf>, Box<dyn Error>> {
    let mut path = None;
    let mut index = 1;
    while index < arguments.len() {
        let argument = arguments[index].to_string_lossy();
        if argument == "--profile" {
            let value = arguments.get(index + 1).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidInput, "--profile requires a path")
            })?;
            path = Some(PathBuf::from(value));
            index += 2;
            continue;
        }
        if let Some(value) = argument.strip_prefix("--profile=") {
            if value.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "--profile requires a path",
                )
                .into());
            }
            path = Some(PathBuf::from(value));
        }
        index += 1;
    }
    Ok(path)
}

fn user_overrides(
    arguments: &[OsString],
    profile_command: &str,
) -> Result<Vec<OsString>, Box<dyn Error>> {
    let mut overrides = Vec::new();
    let mut command_seen = false;
    let mut index = 1;
    while index < arguments.len() {
        let argument = &arguments[index];
        let text = argument.to_string_lossy();
        if text == "--profile" {
            index += 2;
            continue;
        }
        if text.starts_with("--profile=") {
            index += 1;
            continue;
        }
        if COMMANDS.contains(&text.as_ref()) {
            if text != profile_command {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "profile selects command {profile_command}, but the command line selects {text}"
                    ),
                )
                .into());
            }
            if command_seen {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("command {profile_command} was supplied more than once"),
                )
                .into());
            }
            command_seen = true;
            index += 1;
            continue;
        }
        overrides.push(argument.clone());
        index += 1;
    }
    Ok(overrides)
}

fn command_args(profile: &RunProfile) -> Vec<OsString> {
    let mut args = vec![OsString::from(profile.command_name())];
    match &profile.run {
        CommandProfile::Process(value) => {
            option(&mut args, "--input-dir", &value.input_dir);
            option(&mut args, "--output-file", &value.output_file);
            number(&mut args, "--threads", value.threads);
            number(&mut args, "--min-depth", value.min_depth);
            optional_number(&mut args, "--kmer-dedup", value.kmer_dedup);
        }
        CommandProfile::Distrib(value) => {
            option(&mut args, "--markers-table", &value.markers_table);
            option(&mut args, "--popmap", &value.popmap);
            option(&mut args, "--output-file", &value.output_file);
            number(&mut args, "--min-depth", value.min_depth);
            groups(&mut args, &value.groups);
            number(&mut args, "--signif-threshold", value.signif_threshold);
            flag(&mut args, "--disable-correction", value.disable_correction);
            option(&mut args, "--correction", &value.correction);
            option(&mut args, "--test", &value.test_method);
            flag(&mut args, "--bayes", value.output_bayes);
            bayes_model(&mut args, &value.bayes_model);
        }
        CommandProfile::Signif(value) => {
            option(&mut args, "--markers-table", &value.markers_table);
            option(&mut args, "--popmap", &value.popmap);
            option(&mut args, "--output-file", &value.output_file);
            number(&mut args, "--min-depth", value.min_depth);
            groups(&mut args, &value.groups);
            number(&mut args, "--signif-threshold", value.signif_threshold);
            option(&mut args, "--correction", &value.correction);
            option(&mut args, "--test", &value.test_method);
            option(&mut args, "--backend", &value.backend);
            flag(&mut args, "--output-fasta", value.output_fasta);
            flag(&mut args, "--bayes", value.output_bayes);
            bayes_model(&mut args, &value.bayes_model);
        }
        CommandProfile::Triage(value) => {
            option(&mut args, "--markers-table", &value.markers_table);
            option(&mut args, "--popmap", &value.popmap);
            option(&mut args, "--output-file", &value.output_file);
            number(&mut args, "--min-depth", value.min_depth);
            groups(&mut args, &value.groups);
            number(&mut args, "--signif-threshold", value.signif_threshold);
            number(
                &mut args,
                "--posterior-threshold",
                value.posterior_threshold,
            );
            number(
                &mut args,
                "--bayes-factor-threshold",
                value.bayes_factor_threshold,
            );
            bayes_model(&mut args, &value.bayes_model);
        }
        CommandProfile::Freq(value) => {
            option(&mut args, "--markers-table", &value.markers_table);
            option(&mut args, "--output-file", &value.output_file);
            number(&mut args, "--min-depth", value.min_depth);
        }
        CommandProfile::Depth(value) => {
            option(&mut args, "--markers-table", &value.markers_table);
            option(&mut args, "--popmap", &value.popmap);
            option(&mut args, "--output-file", &value.output_file);
            number(&mut args, "--min-frequency", value.min_frequency);
            option(&mut args, "--streaming-mode", value.streaming_mode.as_str());
            number(
                &mut args,
                "--streaming-threshold-bytes",
                value.streaming_threshold_bytes,
            );
        }
        CommandProfile::Map(value) => {
            option(&mut args, "--markers-file", &value.markers_file);
            option(&mut args, "--output-file", &value.output_file);
            option(&mut args, "--popmap", &value.popmap);
            option(&mut args, "--genome-file", &value.genome_file);
            number(&mut args, "--min-depth", value.min_depth);
            groups(&mut args, &value.groups);
            number(&mut args, "--min-quality", value.min_quality);
            number(&mut args, "--min-frequency", value.min_frequency);
            number(&mut args, "--signif-threshold", value.signif_threshold);
            flag(&mut args, "--disable-correction", value.disable_correction);
        }
        CommandProfile::Subset(value) => {
            option(&mut args, "--markers-table", &value.markers_table);
            option(&mut args, "--popmap", &value.popmap);
            option(&mut args, "--output-file", &value.output_file);
            number(&mut args, "--min-depth", value.min_depth);
            groups(&mut args, &value.groups);
            number(&mut args, "--signif-threshold", value.signif_threshold);
            flag(&mut args, "--disable-correction", value.disable_correction);
            flag(&mut args, "--output-fasta", value.output_fasta);
            number(&mut args, "--min-group1", value.min_group1);
            number(&mut args, "--min-group2", value.min_group2);
            number(&mut args, "--max-group1", value.max_group1);
            number(&mut args, "--max-group2", value.max_group2);
            number(&mut args, "--min-individuals", value.min_individuals);
            number(&mut args, "--max-individuals", value.max_individuals);
        }
        CommandProfile::Merge(value) => {
            args.extend(
                value
                    .input_files
                    .iter()
                    .map(|path| OsString::from(path.as_str())),
            );
            option(&mut args, "--output-file", &value.output_file);
            number(&mut args, "--buffer-size", value.buffer_size);
            flag(&mut args, "--output-parquet", value.output_parquet);
        }
        CommandProfile::Pca(value) => {
            option(&mut args, "--markers-table", &value.markers_table);
            option(&mut args, "--output-dir", &value.output_dir);
            number(&mut args, "--min-depth", value.min_depth);
            optional_number(&mut args, "--components", value.components);
        }
    }
    args
}

fn option(args: &mut Vec<OsString>, name: &str, value: impl AsRef<OsStr>) {
    args.push(OsString::from(name));
    args.push(value.as_ref().to_owned());
}

fn number(args: &mut Vec<OsString>, name: &str, value: impl Display) {
    option(args, name, value.to_string());
}

fn optional_number(args: &mut Vec<OsString>, name: &str, value: Option<impl Display>) {
    if let Some(value) = value {
        number(args, name, value);
    }
}

fn groups(args: &mut Vec<OsString>, value: &Option<Vec<String>>) {
    if let Some(value) = value {
        option(args, "--groups", value.join(","));
    }
}

fn flag(args: &mut Vec<OsString>, name: &str, enabled: bool) {
    if enabled {
        args.push(OsString::from(name));
    }
}

fn bayes_model(args: &mut Vec<OsString>, model: &rsx_core::bayes_profile::ModelProfile) {
    number(args, "--prior-probability", model.linkage_prior);
    number(args, "--linked-probability", model.linked_prevalence);
    number(args, "--null-prevalence", model.null_prevalence);
    number(args, "--group1-linked-weight", model.group1_linked_weight);
    if let Some(posterior) = model.posterior {
        posterior_prior(args, "linked", posterior.linked);
        posterior_prior(args, "null", posterior.null);
    }
    number(
        args,
        "--bf-group1-alpha",
        model.bayes_factor.alternative_group1.alpha,
    );
    number(
        args,
        "--bf-group1-beta",
        model.bayes_factor.alternative_group1.beta,
    );
    number(
        args,
        "--bf-group2-alpha",
        model.bayes_factor.alternative_group2.alpha,
    );
    number(
        args,
        "--bf-group2-beta",
        model.bayes_factor.alternative_group2.beta,
    );
    number(args, "--bf-null-alpha", model.bayes_factor.null.alpha);
    number(args, "--bf-null-beta", model.bayes_factor.null.beta);
}

fn posterior_prior(
    args: &mut Vec<OsString>,
    role: &str,
    prior: rsx_core::bayes_profile::PrevalencePriorProfile,
) {
    use rsx_core::bayes_profile::PrevalencePriorProfile;

    match prior {
        PrevalencePriorProfile::Fixed { probability } => {
            option(args, &format!("--posterior-{role}-family"), "fixed");
            number(
                args,
                &format!("--posterior-{role}-probability"),
                probability,
            );
        }
        PrevalencePriorProfile::Beta { alpha, beta } => {
            option(args, &format!("--posterior-{role}-family"), "beta");
            number(args, &format!("--posterior-{role}-alpha"), alpha);
            number(args, &format!("--posterior-{role}-beta"), beta);
        }
    }
}
