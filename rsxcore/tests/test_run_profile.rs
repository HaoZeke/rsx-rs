// GPL-3.0-or-later

use rsx_core::run_profile::{CommandProfile, RunProfile};
use rsx_core::stats::{bayes_factor_2x2, bayes_factor_2x2_with_model};

const PROFILE_HEADER: &str = r#"
schema_version = 1
profile_name = "test-profile"
reproducibility_archive = "run-repro.zip"
write_hydrated_profile = "run.hydrated.toml"
"#;

#[test]
fn every_analysis_command_has_a_strict_toml_variant() {
    let commands = [
        (
            "process",
            r#"
[run]
command = "process"
input_dir = "reads"
output_file = "markers.tsv"
threads = 4
min_depth = 2
kmer_dedup = 31
"#,
        ),
        (
            "distrib",
            r#"
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
output_bayes = false

[run.bayes_model]
linkage_prior = 0.01
linked_prevalence = 0.9
null_prevalence = 0.5
group1_linked_weight = 0.5
"#,
        ),
        (
            "signif",
            r#"
[run]
command = "signif"
markers_table = "markers.tsv"
popmap = "popmap.tsv"
output_file = "signif.tsv"
min_depth = 1
groups = ["M", "F"]
signif_threshold = 0.05
correction = "fdr"
test_method = "fisher"
backend = "cpu"
output_fasta = false
output_bayes = false

[run.bayes_model]
linkage_prior = 0.01
linked_prevalence = 0.9
null_prevalence = 0.5
group1_linked_weight = 0.5
"#,
        ),
        (
            "triage",
            r#"
[run]
command = "triage"
markers_table = "markers.tsv"
popmap = "popmap.tsv"
output_file = "triage.tsv"
min_depth = 1
groups = ["M", "F"]
signif_threshold = 0.05
posterior_threshold = 0.9
bayes_factor_threshold = 10.0

[run.bayes_model]
linkage_prior = 0.01
linked_prevalence = 0.9
null_prevalence = 0.5
group1_linked_weight = 0.5
"#,
        ),
        (
            "freq",
            r#"
[run]
command = "freq"
markers_table = "markers.tsv"
output_file = "freq.tsv"
min_depth = 1
"#,
        ),
        (
            "depth",
            r#"
[run]
command = "depth"
markers_table = "markers.tsv"
popmap = "popmap.tsv"
output_file = "depth.tsv"
min_frequency = 0.75
"#,
        ),
        (
            "map",
            r#"
[run]
command = "map"
markers_file = "markers.tsv"
output_file = "map.tsv"
popmap = "popmap.tsv"
genome_file = "genome.fa"
min_depth = 1
groups = ["M", "F"]
min_quality = 20
min_frequency = 0.1
signif_threshold = 0.05
disable_correction = false
"#,
        ),
        (
            "subset",
            r#"
[run]
command = "subset"
markers_table = "markers.tsv"
popmap = "popmap.tsv"
output_file = "subset.tsv"
min_depth = 1
groups = ["M", "F"]
signif_threshold = 0.05
disable_correction = false
output_fasta = false
min_group1 = 0
min_group2 = 0
max_group1 = 9999
max_group2 = 9999
min_individuals = 0
max_individuals = 9999
"#,
        ),
        (
            "merge",
            r#"
[run]
command = "merge"
input_files = ["a.tsv", "b.tsv"]
output_file = "merged.tsv"
buffer_size = 2000000
output_parquet = false
"#,
        ),
        (
            "pca",
            r#"
[run]
command = "pca"
markers_table = "markers.tsv"
output_dir = "pca"
min_depth = 1
components = 3
"#,
        ),
    ];

    for (expected, body) in commands {
        let profile = RunProfile::parse_toml(&format!("{PROFILE_HEADER}{body}"))
            .unwrap_or_else(|error| panic!("{expected} profile did not parse: {error}"));
        assert_eq!(profile.command_name(), expected);

        let encoded = profile.to_toml().unwrap();
        let decoded = RunProfile::parse_toml(&encoded).unwrap();
        assert_eq!(decoded, profile);
    }
}

#[test]
fn run_profile_rejects_unknown_command_fields() {
    let input = format!(
        "{PROFILE_HEADER}\n[run]\ncommand = \"freq\"\nmarkers_table = \"m.tsv\"\noutput_file = \"f.tsv\"\nmin_depth = 1\nhidden_default = 7\n"
    );
    let error = RunProfile::parse_toml(&input).unwrap_err();
    assert!(error.to_string().contains("hidden_default"));
}

#[test]
fn run_profile_hydrates_explicit_beta_priors_for_bayesian_commands() {
    let input = format!(
        r#"{PROFILE_HEADER}
[run]
command = "triage"
markers_table = "markers.tsv"
popmap = "popmap.tsv"
output_file = "triage.tsv"
min_depth = 1
groups = ["M", "F"]
signif_threshold = 0.05
posterior_threshold = 0.9
bayes_factor_threshold = 10.0

[run.bayes_model]
linkage_prior = 0.01
linked_prevalence = 0.9
null_prevalence = 0.5
group1_linked_weight = 0.5

[run.bayes_model.bayes_factor.alternative_group1]
alpha = 8.0
beta = 2.0

[run.bayes_model.bayes_factor.alternative_group2]
alpha = 2.0
beta = 8.0

[run.bayes_model.bayes_factor.null]
alpha = 10.0
beta = 10.0
"#
    );
    let profile = RunProfile::parse_toml(&input).unwrap();
    let model = match &profile.run {
        CommandProfile::Triage(run) => run.bayes_model.to_runtime().unwrap(),
        _ => panic!("triage profile changed command variants"),
    };
    let configured = bayes_factor_2x2_with_model(9, 1, 10, 10, &model.bayes_factor).unwrap();
    assert!(configured > bayes_factor_2x2(9, 1, 10, 10));
}

#[test]
fn legacy_directional_profile_serializes_its_uniform_prior_defaults() {
    let body = r#"
[run]
command = "triage"
markers_table = "markers.tsv"
popmap = "popmap.tsv"
output_file = "triage.tsv"
min_depth = 1
groups = ["M", "F"]
signif_threshold = 0.05
posterior_threshold = 0.9
bayes_factor_threshold = 10.0

[run.bayes_model]
linkage_prior = 0.01
linked_prevalence = 0.9
null_prevalence = 0.5
group1_linked_weight = 0.5
"#;
    let profile = RunProfile::parse_toml(&format!("{PROFILE_HEADER}{body}")).unwrap();
    let hydrated = profile.to_toml().unwrap();

    assert!(hydrated.contains("[run.bayes_model.bayes_factor.alternative_group1]"));
    assert!(hydrated.contains("alpha = 1.0"));
    assert!(hydrated.contains("beta = 1.0"));
}
