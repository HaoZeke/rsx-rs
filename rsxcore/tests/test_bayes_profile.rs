// GPL-3.0-or-later

use rsx_core::bayes_profile::{BayesProfileInput, ParameterSource, ProfileOverrides};
use rsx_core::stats::{DirectionalModel, posterior_sex_linked_with_model};

const COMPLETE_PROFILE: &str = r#"
schema_version = 1
profile_name = "directional-screening-v1"
model_id = "directional-binomial-mixture"
model_version = 1

[observation]
group1 = "M"
group2 = "F"
min_depth = 1

[model]
linkage_prior = 0.01
linked_prevalence = 0.90
null_prevalence = 0.50
group1_linked_weight = 0.50

[estimation]
linkage_prior = "fixed"
max_iterations = 100
tolerance = 1.0e-8
minimum_probability = 1.0e-12
maximum_probability = 0.999999999999

[decision]
posterior_threshold = 0.90
bayes_factor_threshold = 10.0
significance_threshold = 0.05
correction = "bonferroni"

[execution]
test = "chisq-yates"
pvalue_backend = "cpu"
"#;

#[test]
fn complete_profile_hydrates_with_stable_digest() {
    let input = BayesProfileInput::parse_toml(COMPLETE_PROFILE).unwrap();
    let hydrated = input.hydrate(&ProfileOverrides::default()).unwrap();

    assert_eq!(hydrated.profile.schema_version, 1);
    assert_eq!(hydrated.profile.model.linkage_prior, 0.01);
    assert_eq!(hydrated.profile.model.linked_prevalence, 0.9);
    assert_eq!(
        hydrated.source_for("model.linked_prevalence"),
        Some(ParameterSource::InputProfile)
    );
    assert_eq!(hydrated.digest_sha256().len(), 64);

    let encoded = hydrated.to_toml().unwrap();
    let decoded = BayesProfileInput::parse_toml(&encoded)
        .unwrap()
        .hydrate(&ProfileOverrides::default())
        .unwrap();
    assert_eq!(decoded.digest_sha256(), hydrated.digest_sha256());
}

#[test]
fn directional_posterior_uses_null_prevalence_and_direction_weight() {
    let group1_model = DirectionalModel {
        linkage_prior: 0.5,
        linked_prevalence: 0.9,
        null_prevalence: 0.3,
        group1_linked_weight: 0.99,
        ..DirectionalModel::directional_screening_v1()
    };
    let group2_model = DirectionalModel {
        group1_linked_weight: 0.01,
        ..group1_model
    };

    let group1_posterior = posterior_sex_linked_with_model(10, 0, 10, 10, &group1_model);
    let group2_posterior = posterior_sex_linked_with_model(10, 0, 10, 10, &group2_model);

    assert!(group1_posterior > group2_posterior);
    assert!(group1_posterior > 0.9);
}

#[test]
fn posterior_prevalence_families_hydrate_and_change_the_calculation() {
    let profile_toml = COMPLETE_PROFILE.replace(
        "[estimation]",
        r#"
[model.posterior.linked]
family = "beta"
alpha = 9.0
beta = 1.0

[model.posterior.null]
family = "beta"
alpha = 5.0
beta = 5.0

[estimation]"#,
    );
    let hydrated = BayesProfileInput::parse_toml(&profile_toml)
        .unwrap()
        .hydrate(&ProfileOverrides::default())
        .unwrap();
    assert_eq!(
        hydrated.source_for("model.posterior.linked.alpha"),
        Some(ParameterSource::InputProfile)
    );

    let beta_model = hydrated.profile.model.to_runtime().unwrap();
    let beta_posterior = posterior_sex_linked_with_model(8, 2, 10, 10, &beta_model);
    let fixed_posterior = posterior_sex_linked_with_model(
        8,
        2,
        10,
        10,
        &DirectionalModel {
            linkage_prior: beta_model.linkage_prior,
            linked_prevalence: beta_model.linked_prevalence,
            null_prevalence: beta_model.null_prevalence,
            group1_linked_weight: beta_model.group1_linked_weight,
            ..DirectionalModel::directional_screening_v1()
        },
    );

    assert_ne!(beta_posterior, fixed_posterior);
    let encoded = hydrated.to_toml().unwrap();
    assert!(encoded.contains("family = \"beta\""));
}
