// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

#[cfg(feature = "cuda")]
use rsx_core::compute_backend::compute_chi_squared_batch_with_metrics;
use rsx_core::compute_backend::{AssociationCounts, PValueBackend, compute_chi_squared_batch};

#[test]
fn backend_names_are_explicit() {
    assert_eq!(PValueBackend::parse_str("cpu").unwrap(), PValueBackend::Cpu);
    assert_eq!(
        PValueBackend::parse_str("cuda").unwrap(),
        PValueBackend::Cuda
    );
    assert!(PValueBackend::parse_str("auto").is_err());
}

#[test]
fn cpu_batch_matches_scalar_reference() {
    let total_group1 = 8;
    let total_group2 = 7;
    let counts: Vec<_> = (0..=total_group1)
        .flat_map(|group1| {
            (0..=total_group2).map(move |group2| AssociationCounts { group1, group2 })
        })
        .collect();

    let batch =
        compute_chi_squared_batch(PValueBackend::Cpu, &counts, total_group1, total_group2).unwrap();

    for (counts, observed) in counts.iter().zip(batch) {
        let expected = rsx_core::stats::p_association(
            counts.group1,
            counts.group2,
            total_group1,
            total_group2,
        );
        assert_eq!(observed, expected);
    }
}

#[cfg(feature = "cuda")]
#[test]
fn cuda_batch_matches_cpu_reference() {
    let total_group1 = 16;
    let total_group2 = 13;
    let counts: Vec<_> = (0..=total_group1)
        .flat_map(|group1| {
            (0..=total_group2).map(move |group2| AssociationCounts { group1, group2 })
        })
        .cycle()
        .take(65_537)
        .collect();

    let cpu =
        compute_chi_squared_batch(PValueBackend::Cpu, &counts, total_group1, total_group2).unwrap();
    let cuda = compute_chi_squared_batch(PValueBackend::Cuda, &counts, total_group1, total_group2)
        .unwrap();

    for (index, (expected, observed)) in cpu.iter().zip(cuda).enumerate() {
        let error = (expected - observed).abs();
        assert!(
            error <= 2.0e-15,
            "marker {index}: CPU={expected:.17e}, CUDA={observed:.17e}, abs_error={error:.3e}"
        );
    }
}

#[cfg(feature = "cuda")]
#[test]
fn cuda_reuses_compiled_kernel() {
    let counts = [AssociationCounts {
        group1: 8,
        group2: 2,
    }];
    compute_chi_squared_batch_with_metrics(PValueBackend::Cuda, &counts, 10, 10).unwrap();
    let repeated =
        compute_chi_squared_batch_with_metrics(PValueBackend::Cuda, &counts, 10, 10).unwrap();

    assert_eq!(repeated.metrics.setup_seconds, 0.0);
    assert!(repeated.p_values.is_page_locked());
    assert!(repeated.metrics.output_buffer_reused);
    assert_eq!(repeated.metrics.host_staging_bytes, 0);
}

fn evidence_models() -> Vec<(&'static str, rsx_core::stats::DirectionalModel)> {
    use rsx_core::stats::{BetaPrior, DirectionalModel, PosteriorModel, PrevalencePrior};

    let screening = DirectionalModel::directional_screening_v1();
    let mut beta_posterior = screening;
    beta_posterior.posterior = PosteriorModel {
        linked: PrevalencePrior::Beta(BetaPrior {
            alpha: 8.0,
            beta: 2.0,
        }),
        null: PrevalencePrior::Beta(BetaPrior {
            alpha: 2.0,
            beta: 2.0,
        }),
    };
    let mut skewed = screening;
    skewed.group1_linked_weight = 0.25;
    skewed.linkage_prior = 0.2;
    skewed.bayes_factor.alternative_group1 = BetaPrior {
        alpha: 3.0,
        beta: 5.0,
    };
    skewed.bayes_factor.null = BetaPrior {
        alpha: 0.5,
        beta: 1.5,
    };

    vec![
        ("fixed screening", screening),
        ("beta prevalence", beta_posterior),
        ("skewed directional", skewed),
    ]
}

#[test]
fn cpu_bayes_evidence_matches_the_scalar_path() {
    use rsx_core::compute_backend::compute_bayes_evidence_batch;

    let total_group1 = 9;
    let total_group2 = 6;
    let counts: Vec<_> = (0..=total_group1)
        .flat_map(|group1| {
            (0..=total_group2).map(move |group2| AssociationCounts { group1, group2 })
        })
        .collect();

    for (label, model) in evidence_models() {
        let (factors, posteriors) = compute_bayes_evidence_batch(
            PValueBackend::Cpu,
            &counts,
            total_group1,
            total_group2,
            &model,
        )
        .unwrap();

        for (index, entry) in counts.iter().enumerate() {
            let expected_posterior = rsx_core::stats::posterior_sex_linked_with_model(
                entry.group1,
                entry.group2,
                total_group1,
                total_group2,
                &model,
            );
            assert_eq!(
                posteriors[index], expected_posterior,
                "{label} marker {index}"
            );
            assert!(
                factors[index].is_finite() && factors[index] >= 0.0,
                "{label} marker {index}: bayes factor {}",
                factors[index]
            );
        }
    }
}

#[cfg(feature = "cuda")]
#[test]
fn cuda_bayes_evidence_matches_cpu_reference() {
    use rsx_core::compute_backend::compute_bayes_evidence_batch;

    let total_group1 = 16;
    let total_group2 = 13;
    let counts: Vec<_> = (0..=total_group1)
        .flat_map(|group1| {
            (0..=total_group2).map(move |group2| AssociationCounts { group1, group2 })
        })
        .cycle()
        .take(65_537)
        .collect();

    for (label, model) in evidence_models() {
        let (cpu_factors, cpu_posteriors) = compute_bayes_evidence_batch(
            PValueBackend::Cpu,
            &counts,
            total_group1,
            total_group2,
            &model,
        )
        .unwrap();
        let (cuda_factors, cuda_posteriors) = compute_bayes_evidence_batch(
            PValueBackend::Cuda,
            &counts,
            total_group1,
            total_group2,
            &model,
        )
        .unwrap();

        for index in 0..counts.len() {
            let posterior_error = (cpu_posteriors[index] - cuda_posteriors[index]).abs();
            assert!(
                posterior_error <= 1.0e-12,
                "{label} marker {index}: posterior CPU={:.17e} CUDA={:.17e} abs_error={posterior_error:.3e}",
                cpu_posteriors[index],
                cuda_posteriors[index]
            );

            let expected = cpu_factors[index];
            let observed = cuda_factors[index];
            let relative = if expected == 0.0 {
                observed.abs()
            } else {
                (expected - observed).abs() / expected.abs()
            };
            assert!(
                relative <= 1.0e-12,
                "{label} marker {index}: bayes factor CPU={expected:.17e} CUDA={observed:.17e} rel_error={relative:.3e}"
            );
        }
    }
}

#[cfg(feature = "cuda")]
#[test]
fn cuda_bayes_evidence_shares_the_compiled_module() {
    use rsx_core::compute_backend::compute_bayes_evidence_batch_with_metrics;

    let counts = [AssociationCounts {
        group1: 8,
        group2: 2,
    }];
    let model = rsx_core::stats::DirectionalModel::directional_screening_v1();
    compute_bayes_evidence_batch_with_metrics(PValueBackend::Cuda, &counts, 10, 10, &model)
        .unwrap();
    let repeated =
        compute_bayes_evidence_batch_with_metrics(PValueBackend::Cuda, &counts, 10, 10, &model)
            .unwrap();

    assert_eq!(repeated.metrics.setup_seconds, 0.0);
    assert_eq!(repeated.metrics.markers, 1);
    assert!(repeated.metrics.kernel_seconds > 0.0);
}

#[cfg(feature = "cuda")]
#[test]
fn cuda_gram_accumulation_matches_cpu_reference() {
    use rsx_core::compute_backend::GramAccumulator;

    let individuals = 23usize;
    let markers = 40_001usize;
    // A deterministic spread of depths, including the zeros real tables carry.
    let rows: Vec<Vec<u16>> = (0..markers)
        .map(|marker| {
            (0..individuals)
                .map(|individual| {
                    let value = (marker * 7 + individual * 13) % 41;
                    if value < 9 { 0 } else { value as u16 }
                })
                .collect()
        })
        .collect();

    let mut cpu = GramAccumulator::new(PValueBackend::Cpu, individuals).unwrap();
    let mut cuda = GramAccumulator::new(PValueBackend::Cuda, individuals).unwrap();
    for row in &rows {
        cpu.push(row).unwrap();
        cuda.push(row).unwrap();
    }
    let (cpu_gram, cpu_mean, cpu_markers) = cpu.finish().unwrap();
    let (cuda_gram, cuda_mean, cuda_markers) = cuda.finish().unwrap();

    assert_eq!(cpu_markers, markers as u64);
    assert_eq!(cuda_markers, cpu_markers);
    assert_eq!(
        cuda_mean, cpu_mean,
        "per-individual sums must match exactly"
    );

    // Both sides sum the same integer products in the same order per entry.
    for i in 0..individuals {
        for j in i..individuals {
            let index = i * individuals + j;
            assert_eq!(
                cuda_gram[index], cpu_gram[index],
                "gram entry ({i},{j}) CPU={} CUDA={}",
                cpu_gram[index], cuda_gram[index]
            );
        }
    }
}

#[cfg(feature = "cuda")]
#[test]
fn cuda_matches_cpu_for_every_association_test() {
    use rsx_core::compute_backend::compute_p_batch;
    use rsx_core::test_method::TestMethod;

    let total_group1 = 17;
    let total_group2 = 14;
    let counts: Vec<_> = (0..=total_group1)
        .flat_map(|group1| {
            (0..=total_group2).map(move |group2| AssociationCounts { group1, group2 })
        })
        .cycle()
        .take(20_003)
        .collect();

    for test in [
        TestMethod::ChiSquared,
        TestMethod::Fisher,
        TestMethod::GTest,
    ] {
        let cpu = compute_p_batch(
            PValueBackend::Cpu,
            test,
            &counts,
            total_group1,
            total_group2,
        )
        .unwrap();
        let cuda = compute_p_batch(
            PValueBackend::Cuda,
            test,
            &counts,
            total_group1,
            total_group2,
        )
        .unwrap();

        let mut worst = 0.0f64;
        for (index, (expected, observed)) in cpu.iter().zip(&cuda).enumerate() {
            let error = (expected - observed).abs();
            worst = worst.max(error);
            assert!(
                error <= 1.0e-12,
                "{test:?} marker {index}: CPU={expected:.17e} CUDA={observed:.17e} abs_error={error:.3e}"
            );
        }
        println!("{test:?} worst absolute difference {worst:.3e}");
    }
}
