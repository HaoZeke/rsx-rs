// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

#[cfg(feature = "cuda")]
use rsx_core::compute_backend::compute_chi_squared_batch_with_metrics;
use rsx_core::compute_backend::{compute_chi_squared_batch, AssociationCounts, PValueBackend};

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
