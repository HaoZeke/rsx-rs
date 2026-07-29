// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Raw per-repetition measurements for every kernel the device backend covers.
//!
//! Usage: `benchmark_cuda_coverage [markers,...] [individuals] [repetitions]`
//!
//! One row per (kernel, size, repetition). Summaries are left to the analysis
//! step so the medians reported anywhere can be recomputed from this output.
//!
//! Agreement is reported both ways because the quantities differ in scale: a
//! p-value and a posterior live in [0, 1] where the absolute difference is the
//! meaningful one, while a Bayes factor is unbounded and only its relative
//! difference is interpretable.
//! Repetition 1 carries CUDA context creation and the runtime kernel
//! compilation; it is recorded rather than dropped so the cost stays visible.

use std::time::Instant;

use rsx_core::compute_backend::{
    AssociationCounts, GramAccumulator, PValueBackend, compute_bayes_evidence_batch_with_metrics,
    compute_p_batch_with_metrics,
};
use rsx_core::stats::DirectionalModel;
use rsx_core::test_method::TestMethod;

fn argument<T: std::str::FromStr>(position: usize, fallback: &str) -> Result<T, T::Err> {
    std::env::args()
        .nth(position)
        .unwrap_or_else(|| fallback.to_owned())
        .parse()
}

/// Worst absolute and worst relative difference between host and device values.
fn agreement(host_values: &[f64], device_values: &[f64]) -> (f64, f64) {
    let mut absolute = 0.0f64;
    let mut relative = 0.0f64;
    for (expected, observed) in host_values.iter().zip(device_values) {
        let difference = (expected - observed).abs();
        absolute = absolute.max(difference);
        if *expected != 0.0 {
            relative = relative.max(difference / expected.abs());
        }
    }
    (absolute, relative)
}

fn host() -> String {
    std::fs::read_to_string("/proc/sys/kernel/hostname")
        .map(|name| name.trim().to_owned())
        .unwrap_or_else(|_| "unknown".to_owned())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let sizes: Vec<usize> = argument::<String>(1, "1000,10000,100000,1000000,10000000")?
        .split(',')
        .map(str::parse::<usize>)
        .collect::<Result<_, _>>()?;
    let individuals: usize = argument(2, "48")?;
    let repetitions: usize = argument(3, "5")?;

    let total_group1 = (individuals / 2) as u32;
    let total_group2 = (individuals - individuals / 2) as u32;
    let model = DirectionalModel::directional_screening_v1();
    let host = host();

    println!(
        "host,device,kernel,markers,individuals,repetition,cpu_s,cuda_setup_s,cuda_kernel_s,cuda_total_s,h2d_bytes,d2h_bytes,total_speedup,max_abs_diff,max_rel_diff"
    );

    for &markers in &sizes {
        let counts: Vec<_> = (0..markers)
            .map(|index| AssociationCounts {
                group1: (index as u32).wrapping_mul(17) % (total_group1 + 1),
                group2: (index as u32).wrapping_mul(29) % (total_group2 + 1),
            })
            .collect();

        for test in [
            TestMethod::ChiSquared,
            TestMethod::Fisher,
            TestMethod::GTest,
        ] {
            for repetition in 1..=repetitions {
                let cpu = compute_p_batch_with_metrics(
                    PValueBackend::Cpu,
                    test,
                    &counts,
                    total_group1,
                    total_group2,
                )?;
                let cuda = compute_p_batch_with_metrics(
                    PValueBackend::Cuda,
                    test,
                    &counts,
                    total_group1,
                    total_group2,
                )?;
                let (difference, relative) =
                    agreement(cpu.p_values.try_as_slice()?, cuda.p_values.try_as_slice()?);
                println!(
                    "{host},{},{test:?},{markers},{individuals},{repetition},{:.9},{:.9},{:.9},{:.9},{},{},{:.4},{difference:.3e},{relative:.3e}",
                    cuda.metrics.device,
                    cpu.metrics.total_seconds,
                    cuda.metrics.setup_seconds,
                    cuda.metrics.kernel_seconds,
                    cuda.metrics.total_seconds,
                    cuda.metrics.host_to_device_bytes,
                    cuda.metrics.device_to_host_bytes,
                    cpu.metrics.total_seconds / cuda.metrics.total_seconds,
                );
            }
        }

        for repetition in 1..=repetitions {
            let cpu = compute_bayes_evidence_batch_with_metrics(
                PValueBackend::Cpu,
                &counts,
                total_group1,
                total_group2,
                &model,
            )?;
            let cuda = compute_bayes_evidence_batch_with_metrics(
                PValueBackend::Cuda,
                &counts,
                total_group1,
                total_group2,
                &model,
            )?;
            for (label, host_values, device_values) in [
                ("Posterior", &cpu.posteriors, &cuda.posteriors),
                ("BayesFactor", &cpu.bayes_factors, &cuda.bayes_factors),
            ] {
                let (difference, relative) = agreement(host_values, device_values);
                println!(
                    "{host},{},{label},{markers},{individuals},{repetition},{:.9},{:.9},{:.9},{:.9},{},{},{:.4},{difference:.3e},{relative:.3e}",
                    cuda.metrics.device,
                    cpu.metrics.total_seconds,
                    cuda.metrics.setup_seconds,
                    cuda.metrics.kernel_seconds,
                    cuda.metrics.total_seconds,
                    cuda.metrics.host_to_device_bytes,
                    cuda.metrics.device_to_host_bytes,
                    cpu.metrics.total_seconds / cuda.metrics.total_seconds,
                );
            }
        }

        // The Gram accumulation streams markers, so it is timed end to end and
        // capped to keep the host reference tractable at the largest sizes.
        let gram_markers = markers.min(2_000_000);
        let rows: Vec<Vec<u16>> = (0..gram_markers)
            .map(|marker| {
                (0..individuals)
                    .map(|individual| ((marker * 7 + individual * 13) % 41) as u16)
                    .collect()
            })
            .collect();
        for repetition in 1..=repetitions {
            let mut seconds = [0.0f64; 2];
            let mut totals = Vec::new();
            for (slot, backend) in [PValueBackend::Cpu, PValueBackend::Cuda].iter().enumerate() {
                let started = Instant::now();
                let mut accumulator = GramAccumulator::new(*backend, individuals)?;
                for row in &rows {
                    accumulator.push(row)?;
                }
                let (gram, _, _) = accumulator.finish()?;
                seconds[slot] = started.elapsed().as_secs_f64();
                totals.push(gram);
            }
            let (difference, relative) = agreement(&totals[0], &totals[1]);
            println!(
                "{host},cuda,Gram,{gram_markers},{individuals},{repetition},{:.9},0.0,0.0,{:.9},{},{},{:.4},{difference:.3e},{relative:.3e}",
                seconds[0],
                seconds[1],
                gram_markers * individuals * 2,
                individuals * individuals * 8,
                seconds[0] / seconds[1],
            );
        }
    }
    Ok(())
}
