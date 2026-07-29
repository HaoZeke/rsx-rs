// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Measure the CPU and CUDA paths for the marker evidence and Gram kernels.
//!
//! Usage: `benchmark_bayes_and_gram [markers,...] [individuals] [repetitions]`.

use std::time::Instant;

use rsx_core::compute_backend::{
    AssociationCounts, GramAccumulator, PValueBackend, compute_bayes_evidence_batch,
};
use rsx_core::stats::DirectionalModel;

fn argument<T: std::str::FromStr>(position: usize, fallback: &str) -> Result<T, T::Err> {
    std::env::args()
        .nth(position)
        .unwrap_or_else(|| fallback.to_owned())
        .parse()
}

fn median(mut values: Vec<f64>) -> f64 {
    values.sort_by(f64::total_cmp);
    let middle = values.len() / 2;
    if values.len() % 2 == 1 {
        values[middle]
    } else {
        (values[middle - 1] + values[middle]) / 2.0
    }
}

fn time<T>(mut action: impl FnMut() -> T) -> (T, f64) {
    let started = Instant::now();
    let value = action();
    (value, started.elapsed().as_secs_f64())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let sizes: Vec<usize> = argument::<String>(1, "10000,100000,1000000,10000000")?
        .split(',')
        .map(|item| item.parse::<usize>())
        .collect::<Result<_, _>>()?;
    let individuals: usize = argument(2, "48")?;
    let repetitions: usize = argument(3, "5")?;

    let total_group1 = (individuals / 2) as u32;
    let total_group2 = (individuals - individuals / 2) as u32;
    let model = DirectionalModel::directional_screening_v1();

    println!("kernel,markers,individuals,cpu_median_s,cuda_median_s,speedup");
    for &markers in &sizes {
        let counts: Vec<_> = (0..markers)
            .map(|index| AssociationCounts {
                group1: (index as u32).wrapping_mul(17) % (total_group1 + 1),
                group2: (index as u32).wrapping_mul(29) % (total_group2 + 1),
            })
            .collect();

        let mut cpu = Vec::new();
        let mut cuda = Vec::new();
        for _ in 0..repetitions {
            cpu.push(
                time(|| {
                    compute_bayes_evidence_batch(
                        PValueBackend::Cpu,
                        &counts,
                        total_group1,
                        total_group2,
                        &model,
                    )
                })
                .1,
            );
            cuda.push(
                time(|| {
                    compute_bayes_evidence_batch(
                        PValueBackend::Cuda,
                        &counts,
                        total_group1,
                        total_group2,
                        &model,
                    )
                })
                .1,
            );
        }
        let (cpu_median, cuda_median) = (median(cpu), median(cuda));
        println!(
            "bayes_evidence,{markers},{individuals},{cpu_median:.6},{cuda_median:.6},{:.2}",
            cpu_median / cuda_median
        );

        let rows: Vec<Vec<u16>> = (0..markers.min(2_000_000))
            .map(|marker| {
                (0..individuals)
                    .map(|individual| ((marker * 7 + individual * 13) % 41) as u16)
                    .collect()
            })
            .collect();
        let mut cpu = Vec::new();
        let mut cuda = Vec::new();
        for _ in 0..repetitions {
            for backend in [PValueBackend::Cpu, PValueBackend::Cuda] {
                let (_, seconds) = time(|| -> Result<(), Box<dyn std::error::Error>> {
                    let mut accumulator = GramAccumulator::new(backend, individuals)?;
                    for row in &rows {
                        accumulator.push(row)?;
                    }
                    accumulator.finish()?;
                    Ok(())
                });
                match backend {
                    PValueBackend::Cpu => cpu.push(seconds),
                    PValueBackend::Cuda => cuda.push(seconds),
                }
            }
        }
        let (cpu_median, cuda_median) = (median(cpu), median(cuda));
        println!(
            "gram,{},{individuals},{cpu_median:.6},{cuda_median:.6},{:.2}",
            rows.len(),
            cpu_median / cuda_median
        );
    }
    Ok(())
}
