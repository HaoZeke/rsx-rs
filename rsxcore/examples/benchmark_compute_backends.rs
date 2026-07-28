// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

use rsx_core::compute_backend::{
    AssociationCounts, PValueBackend, compute_chi_squared_batch_with_metrics,
};

fn parse_sizes() -> Result<Vec<usize>, Box<dyn std::error::Error>> {
    let value = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "1000,10000,100000,1000000,10000000".to_string());
    value
        .split(',')
        .map(|item| Ok(item.parse::<usize>()?))
        .collect()
}

fn parse_repetitions() -> Result<usize, Box<dyn std::error::Error>> {
    Ok(std::env::args()
        .nth(2)
        .unwrap_or_else(|| "1".to_string())
        .parse()?)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let total_group1 = 24;
    let total_group2 = 24;
    let sizes = parse_sizes()?;

    for repetition in 1..=parse_repetitions()? {
        println!("# repetition={repetition}");
        println!(
            "markers,cpu_total_s,cuda_setup_s,cuda_h2d_s,cuda_kernel_s,cuda_d2h_s,cuda_total_s,h2d_bytes,d2h_bytes,h2d_gb_s,d2h_gb_s,kernel_speedup,total_speedup,output_buffer_reused,max_abs_error,device"
        );
        for &markers in &sizes {
            let counts: Vec<_> = (0..markers)
                .map(|index| AssociationCounts {
                    group1: (index as u32).wrapping_mul(17) % (total_group1 + 1),
                    group2: (index as u32).wrapping_mul(29) % (total_group2 + 1),
                })
                .collect();
            let cpu = compute_chi_squared_batch_with_metrics(
                PValueBackend::Cpu,
                &counts,
                total_group1,
                total_group2,
            )?;
            let cuda = compute_chi_squared_batch_with_metrics(
                PValueBackend::Cuda,
                &counts,
                total_group1,
                total_group2,
            )?;
            let max_abs_error = cpu
                .p_values
                .try_as_slice()?
                .iter()
                .zip(cuda.p_values.try_as_slice()?)
                .map(|(expected, observed)| (expected - observed).abs())
                .fold(0.0f64, f64::max);
            if max_abs_error > 2.0e-15 {
                return Err(format!(
                "CPU/CUDA p-value mismatch for {markers} markers: max_abs_error={max_abs_error:.3e}"
            )
            .into());
            }

            let cpu_metrics = cpu.metrics;
            let cuda_metrics = cuda.metrics;
            let h2d_gb_s = cuda_metrics.host_to_device_bytes as f64
                / cuda_metrics.host_to_device_seconds
                / 1.0e9;
            let d2h_gb_s = cuda_metrics.device_to_host_bytes as f64
                / cuda_metrics.device_to_host_seconds
                / 1.0e9;
            let kernel_speedup = cpu_metrics.kernel_seconds / cuda_metrics.kernel_seconds;
            let total_speedup = cpu_metrics.total_seconds / cuda_metrics.total_seconds;
            println!(
                "{},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{},{},{:.6},{:.6},{:.6},{:.6},{},{:.3e},{}",
                markers,
                cpu_metrics.total_seconds,
                cuda_metrics.setup_seconds,
                cuda_metrics.host_to_device_seconds,
                cuda_metrics.kernel_seconds,
                cuda_metrics.device_to_host_seconds,
                cuda_metrics.total_seconds,
                cuda_metrics.host_to_device_bytes,
                cuda_metrics.device_to_host_bytes,
                h2d_gb_s,
                d2h_gb_s,
                kernel_speedup,
                total_speedup,
                u8::from(cuda_metrics.output_buffer_reused),
                max_abs_error,
                cuda_metrics.device.replace(',', " "),
            );
        }
    }
    Ok(())
}
