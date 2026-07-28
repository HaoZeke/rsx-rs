// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Batched chi-square p-value evaluation on CPU or CUDA.

use std::time::Instant;

/// Marker-presence counts for the two groups under comparison.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AssociationCounts {
    pub group1: u32,
    pub group2: u32,
}

/// Execution backend for batched chi-square p-values.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum PValueBackend {
    #[default]
    Cpu,
    Cuda,
}

impl PValueBackend {
    pub fn parse_str(value: &str) -> Result<Self, String> {
        match value.to_ascii_lowercase().as_str() {
            "cpu" => Ok(Self::Cpu),
            "cuda" | "gpu" => Ok(Self::Cuda),
            _ => Err(format!(
                "Unknown p-value backend: {value}. Options: cpu, cuda"
            )),
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::Cpu => "cpu",
            Self::Cuda => "cuda",
        }
    }
}

/// Timing and transfer accounting for one batch evaluation.
#[derive(Debug, Clone)]
pub struct BatchMetrics {
    pub backend: PValueBackend,
    pub device: String,
    pub markers: usize,
    pub host_to_device_bytes: usize,
    pub device_to_host_bytes: usize,
    pub setup_seconds: f64,
    pub host_to_device_seconds: f64,
    pub kernel_seconds: f64,
    pub device_to_host_seconds: f64,
    pub total_seconds: f64,
}

/// P-values and backend measurements for one batch.
#[derive(Debug, Clone)]
pub struct BatchResult {
    pub p_values: Vec<f64>,
    pub metrics: BatchMetrics,
}

pub fn compute_chi_squared_batch(
    backend: PValueBackend,
    counts: &[AssociationCounts],
    total_group1: u32,
    total_group2: u32,
) -> Result<Vec<f64>, Box<dyn std::error::Error>> {
    Ok(
        compute_chi_squared_batch_with_metrics(backend, counts, total_group1, total_group2)?
            .p_values,
    )
}

pub fn compute_chi_squared_batch_with_metrics(
    backend: PValueBackend,
    counts: &[AssociationCounts],
    total_group1: u32,
    total_group2: u32,
) -> Result<BatchResult, Box<dyn std::error::Error>> {
    match backend {
        PValueBackend::Cpu => Ok(compute_cpu(counts, total_group1, total_group2)),
        PValueBackend::Cuda => compute_cuda(counts, total_group1, total_group2),
    }
}

fn compute_cpu(counts: &[AssociationCounts], total_group1: u32, total_group2: u32) -> BatchResult {
    let started = Instant::now();
    let p_values = counts
        .iter()
        .map(|counts| {
            crate::stats::p_association(counts.group1, counts.group2, total_group1, total_group2)
        })
        .collect();
    let total_seconds = started.elapsed().as_secs_f64();
    BatchResult {
        p_values,
        metrics: BatchMetrics {
            backend: PValueBackend::Cpu,
            device: "host".to_string(),
            markers: counts.len(),
            host_to_device_bytes: 0,
            device_to_host_bytes: 0,
            setup_seconds: 0.0,
            host_to_device_seconds: 0.0,
            kernel_seconds: total_seconds,
            device_to_host_seconds: 0.0,
            total_seconds,
        },
    }
}

#[cfg(not(feature = "cuda"))]
fn compute_cuda(
    _counts: &[AssociationCounts],
    _total_group1: u32,
    _total_group2: u32,
) -> Result<BatchResult, Box<dyn std::error::Error>> {
    Err(std::io::Error::new(
        std::io::ErrorKind::Unsupported,
        "CUDA backend requested, but rsxcore was built without the `cuda` feature",
    )
    .into())
}

#[cfg(feature = "cuda")]
fn compute_cuda(
    counts: &[AssociationCounts],
    total_group1: u32,
    total_group2: u32,
) -> Result<BatchResult, Box<dyn std::error::Error>> {
    use cudarc::driver::{CudaContext, LaunchArgs, LaunchConfig};
    use cudarc::nvrtc::compile_ptx;

    const KERNEL: &str = r#"
        #include <math.h>
        #include <stdint.h>

        extern "C" __global__ void chi_squared_p_values(
            const uint32_t *group1,
            const uint32_t *group2,
            uint32_t total_group1,
            uint32_t total_group2,
            double *p_values,
            uint64_t length
        ) {
            uint64_t index = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;
            if (index >= length) return;

            uint64_t n1 = group1[index];
            uint64_t n2 = group2[index];
            uint64_t t1 = total_group1;
            uint64_t t2 = total_group2;
            double n = (double)(t1 + t2);
            double present = (double)(n1 + n2);
            double absent = n - present;

            if (t1 == 0 || t2 == 0 || present == 0.0 || absent == 0.0) {
                p_values[index] = 1.0;
                return;
            }

            uint64_t ad = n1 * t2;
            uint64_t bc = n2 * t1;
            double difference = (double)(ad > bc ? ad - bc : bc - ad);
            double yates = fmax(difference - n / 2.0, 0.0);
            double chi_squared = n * yates * yates /
                ((double)t1 * (double)t2 * present * absent);
            double p = chi_squared > 0.0 ? erfc(sqrt(chi_squared / 2.0)) : 1.0;
            p_values[index] = fmax(fmin(p, 1.0), 1.0e-16);
        }
    "#;

    let total_started = Instant::now();
    let setup_started = Instant::now();
    let context = CudaContext::new(0)?;
    let device = context.name()?;
    let stream = context.default_stream();
    let ptx = compile_ptx(KERNEL)?;
    let module = context.load_module(ptx)?;
    let function = module.load_function("chi_squared_p_values")?;
    let setup_seconds = setup_started.elapsed().as_secs_f64();

    if counts.is_empty() {
        return Ok(BatchResult {
            p_values: Vec::new(),
            metrics: BatchMetrics {
                backend: PValueBackend::Cuda,
                device,
                markers: 0,
                host_to_device_bytes: 0,
                device_to_host_bytes: 0,
                setup_seconds,
                host_to_device_seconds: 0.0,
                kernel_seconds: 0.0,
                device_to_host_seconds: 0.0,
                total_seconds: total_started.elapsed().as_secs_f64(),
            },
        });
    }

    let group1: Vec<u32> = counts.iter().map(|counts| counts.group1).collect();
    let group2: Vec<u32> = counts.iter().map(|counts| counts.group2).collect();
    let host_to_device_bytes = counts.len() * 2 * std::mem::size_of::<u32>();
    let device_to_host_bytes = counts.len() * std::mem::size_of::<f64>();

    let transfer_started = Instant::now();
    let device_group1 = stream.clone_htod(&group1)?;
    let device_group2 = stream.clone_htod(&group2)?;
    let mut device_p_values = stream.alloc_zeros::<f64>(counts.len())?;
    stream.synchronize()?;
    let host_to_device_seconds = transfer_started.elapsed().as_secs_f64();

    let kernel_started = Instant::now();
    let length = counts.len() as u64;
    let config = LaunchConfig::for_num_elems(counts.len() as u32);
    unsafe {
        stream
            .launch_builder(&function)
            .arg(&device_group1)
            .arg(&device_group2)
            .arg(&total_group1)
            .arg(&total_group2)
            .arg(&mut device_p_values)
            .arg(&length)
            .launch(config)?;
    }
    stream.synchronize()?;
    let kernel_seconds = kernel_started.elapsed().as_secs_f64();

    let return_started = Instant::now();
    let p_values = stream.clone_dtoh(&device_p_values)?;
    stream.synchronize()?;
    let device_to_host_seconds = return_started.elapsed().as_secs_f64();
    let total_seconds = total_started.elapsed().as_secs_f64();

    log::info!(
        "pvalue backend=cuda device={device:?} markers={} h2d_bytes={} d2h_bytes={} setup_s={setup_seconds:.9} h2d_s={host_to_device_seconds:.9} kernel_s={kernel_seconds:.9} d2h_s={device_to_host_seconds:.9} total_s={total_seconds:.9}",
        counts.len(),
        host_to_device_bytes,
        device_to_host_bytes,
    );

    Ok(BatchResult {
        p_values,
        metrics: BatchMetrics {
            backend: PValueBackend::Cuda,
            device,
            markers: counts.len(),
            host_to_device_bytes,
            device_to_host_bytes,
            setup_seconds,
            host_to_device_seconds,
            kernel_seconds,
            device_to_host_seconds,
            total_seconds,
        },
    })
}
