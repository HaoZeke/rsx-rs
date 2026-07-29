// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Batched marker evidence on CPU or CUDA.
//!
//! Two kernels share one module and one retained context: the Yates
//! chi-squared p-value used by `signif`, and the Bayes factor with directional
//! posterior used by `triage`. Both consume the same per-marker presence
//! counts, so a caller that already staged counts for one can evaluate the
//! other without re-reading the marker table.

#[cfg(feature = "cuda")]
use std::sync::{Mutex, OnceLock};
use std::time::Instant;

#[cfg(feature = "cuda")]
const CUDA_KERNEL_SOURCE: &str = r#"
    struct AssociationCounts {
        unsigned int group1;
        unsigned int group2;
    };

    struct PrevalencePrior {
        unsigned int kind;
        unsigned int padding;
        double first;
        double second;
    };

    struct BetaPrior {
        double alpha;
        double beta;
    };

    struct DirectionalModel {
        double log_prior_odds;
        double log_group1_weight;
        double log_group2_weight;
        struct PrevalencePrior posterior_linked;
        struct PrevalencePrior posterior_null;
        struct BetaPrior bayes_group1;
        struct BetaPrior bayes_group2;
        struct BetaPrior bayes_null;
    };

    __device__ double rsx_log_beta(double alpha, double beta) {
        return lgamma(alpha) + lgamma(beta) - lgamma(alpha + beta);
    }

    __device__ double rsx_log_beta_binomial(
        double k, double n, double alpha, double beta
    ) {
        return rsx_log_beta(k + alpha, (n - k) + beta) - rsx_log_beta(alpha, beta);
    }

    __device__ double rsx_log_prevalence_marginal(
        double k, double n, struct PrevalencePrior prior
    ) {
        if (prior.kind == 0u) {
            return k * log(prior.first) + (n - k) * log(1.0 - prior.first);
        }
        return rsx_log_beta_binomial(k, n, prior.first, prior.second);
    }

    __device__ double rsx_logsumexp2(double a, double b) {
        double maximum = fmax(a, b);
        if (isinf(maximum)) return maximum;
        return maximum + log(exp(a - maximum) + exp(b - maximum));
    }

    extern "C" __global__ void bayes_evidence(
        const struct AssociationCounts *counts,
        unsigned int total_group1,
        unsigned int total_group2,
        struct DirectionalModel model,
        double *bayes_factors,
        double *posteriors,
        unsigned long long length
    ) {
        unsigned long long index =
            (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
        if (index >= length) return;

        double n1 = (double)counts[index].group1;
        double n2 = (double)counts[index].group2;
        double t1 = (double)total_group1;
        double t2 = (double)total_group2;
        double total = t1 + t2;

        double log_alternative =
            rsx_log_beta_binomial(n1, t1, model.bayes_group1.alpha, model.bayes_group1.beta) +
            rsx_log_beta_binomial(n2, t2, model.bayes_group2.alpha, model.bayes_group2.beta);
        double log_null =
            rsx_log_beta_binomial(n1 + n2, total, model.bayes_null.alpha, model.bayes_null.beta);
        bayes_factors[index] = exp(log_alternative - log_null);

        double linked_group1 =
            rsx_log_prevalence_marginal(n1 + (t2 - n2), total, model.posterior_linked);
        double linked_group2 =
            rsx_log_prevalence_marginal((t1 - n1) + n2, total, model.posterior_linked);
        double log_linked = rsx_logsumexp2(
            linked_group1 + model.log_group1_weight,
            linked_group2 + model.log_group2_weight
        );
        double log_null_prevalence =
            rsx_log_prevalence_marginal(n1 + n2, total, model.posterior_null);

        double log_odds = log_linked - log_null_prevalence + model.log_prior_odds;
        if (log_odds > 20.0) {
            posteriors[index] = 1.0;
        } else if (log_odds < -20.0) {
            posteriors[index] = 0.0;
        } else {
            posteriors[index] = 1.0 / (1.0 + exp(-log_odds));
        }
    }

    // Markers are split across the grid's second dimension so the launch is
    // not limited to the few thousand upper-triangle entries. Products of two
    // depths are integers held exactly in binary64, so the partial sums may be
    // combined in any order without changing the result.
    extern "C" __global__ void gram_accumulate(
        const unsigned short *depths,
        unsigned int individuals,
        unsigned long long markers,
        unsigned long long markers_per_chunk,
        double *gram,
        double *mean
    ) {
        unsigned long long pair =
            (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
        unsigned long long n = individuals;
        if (pair >= n * n) return;

        unsigned long long i = pair / n;
        unsigned long long j = pair % n;
        if (j < i) return;

        unsigned long long first = (unsigned long long)blockIdx.y * markers_per_chunk;
        if (first >= markers) return;
        unsigned long long last = first + markers_per_chunk;
        if (last > markers) last = markers;

        double total = 0.0;
        for (unsigned long long m = first; m < last; ++m) {
            const unsigned short *row = depths + m * n;
            total += (double)row[i] * (double)row[j];
        }
        if (total != 0.0) atomicAdd(&gram[i * n + j], total);

        if (i == j) {
            double sum = 0.0;
            for (unsigned long long m = first; m < last; ++m) {
                sum += (double)depths[m * n + i];
            }
            if (sum != 0.0) atomicAdd(&mean[i], sum);
        }
    }

    __device__ double rsx_lfact(double value) {
        return lgamma(value + 1.0);
    }

    __device__ double rsx_log_hypergeometric(
        double a, double b, double c, double d,
        double n, double row1, double col1
    ) {
        double col2 = n - col1;
        double row2 = n - row1;
        return rsx_lfact(row1) + rsx_lfact(row2) + rsx_lfact(col1) + rsx_lfact(col2)
             - rsx_lfact(n) - rsx_lfact(a) - rsx_lfact(b) - rsx_lfact(c) - rsx_lfact(d);
    }

    // Two-sided Fisher by the probability method: sum the hypergeometric
    // density over tables no more likely than the observed one. The tail walk
    // is at most one step per individual, which is why this is the test that
    // gains most from the device.
    extern "C" __global__ void fisher_exact_p_values(
        const struct AssociationCounts *counts,
        unsigned int total_group1,
        unsigned int total_group2,
        double *p_values,
        unsigned long long length
    ) {
        unsigned long long index =
            (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
        if (index >= length) return;

        unsigned int n1 = counts[index].group1;
        unsigned int n2 = counts[index].group2;
        if (n1 > total_group1 || n2 > total_group2) {
            p_values[index] = 1.0;
            return;
        }

        unsigned int a = n1;
        unsigned int b = total_group1 - n1;
        unsigned int c = n2;
        unsigned int d = total_group2 - n2;
        unsigned int n = total_group1 + total_group2;
        unsigned int row1 = a + c;
        unsigned int row2 = n - row1;
        unsigned int col1 = a + b;

        double observed = rsx_log_hypergeometric(a, b, c, d, n, row1, col1);
        unsigned int lowest = col1 > row2 ? col1 - row2 : 0u;
        unsigned int highest = row1 < col1 ? row1 : col1;

        // `maximum` is only read once `any` is set, so it needs no sentinel:
        // nvrtc compiles without headers and has no INFINITY.
        double maximum = 0.0;
        double relative = 0.0;
        bool any = false;
        for (unsigned int ai = lowest; ai <= highest; ++ai) {
            unsigned int bi = col1 - ai;
            unsigned int ci = row1 - ai;
            unsigned int di = row2 - bi;
            double density = rsx_log_hypergeometric(ai, bi, ci, di, n, row1, col1);
            if (density > observed + 1e-9) continue;
            if (!any) {
                maximum = density;
                relative = 1.0;
                any = true;
            } else if (density > maximum) {
                relative = relative * exp(maximum - density) + 1.0;
                maximum = density;
            } else {
                relative += exp(density - maximum);
            }
        }

        double p = any ? exp(maximum + log(relative)) : 0.0;
        p_values[index] = fmax(fmin(p, 1.0), 1.0e-16);
    }

    // Log-likelihood-ratio test against the same one-degree chi-squared tail.
    extern "C" __global__ void g_test_p_values(
        const struct AssociationCounts *counts,
        unsigned int total_group1,
        unsigned int total_group2,
        double *p_values,
        unsigned long long length
    ) {
        unsigned long long index =
            (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
        if (index >= length) return;

        unsigned int n1 = counts[index].group1;
        unsigned int n2 = counts[index].group2;
        if (n1 > total_group1 || n2 > total_group2) {
            p_values[index] = 1.0;
            return;
        }

        double a = (double)n1;
        double b = (double)(total_group1 - n1);
        double c = (double)n2;
        double d = (double)(total_group2 - n2);
        double n = (double)(total_group1 + total_group2);

        double row1 = a + c;
        double row2 = b + d;
        double col1 = a + b;
        double col2 = c + d;

        double observations[4] = {a, b, c, d};
        double rows[4] = {row1, row2, row1, row2};
        double columns[4] = {col1, col1, col2, col2};

        double g = 0.0;
        for (int cell = 0; cell < 4; ++cell) {
            double observed = observations[cell];
            if (observed <= 0.0) continue;
            double expected = rows[cell] * columns[cell] / n;
            if (expected > 0.0) {
                g += observed * log(observed / expected);
            }
        }
        g *= 2.0;

        if (isnan(g) || g <= 0.0) {
            p_values[index] = 1.0;
            return;
        }
        double p = erfc(sqrt(g / 2.0));
        p_values[index] = fmax(fmin(p, 1.0), 1.0e-16);
    }

    extern "C" __global__ void chi_squared_p_values(
        const AssociationCounts *counts,
        unsigned int total_group1,
        unsigned int total_group2,
        double *p_values,
        unsigned long long length
    ) {
        unsigned long long index =
            (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
        if (index >= length) return;

        unsigned long long n1 = counts[index].group1;
        unsigned long long n2 = counts[index].group2;
        unsigned long long t1 = total_group1;
        unsigned long long t2 = total_group2;
        double n = (double)(t1 + t2);
        double present = (double)(n1 + n2);
        double absent = n - present;

        if (t1 == 0 || t2 == 0 || present == 0.0 || absent == 0.0) {
            p_values[index] = 1.0;
            return;
        }

        unsigned long long ad = n1 * t2;
        unsigned long long bc = n2 * t1;
        double difference = (double)(ad > bc ? ad - bc : bc - ad);
        double yates = fmax(difference - n / 2.0, 0.0);
        double chi_squared = n * yates * yates /
            ((double)t1 * (double)t2 * present * absent);
        double p = chi_squared > 0.0 ? erfc(sqrt(chi_squared / 2.0)) : 1.0;
        p_values[index] = fmax(fmin(p, 1.0), 1.0e-16);
    }
"#;

#[cfg(feature = "cuda")]
struct CachedCudaModule {
    context: std::sync::Arc<cudarc::driver::CudaContext>,
    chi_squared: cudarc::driver::CudaFunction,
    fisher_exact: cudarc::driver::CudaFunction,
    g_test: cudarc::driver::CudaFunction,
    bayes_evidence: cudarc::driver::CudaFunction,
    gram_accumulate: cudarc::driver::CudaFunction,
}

#[cfg(feature = "cuda")]
static CUDA_MODULE: OnceLock<CachedCudaModule> = OnceLock::new();

#[cfg(feature = "cuda")]
static CUDA_PINNED_RESULTS: Mutex<Vec<cudarc::driver::PinnedHostSlice<f64>>> =
    Mutex::new(Vec::new());

#[cfg(feature = "cuda")]
#[derive(Debug)]
pub struct PooledPinnedResult {
    values: Option<cudarc::driver::PinnedHostSlice<f64>>,
    len: usize,
}

#[cfg(feature = "cuda")]
impl PooledPinnedResult {
    fn acquire(
        context: &std::sync::Arc<cudarc::driver::CudaContext>,
        len: usize,
    ) -> Result<(Self, bool), Box<dyn std::error::Error>> {
        let mut pool = CUDA_PINNED_RESULTS
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if let Some(index) = pool.iter().position(|values| values.len() >= len) {
            let values = pool.swap_remove(index);
            return Ok((
                Self {
                    values: Some(values),
                    len,
                },
                true,
            ));
        }
        drop(pool);

        let values = unsafe { context.alloc_pinned::<f64>(len)? };
        Ok((
            Self {
                values: Some(values),
                len,
            },
            false,
        ))
    }

    fn try_as_mut_slice(&mut self) -> Result<&mut [f64], cudarc::driver::DriverError> {
        let values = self
            .values
            .as_mut()
            .expect("pooled CUDA result must own its allocation");
        Ok(&mut values.as_mut_slice()?[..self.len])
    }

    fn try_as_slice(&self) -> Result<&[f64], cudarc::driver::DriverError> {
        Ok(&self
            .values
            .as_ref()
            .expect("pooled CUDA result must own its allocation")
            .as_slice()?[..self.len])
    }
}

#[cfg(feature = "cuda")]
impl Drop for PooledPinnedResult {
    fn drop(&mut self) {
        let Some(values) = self.values.take() else {
            return;
        };
        let mut pool = CUDA_PINNED_RESULTS
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if pool.iter().all(|cached| cached.len() < values.len()) {
            pool.clear();
            pool.push(values);
        }
    }
}

#[cfg(feature = "cuda")]
fn cuda_module() -> Result<(&'static CachedCudaModule, f64), Box<dyn std::error::Error>> {
    use cudarc::driver::CudaContext;
    use cudarc::nvrtc::compile_ptx_with_opts;

    if let Some(module) = CUDA_MODULE.get() {
        return Ok((module, 0.0));
    }

    let started = Instant::now();
    let context = CudaContext::new(0)?;

    // Compile for the device that will run it. The Gram kernel adds doubles
    // atomically, which needs compute capability 6.0 or later, and nvrtc's
    // default target is older than that on some toolkits.
    let (major, minor) = context.compute_capability()?;
    let architecture: &'static str = Box::leak(format!("compute_{major}{minor}").into_boxed_str());
    let ptx = compile_ptx_with_opts(
        CUDA_KERNEL_SOURCE,
        cudarc::nvrtc::CompileOptions {
            arch: Some(architecture),
            ..Default::default()
        },
    )?;
    let module = context.load_module(ptx)?;
    let chi_squared = module.load_function("chi_squared_p_values")?;
    let fisher_exact = module.load_function("fisher_exact_p_values")?;
    let g_test = module.load_function("g_test_p_values")?;
    let bayes_evidence = module.load_function("bayes_evidence")?;
    let gram_accumulate = module.load_function("gram_accumulate")?;
    let _ = CUDA_MODULE.set(CachedCudaModule {
        context,
        chi_squared,
        fisher_exact,
        g_test,
        bayes_evidence,
        gram_accumulate,
    });

    Ok((
        CUDA_MODULE
            .get()
            .expect("CUDA module cache must be initialized"),
        started.elapsed().as_secs_f64(),
    ))
}

/// Marker-presence counts for the two groups under comparison.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AssociationCounts {
    pub group1: u32,
    pub group2: u32,
}

#[cfg(feature = "cuda")]
unsafe impl cudarc::driver::DeviceRepr for AssociationCounts {}

/// One prevalence prior in the flat form the kernel reads.
///
/// `kind` 0 carries a fixed probability in `first`; kind 1 carries Beta shapes
/// in `first` and `second`.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DevicePrevalencePrior {
    kind: u32,
    padding: u32,
    first: f64,
    second: f64,
}

impl From<crate::stats::PrevalencePrior> for DevicePrevalencePrior {
    fn from(prior: crate::stats::PrevalencePrior) -> Self {
        match prior {
            crate::stats::PrevalencePrior::Fixed { probability } => Self {
                kind: 0,
                padding: 0,
                first: probability,
                second: 0.0,
            },
            crate::stats::PrevalencePrior::Beta(beta) => Self {
                kind: 1,
                padding: 0,
                first: beta.alpha,
                second: beta.beta,
            },
        }
    }
}

/// Beta shapes for one Bayes-factor hypothesis.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DeviceBetaPrior {
    alpha: f64,
    beta: f64,
}

impl From<crate::stats::BetaPrior> for DeviceBetaPrior {
    fn from(prior: crate::stats::BetaPrior) -> Self {
        Self {
            alpha: prior.alpha,
            beta: prior.beta,
        }
    }
}

/// The directional model as the kernel consumes it.
///
/// The three logarithms are taken on the host from the same expressions the
/// scalar path uses, so the device never recomputes a value that would drift
/// from the CPU result.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DeviceDirectionalModel {
    log_prior_odds: f64,
    log_group1_weight: f64,
    log_group2_weight: f64,
    posterior_linked: DevicePrevalencePrior,
    posterior_null: DevicePrevalencePrior,
    bayes_group1: DeviceBetaPrior,
    bayes_group2: DeviceBetaPrior,
    bayes_null: DeviceBetaPrior,
}

impl From<&crate::stats::DirectionalModel> for DeviceDirectionalModel {
    fn from(model: &crate::stats::DirectionalModel) -> Self {
        Self {
            log_prior_odds: (model.linkage_prior / (1.0 - model.linkage_prior)).ln(),
            log_group1_weight: model.group1_linked_weight.ln(),
            log_group2_weight: (1.0 - model.group1_linked_weight).ln(),
            posterior_linked: model.posterior.linked.into(),
            posterior_null: model.posterior.null.into(),
            bayes_group1: model.bayes_factor.alternative_group1.into(),
            bayes_group2: model.bayes_factor.alternative_group2.into(),
            bayes_null: model.bayes_factor.null.into(),
        }
    }
}

#[cfg(feature = "cuda")]
unsafe impl cudarc::driver::DeviceRepr for DeviceDirectionalModel {}

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
    pub output_buffer_reused: bool,
    pub host_staging_bytes: usize,
}

/// Host storage for computed p-values.
#[derive(Debug)]
pub enum PValueBuffer {
    Owned(Vec<f64>),
    #[cfg(feature = "cuda")]
    PageLocked(PooledPinnedResult),
}

impl PValueBuffer {
    pub fn try_as_slice(&self) -> Result<&[f64], Box<dyn std::error::Error>> {
        match self {
            Self::Owned(values) => Ok(values),
            #[cfg(feature = "cuda")]
            Self::PageLocked(values) => Ok(values.try_as_slice()?),
        }
    }

    pub fn is_page_locked(&self) -> bool {
        match self {
            Self::Owned(_) => false,
            #[cfg(feature = "cuda")]
            Self::PageLocked(_) => true,
        }
    }

    fn into_vec(self) -> Result<Vec<f64>, Box<dyn std::error::Error>> {
        match self {
            Self::Owned(values) => Ok(values),
            #[cfg(feature = "cuda")]
            Self::PageLocked(values) => Ok(values.try_as_slice()?.to_vec()),
        }
    }
}

/// P-values and backend measurements for one batch.
#[derive(Debug)]
pub struct BatchResult {
    pub p_values: PValueBuffer,
    pub metrics: BatchMetrics,
}

pub fn compute_chi_squared_batch(
    backend: PValueBackend,
    counts: &[AssociationCounts],
    total_group1: u32,
    total_group2: u32,
) -> Result<Vec<f64>, Box<dyn std::error::Error>> {
    compute_chi_squared_batch_with_metrics(backend, counts, total_group1, total_group2)?
        .p_values
        .into_vec()
}

pub fn compute_chi_squared_batch_with_metrics(
    backend: PValueBackend,
    counts: &[AssociationCounts],
    total_group1: u32,
    total_group2: u32,
) -> Result<BatchResult, Box<dyn std::error::Error>> {
    compute_p_batch_with_metrics(
        backend,
        crate::test_method::TestMethod::ChiSquared,
        counts,
        total_group1,
        total_group2,
    )
}

/// Batched p-values for any supported association test.
pub fn compute_p_batch(
    backend: PValueBackend,
    test: crate::test_method::TestMethod,
    counts: &[AssociationCounts],
    total_group1: u32,
    total_group2: u32,
) -> Result<Vec<f64>, Box<dyn std::error::Error>> {
    compute_p_batch_with_metrics(backend, test, counts, total_group1, total_group2)?
        .p_values
        .into_vec()
}

pub fn compute_p_batch_with_metrics(
    backend: PValueBackend,
    test: crate::test_method::TestMethod,
    counts: &[AssociationCounts],
    total_group1: u32,
    total_group2: u32,
) -> Result<BatchResult, Box<dyn std::error::Error>> {
    match backend {
        PValueBackend::Cpu => Ok(compute_cpu(test, counts, total_group1, total_group2)),
        PValueBackend::Cuda => compute_cuda(test, counts, total_group1, total_group2),
    }
}

fn compute_cpu(
    test: crate::test_method::TestMethod,
    counts: &[AssociationCounts],
    total_group1: u32,
    total_group2: u32,
) -> BatchResult {
    let started = Instant::now();
    let p_values = counts
        .iter()
        .map(|counts| {
            crate::test_method::compute_p(
                test,
                counts.group1,
                counts.group2,
                total_group1,
                total_group2,
            )
        })
        .collect();
    let total_seconds = started.elapsed().as_secs_f64();
    BatchResult {
        p_values: PValueBuffer::Owned(p_values),
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
            output_buffer_reused: false,
            host_staging_bytes: 0,
        },
    }
}

#[cfg(not(feature = "cuda"))]
fn compute_cuda(
    _test: crate::test_method::TestMethod,
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
    test: crate::test_method::TestMethod,
    counts: &[AssociationCounts],
    total_group1: u32,
    total_group2: u32,
) -> Result<BatchResult, Box<dyn std::error::Error>> {
    use cudarc::driver::PushKernelArg;

    let total_started = Instant::now();
    let (module, setup_seconds) = cuda_module()?;
    let kernel = match test {
        crate::test_method::TestMethod::ChiSquared => &module.chi_squared,
        crate::test_method::TestMethod::Fisher => &module.fisher_exact,
        crate::test_method::TestMethod::GTest => &module.g_test,
    };
    let context = &module.context;
    let device = context.name()?;
    let stream = context.default_stream();

    if counts.is_empty() {
        return Ok(BatchResult {
            p_values: PValueBuffer::Owned(Vec::new()),
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
                output_buffer_reused: false,
                host_staging_bytes: 0,
            },
        });
    }

    let host_to_device_bytes = std::mem::size_of_val(counts);
    let device_to_host_bytes = counts.len() * std::mem::size_of::<f64>();

    let transfer_started = Instant::now();
    let device_counts = stream.clone_htod(counts)?;
    let mut device_p_values = stream.alloc_zeros::<f64>(counts.len())?;
    stream.synchronize()?;
    let host_to_device_seconds = transfer_started.elapsed().as_secs_f64();

    let kernel_started = Instant::now();
    let length = counts.len() as u64;
    let config = marker_launch_config(counts.len());
    unsafe {
        stream
            .launch_builder(kernel)
            .arg(&device_counts)
            .arg(&total_group1)
            .arg(&total_group2)
            .arg(&mut device_p_values)
            .arg(&length)
            .launch(config)?;
    }
    stream.synchronize()?;
    let kernel_seconds = kernel_started.elapsed().as_secs_f64();

    let return_started = Instant::now();
    let (mut p_values, output_buffer_reused) = PooledPinnedResult::acquire(context, counts.len())?;
    stream.memcpy_dtoh(&device_p_values, p_values.try_as_mut_slice()?)?;
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
        p_values: PValueBuffer::PageLocked(p_values),
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
            output_buffer_reused,
            host_staging_bytes: 0,
        },
    })
}

#[cfg(feature = "cuda")]
/// Threads per block for the per-marker kernels.
///
/// The Fisher tail walk needs enough registers that a full 1024-thread block
/// is refused with CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES, so every per-marker
/// kernel launches at a width all three of them can occupy.
const MARKER_KERNEL_BLOCK: u32 = 256;

#[cfg(feature = "cuda")]
fn marker_launch_config(elements: usize) -> cudarc::driver::LaunchConfig {
    cudarc::driver::LaunchConfig {
        grid_dim: ((elements as u32).div_ceil(MARKER_KERNEL_BLOCK).max(1), 1, 1),
        block_dim: (MARKER_KERNEL_BLOCK, 1, 1),
        shared_mem_bytes: 0,
    }
}

/// Per-marker Bayes factors and directional posteriors with backend timings.
#[derive(Debug)]
pub struct BayesEvidenceResult {
    pub bayes_factors: Vec<f64>,
    pub posteriors: Vec<f64>,
    pub metrics: BatchMetrics,
}

/// Evaluate the Bayes factor and directional posterior for a batch of markers.
///
/// The model is validated once by the caller rather than per marker, matching
/// the scalar hot path in `triage`.
pub fn compute_bayes_evidence_batch(
    backend: PValueBackend,
    counts: &[AssociationCounts],
    total_group1: u32,
    total_group2: u32,
    model: &crate::stats::DirectionalModel,
) -> Result<(Vec<f64>, Vec<f64>), Box<dyn std::error::Error>> {
    let result = compute_bayes_evidence_batch_with_metrics(
        backend,
        counts,
        total_group1,
        total_group2,
        model,
    )?;
    Ok((result.bayes_factors, result.posteriors))
}

pub fn compute_bayes_evidence_batch_with_metrics(
    backend: PValueBackend,
    counts: &[AssociationCounts],
    total_group1: u32,
    total_group2: u32,
    model: &crate::stats::DirectionalModel,
) -> Result<BayesEvidenceResult, Box<dyn std::error::Error>> {
    match backend {
        PValueBackend::Cpu => Ok(bayes_evidence_cpu(
            counts,
            total_group1,
            total_group2,
            model,
        )),
        PValueBackend::Cuda => bayes_evidence_cuda(counts, total_group1, total_group2, model),
    }
}

fn bayes_evidence_cpu(
    counts: &[AssociationCounts],
    total_group1: u32,
    total_group2: u32,
    model: &crate::stats::DirectionalModel,
) -> BayesEvidenceResult {
    let started = Instant::now();
    let mut bayes_factors = Vec::with_capacity(counts.len());
    let mut posteriors = Vec::with_capacity(counts.len());
    for entry in counts {
        bayes_factors.push(crate::stats::bayes_factor_2x2_with_validated_model(
            entry.group1,
            entry.group2,
            total_group1,
            total_group2,
            &model.bayes_factor,
        ));
        posteriors.push(crate::stats::posterior_sex_linked_with_model(
            entry.group1,
            entry.group2,
            total_group1,
            total_group2,
            model,
        ));
    }
    let total_seconds = started.elapsed().as_secs_f64();
    BayesEvidenceResult {
        bayes_factors,
        posteriors,
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
            output_buffer_reused: false,
            host_staging_bytes: 0,
        },
    }
}

#[cfg(not(feature = "cuda"))]
fn bayes_evidence_cuda(
    _counts: &[AssociationCounts],
    _total_group1: u32,
    _total_group2: u32,
    _model: &crate::stats::DirectionalModel,
) -> Result<BayesEvidenceResult, Box<dyn std::error::Error>> {
    Err(std::io::Error::new(
        std::io::ErrorKind::Unsupported,
        "CUDA backend requested, but rsxcore was built without the `cuda` feature",
    )
    .into())
}

#[cfg(feature = "cuda")]
fn bayes_evidence_cuda(
    counts: &[AssociationCounts],
    total_group1: u32,
    total_group2: u32,
    model: &crate::stats::DirectionalModel,
) -> Result<BayesEvidenceResult, Box<dyn std::error::Error>> {
    use cudarc::driver::PushKernelArg;

    let total_started = Instant::now();
    let (module, setup_seconds) = cuda_module()?;
    let context = &module.context;
    let device = context.name()?;
    let stream = context.default_stream();

    if counts.is_empty() {
        return Ok(BayesEvidenceResult {
            bayes_factors: Vec::new(),
            posteriors: Vec::new(),
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
                output_buffer_reused: false,
                host_staging_bytes: 0,
            },
        });
    }

    let device_model = DeviceDirectionalModel::from(model);
    let host_to_device_bytes =
        std::mem::size_of_val(counts) + std::mem::size_of::<DeviceDirectionalModel>();
    let device_to_host_bytes = 2 * counts.len() * std::mem::size_of::<f64>();

    let transfer_started = Instant::now();
    let device_counts = stream.clone_htod(counts)?;
    let mut device_bayes_factors = stream.alloc_zeros::<f64>(counts.len())?;
    let mut device_posteriors = stream.alloc_zeros::<f64>(counts.len())?;
    stream.synchronize()?;
    let host_to_device_seconds = transfer_started.elapsed().as_secs_f64();

    let kernel_started = Instant::now();
    let length = counts.len() as u64;
    let config = marker_launch_config(counts.len());
    unsafe {
        stream
            .launch_builder(&module.bayes_evidence)
            .arg(&device_counts)
            .arg(&total_group1)
            .arg(&total_group2)
            .arg(&device_model)
            .arg(&mut device_bayes_factors)
            .arg(&mut device_posteriors)
            .arg(&length)
            .launch(config)?;
    }
    stream.synchronize()?;
    let kernel_seconds = kernel_started.elapsed().as_secs_f64();

    let return_started = Instant::now();
    let bayes_factors = stream.clone_dtoh(&device_bayes_factors)?;
    let posteriors = stream.clone_dtoh(&device_posteriors)?;
    stream.synchronize()?;
    let device_to_host_seconds = return_started.elapsed().as_secs_f64();
    let total_seconds = total_started.elapsed().as_secs_f64();

    log::info!(
        "bayes backend=cuda device={device:?} markers={} h2d_bytes={host_to_device_bytes} d2h_bytes={device_to_host_bytes} setup_s={setup_seconds:.9} h2d_s={host_to_device_seconds:.9} kernel_s={kernel_seconds:.9} d2h_s={device_to_host_seconds:.9} total_s={total_seconds:.9}",
        counts.len(),
    );

    Ok(BayesEvidenceResult {
        bayes_factors,
        posteriors,
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
            output_buffer_reused: false,
            host_staging_bytes: 0,
        },
    })
}

/// Upper-triangle Gram matrix, per-individual depth sums, and marker count.
pub type GramTotals = (Vec<f64>, Vec<f64>, u64);

/// Streaming accumulation of the marker-by-individual Gram matrix.
///
/// Markers arrive one at a time but the device wants many, so the host fills a
/// tile and hands whole tiles over. The CPU variant applies the same rank-1
/// update directly, which keeps PCA on one code path.
pub struct GramAccumulator {
    backend: PValueBackend,
    individuals: usize,
    tile_markers: usize,
    tile: Vec<u16>,
    gram: Vec<f64>,
    mean: Vec<f64>,
    markers: u64,
    #[cfg(feature = "cuda")]
    device: Option<CudaGramState>,
}

#[cfg(feature = "cuda")]
struct CudaGramState {
    gram: cudarc::driver::CudaSlice<f64>,
    mean: cudarc::driver::CudaSlice<f64>,
}

impl GramAccumulator {
    /// Markers per device transfer. Sized so one tile stays a few megabytes for
    /// realistic individual counts.
    const DEFAULT_TILE_MARKERS: usize = 262_144;
    /// Markers each thread walks before the grid adds another chunk.
    const MARKERS_PER_CHUNK: u64 = 512;

    pub fn new(
        backend: PValueBackend,
        individuals: usize,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let mut accumulator = Self {
            backend,
            individuals,
            tile_markers: Self::DEFAULT_TILE_MARKERS,
            tile: Vec::new(),
            gram: vec![0.0; individuals * individuals],
            mean: vec![0.0; individuals],
            markers: 0,
            #[cfg(feature = "cuda")]
            device: None,
        };
        if backend == PValueBackend::Cuda {
            accumulator.prepare_device()?;
        }
        Ok(accumulator)
    }

    #[cfg(feature = "cuda")]
    fn prepare_device(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        let (module, _) = cuda_module()?;
        let stream = module.context.default_stream();
        self.device = Some(CudaGramState {
            gram: stream.alloc_zeros::<f64>(self.individuals * self.individuals)?,
            mean: stream.alloc_zeros::<f64>(self.individuals)?,
        });
        self.tile.reserve(self.tile_markers * self.individuals);
        Ok(())
    }

    #[cfg(not(feature = "cuda"))]
    fn prepare_device(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        Err(std::io::Error::new(
            std::io::ErrorKind::Unsupported,
            "CUDA backend requested, but rsxcore was built without the `cuda` feature",
        )
        .into())
    }

    /// Fold one marker's per-individual depths into the accumulation.
    pub fn push(&mut self, depths: &[u16]) -> Result<(), Box<dyn std::error::Error>> {
        debug_assert_eq!(depths.len(), self.individuals);
        self.markers += 1;
        match self.backend {
            PValueBackend::Cpu => {
                let n = self.individuals;
                for (i, &di) in depths.iter().enumerate() {
                    let xi = f64::from(di);
                    self.mean[i] += xi;
                    let base = i * n;
                    for (j, &dj) in depths.iter().enumerate().skip(i) {
                        self.gram[base + j] += xi * f64::from(dj);
                    }
                }
                Ok(())
            }
            PValueBackend::Cuda => {
                self.tile.extend_from_slice(depths);
                if self.tile.len() >= self.tile_markers * self.individuals {
                    self.flush()?;
                }
                Ok(())
            }
        }
    }

    #[cfg(feature = "cuda")]
    fn flush(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        use cudarc::driver::{LaunchConfig, PushKernelArg};

        if self.tile.is_empty() {
            return Ok(());
        }
        let (module, _) = cuda_module()?;
        let stream = module.context.default_stream();
        let state = self
            .device
            .as_mut()
            .expect("device state exists whenever the CUDA backend is selected");

        let markers = (self.tile.len() / self.individuals) as u64;
        let device_tile = stream.clone_htod(&self.tile)?;
        let individuals = self.individuals as u32;
        let pairs = (self.individuals * self.individuals) as u32;

        // Spread the tile over enough chunks to fill the device: the pair count
        // alone is only a few thousand threads.
        let block = 128u32;
        let chunks = markers.div_ceil(Self::MARKERS_PER_CHUNK).max(1);
        let markers_per_chunk = markers.div_ceil(chunks).max(1);
        let config = LaunchConfig {
            grid_dim: (pairs.div_ceil(block), chunks as u32, 1),
            block_dim: (block, 1, 1),
            shared_mem_bytes: 0,
        };
        unsafe {
            stream
                .launch_builder(&module.gram_accumulate)
                .arg(&device_tile)
                .arg(&individuals)
                .arg(&markers)
                .arg(&markers_per_chunk)
                .arg(&mut state.gram)
                .arg(&mut state.mean)
                .launch(config)?;
        }
        stream.synchronize()?;
        self.tile.clear();
        Ok(())
    }

    #[cfg(not(feature = "cuda"))]
    fn flush(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        Ok(())
    }

    /// Upper-triangle Gram, per-individual sums, and the marker count.
    pub fn finish(mut self) -> Result<GramTotals, Box<dyn std::error::Error>> {
        if self.backend == PValueBackend::Cuda {
            self.flush()?;
            #[cfg(feature = "cuda")]
            {
                let (module, _) = cuda_module()?;
                let stream = module.context.default_stream();
                let state = self
                    .device
                    .as_ref()
                    .expect("device state exists whenever the CUDA backend is selected");
                self.gram = stream.clone_dtoh(&state.gram)?;
                self.mean = stream.clone_dtoh(&state.mean)?;
            }
        }
        Ok((self.gram, self.mean, self.markers))
    }
}
