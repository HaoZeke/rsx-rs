// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Python bindings for rsx via PyO3.

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;

/// Process demultiplexed reads into a marker depth table.
#[pyfunction]
#[pyo3(signature = (input_dir, output_file, threads=1, min_depth=1, kmer_dedup=None))]
fn process(
    input_dir: &str,
    output_file: &str,
    threads: u32,
    min_depth: u16,
    kmer_dedup: Option<usize>,
) -> PyResult<()> {
    rsx_core::commands::process::run(&rsx_core::commands::process::ProcessParams {
        input_dir_path: input_dir.to_string(),
        output_file_path: output_file.to_string(),
        n_threads: threads,
        min_depth,
        kmer_dedup,
    })
    .map_err(|e| PyRuntimeError::new_err(e.to_string()))
}

/// Compute marker distribution between two groups.
#[pyfunction]
#[pyo3(signature = (table_path, popmap_path, output_file, min_depth=1, signif_threshold=0.05, group1="", group2="", correction="bonferroni", test="chisq"))]
#[allow(clippy::too_many_arguments)]
fn distrib(
    table_path: &str,
    popmap_path: &str,
    output_file: &str,
    min_depth: u16,
    signif_threshold: f32,
    group1: &str,
    group2: &str,
    correction: &str,
    test: &str,
) -> PyResult<()> {
    let corr = rsx_core::test_method::CorrectionMethod::parse_str(correction)
        .map_err(PyRuntimeError::new_err)?;
    let tm = rsx_core::test_method::TestMethod::parse_str(test).map_err(PyRuntimeError::new_err)?;

    rsx_core::commands::distrib::run(&rsx_core::commands::distrib::DistribParams {
        markers_table_path: table_path.to_string(),
        popmap_file_path: popmap_path.to_string(),
        output_file_path: output_file.to_string(),
        min_depth,
        signif_threshold,
        correction: corr,
        test_method: tm,
        output_bayes: false,
        group1: group1.to_string(),
        group2: group2.to_string(),
    })
    .map_err(|e| PyRuntimeError::new_err(e.to_string()))
}

/// Extract markers significantly associated with a group.
#[pyfunction]
#[pyo3(signature = (table_path, popmap_path, output_file, min_depth=1, signif_threshold=0.05, group1="", group2="", correction="bonferroni", test="chisq", output_fasta=false, bayes=false))]
#[allow(clippy::too_many_arguments)]
fn signif(
    table_path: &str,
    popmap_path: &str,
    output_file: &str,
    min_depth: u16,
    signif_threshold: f32,
    group1: &str,
    group2: &str,
    correction: &str,
    test: &str,
    output_fasta: bool,
    bayes: bool,
) -> PyResult<()> {
    let corr = rsx_core::test_method::CorrectionMethod::parse_str(correction)
        .map_err(PyRuntimeError::new_err)?;
    let tm = rsx_core::test_method::TestMethod::parse_str(test).map_err(PyRuntimeError::new_err)?;

    rsx_core::commands::signif::run(&rsx_core::commands::signif::SignifParams {
        markers_table_path: table_path.to_string(),
        popmap_file_path: popmap_path.to_string(),
        output_file_path: output_file.to_string(),
        min_depth,
        signif_threshold,
        correction: corr,
        test_method: tm,
        output_fasta,
        output_bayes: bayes,
        group1: group1.to_string(),
        group2: group2.to_string(),
    })
    .map_err(|e| PyRuntimeError::new_err(e.to_string()))
}

/// Triage strict and Bayesian marker candidates.
#[pyfunction]
#[pyo3(signature = (table_path, popmap_path, output_file, min_depth=1, signif_threshold=0.05, posterior_threshold=0.9, bayes_factor_threshold=10.0, prior_probability=0.01, linked_probability=0.9, group1="", group2=""))]
#[allow(clippy::too_many_arguments)]
fn triage(
    table_path: &str,
    popmap_path: &str,
    output_file: &str,
    min_depth: u16,
    signif_threshold: f32,
    posterior_threshold: f64,
    bayes_factor_threshold: f64,
    prior_probability: f64,
    linked_probability: f64,
    group1: &str,
    group2: &str,
) -> PyResult<()> {
    rsx_core::commands::triage::run(&rsx_core::commands::triage::TriageParams {
        markers_table_path: table_path.to_string(),
        popmap_file_path: popmap_path.to_string(),
        output_file_path: output_file.to_string(),
        min_depth,
        signif_threshold,
        posterior_threshold,
        bayes_factor_threshold,
        prior_probability,
        linked_probability,
        group1: group1.to_string(),
        group2: group2.to_string(),
    })
    .map_err(|e| PyRuntimeError::new_err(e.to_string()))
}

/// Compute marker frequency distribution.
#[pyfunction]
#[pyo3(signature = (table_path, output_file, min_depth=1))]
fn freq(table_path: &str, output_file: &str, min_depth: u16) -> PyResult<()> {
    rsx_core::commands::freq::run(&rsx_core::commands::freq::FreqParams {
        markers_table_path: table_path.to_string(),
        output_file_path: output_file.to_string(),
        min_depth,
    })
    .map_err(|e| PyRuntimeError::new_err(e.to_string()))
}

/// Compute depth statistics per individual.
#[pyfunction]
#[pyo3(signature = (table_path, popmap_path, output_file, min_frequency=0.75))]
fn depth(
    table_path: &str,
    popmap_path: &str,
    output_file: &str,
    min_frequency: f32,
) -> PyResult<()> {
    let file_size = std::fs::metadata(table_path).map(|m| m.len()).unwrap_or(0);
    rsx_core::commands::depth::run(&rsx_core::commands::depth::DepthParams {
        markers_table_path: table_path.to_string(),
        popmap_file_path: popmap_path.to_string(),
        output_file_path: output_file.to_string(),
        min_frequency,
        streaming: file_size > 2_000_000_000,
    })
    .map_err(|e| PyRuntimeError::new_err(e.to_string()))
}

/// Merge multiple marker depth tables.
#[pyfunction]
#[pyo3(signature = (input_files, output_file, buffer_size=2000000, output_parquet=false))]
fn merge(
    input_files: Vec<String>,
    output_file: &str,
    buffer_size: usize,
    output_parquet: bool,
) -> PyResult<()> {
    rsx_core::commands::merge::run(&rsx_core::commands::merge::MergeParams {
        input_files,
        output_file_path: output_file.to_string(),
        buffer_size: Some(buffer_size),
        output_parquet,
    })
    .map_err(|e| PyRuntimeError::new_err(e.to_string()))
}

/// Streaming PCA of the depth matrix.
#[pyfunction]
#[pyo3(signature = (table_path, output_dir, min_depth=1, n_components=None))]
fn pca(
    table_path: &str,
    output_dir: &str,
    min_depth: u16,
    n_components: Option<usize>,
) -> PyResult<()> {
    rsx_core::commands::pca::run(&rsx_core::commands::pca::PcaParams {
        markers_table_path: table_path.to_string(),
        output_dir: output_dir.to_string(),
        min_depth,
        n_components,
    })
    .map_err(|e| PyRuntimeError::new_err(e.to_string()))
}

/// pyrsx: Python bindings for rsx (RAD-seq sex determination toolkit)
#[pymodule]
fn pyrsx(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(process, m)?)?;
    m.add_function(wrap_pyfunction!(distrib, m)?)?;
    m.add_function(wrap_pyfunction!(signif, m)?)?;
    m.add_function(wrap_pyfunction!(triage, m)?)?;
    m.add_function(wrap_pyfunction!(freq, m)?)?;
    m.add_function(wrap_pyfunction!(depth, m)?)?;
    m.add_function(wrap_pyfunction!(merge, m)?)?;
    m.add_function(wrap_pyfunction!(pca, m)?)?;

    // Direct Arrow production (first step toward zero temp files for the high-level API)
    m.add_function(wrap_pyfunction!(triage_to_arrow, m)?)?;
    m.add_function(wrap_pyfunction!(pca_to_arrow, m)?)?;

    Ok(())
}

/// Returns the triage results directly as a pyarrow.Table.
///
/// Calls the real in-memory `run_to_arrow` (pure Rust). The RecordBatch(es)
/// are serialized via Arrow IPC entirely in RAM and reconstructed as a
/// pyarrow.Table on the Python side — zero disk files at any point.
/// pyarrow on the Python side. The temp file is deleted before return,
/// so callers see a direct data-producing API with no visible artifacts.
#[pyfunction]
#[pyo3(signature = (table_path, popmap_path, min_depth=1, posterior_threshold=0.9, prior_probability=0.01, linked_probability=0.9, group1="", group2=""))]
#[allow(clippy::too_many_arguments)]
fn triage_to_arrow(
    py: Python<'_>,
    table_path: &str,
    popmap_path: &str,
    min_depth: u16,
    posterior_threshold: f64,
    prior_probability: f64,
    linked_probability: f64,
    group1: &str,
    group2: &str,
) -> PyResult<PyObject> {
    let params = rsx_core::commands::triage::TriageParams {
        markers_table_path: table_path.to_string(),
        popmap_file_path: popmap_path.to_string(),
        output_file_path: String::new(),
        min_depth,
        signif_threshold: 0.05,
        posterior_threshold,
        bayes_factor_threshold: 10.0,
        prior_probability,
        linked_probability,
        group1: group1.to_string(),
        group2: group2.to_string(),
    };

    let batches = rsx_core::commands::triage::run_to_arrow(&params)
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

    if batches.is_empty() {
        let pyarrow = py.import("pyarrow")?;
        let empty = pyarrow.getattr("Table")?.call_method0("from_batches")?;
        return Ok(empty.into());
    }

    // True zero-disk-temp path: in-memory Arrow IPC → pyarrow Table
    let table = batches_to_pyarrow_table(py, &batches)?;
    Ok(table)
}

/// In-memory Arrow IPC export — zero disk files anywhere.
/// Converts Rust RecordBatch(es) into real pyarrow.Table objects via Arrow IPC in RAM only.
fn batches_to_pyarrow_table(
    py: Python<'_>,
    batches: &[arrow::record_batch::RecordBatch],
) -> PyResult<PyObject> {
    if batches.is_empty() {
        // Cache the pyarrow.Table reference for the (rare) empty case
        let pyarrow = py.import("pyarrow")?;
        return Ok(pyarrow.getattr("Table")?.call_method0("from_batches")?.into());
    }

    let mut buf = Vec::new();
    {
        let schema = batches[0].schema();
        let mut writer = arrow::ipc::writer::StreamWriter::try_new(&mut buf, &schema)
            .map_err(|e| PyRuntimeError::new_err(format!("IPC StreamWriter: {e}")))?;
        for batch in batches {
            writer.write(batch)
                .map_err(|e| PyRuntimeError::new_err(format!("IPC write: {e}")))?;
        }
        writer.finish()
            .map_err(|e| PyRuntimeError::new_err(format!("IPC finish: {e}")))?;
    }

    let py_bytes = pyo3::types::PyBytes::new(py, &buf);

    // One import + getattr sequence per top-level call (acceptable; the heavy
    // work is already done in Rust before we reach this point).
    let pyarrow = py.import("pyarrow")?;
    let ipc = pyarrow.getattr("ipc")?;
    let reader = ipc.call_method1("RecordBatchStreamReader", (py_bytes,))?;
    let all_batches_py = reader.call_method0("read_all")?;
    let table = pyarrow.getattr("Table")?.call_method1("from_batches", (all_batches_py,))?;

    Ok(table.into())
}

/// Returns PCA results as a Python dict of pyarrow.Tables:
/// {
///     "eigenvalues": Table with component, eigenvalue, variance_fraction, cumulative,
///     "loadings": Table with individual + PC1..PCk,
///     "n_markers", "n_individuals", "n_components", "total_variance"
/// }
///
/// Uses the real in-memory `run_to_arrow` in rsx-core (pure IPC handoff).
#[pyfunction]
#[pyo3(signature = (table_path, min_depth=1, n_components=None))]
fn pca_to_arrow(
    py: Python<'_>,
    table_path: &str,
    min_depth: u16,
    n_components: Option<usize>,
) -> PyResult<PyObject> {
    let params = rsx_core::commands::pca::PcaParams {
        markers_table_path: table_path.to_string(),
        output_dir: String::new(),
        min_depth,
        n_components,
    };

    let res = rsx_core::commands::pca::run_to_arrow(&params)
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

    // One-shot IPC: write both logical tables (eigenvalues + loadings) into a
    // single Arrow IPC stream. This gives us one import + one round-trip for PCA.
    let mut buf = Vec::new();
    {
        // First batch (eigenvalues)
        let mut writer = arrow::ipc::writer::StreamWriter::try_new(&mut buf, &res.eigenvalues.schema())
            .map_err(|e| PyRuntimeError::new_err(format!("IPC writer (eigen): {e}")))?;
        writer.write(&res.eigenvalues)
            .map_err(|e| PyRuntimeError::new_err(format!("IPC write eigen: {e}")))?;
        // Second batch (loadings) – different schema is fine in one stream
        writer.write(&res.loadings)
            .map_err(|e| PyRuntimeError::new_err(format!("IPC write loadings: {e}")))?;
        writer.finish()
            .map_err(|e| PyRuntimeError::new_err(format!("IPC finish: {e}")))?;
    }

    let py_bytes = pyo3::types::PyBytes::new(py, &buf);
    let pyarrow = py.import("pyarrow")?;
    let ipc = pyarrow.getattr("ipc")?;
    let reader = ipc.call_method1("RecordBatchStreamReader", (py_bytes,))?;

    // Read the two batches we wrote
    let b0 = reader.call_method0("read_next_batch")?;
    let b1 = reader.call_method0("read_next_batch")?;

    let ev_table = pyarrow.getattr("Table")?.call_method1("from_batches", (vec![b0],))?;
    let ld_table = pyarrow.getattr("Table")?.call_method1("from_batches", (vec![b1],))?;

    let dict = pyo3::types::PyDict::new(py);
    dict.set_item("eigenvalues", ev_table)?;
    dict.set_item("loadings", ld_table)?;
    dict.set_item("n_markers", res.n_markers)?;
    dict.set_item("n_individuals", res.n_individuals)?;
    dict.set_item("n_components", res.n_components)?;
    dict.set_item("total_variance", res.total_variance)?;

    Ok(dict.into())
}
