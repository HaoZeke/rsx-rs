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
/// Calls the real in-memory `run_to_arrow` (pure Rust, identical logic to
/// the classic file writer). The resulting RecordBatch(es) are written to
/// a hidden NamedTempFile using the Arrow Parquet writer, then read by
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
    use std::fs;
    use tempfile::NamedTempFile;

    let params = rsx_core::commands::triage::TriageParams {
        markers_table_path: table_path.to_string(),
        popmap_file_path: popmap_path.to_string(),
        output_file_path: String::new(), // unused — we go through the Arrow path
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
        // Empty result — return a proper 0-row table via Python (schema will be inferred or empty)
        let pyarrow = py.import("pyarrow")?;
        let py_list = pyo3::types::PyList::empty(py);
        let empty_table = pyarrow.getattr("Table")?.call_method1("from_batches", (py_list,))?;
        return Ok(empty_table.into());
    }

    // Write the real Arrow batches to a hidden temp Parquet
    let output_file = NamedTempFile::new()
        .map_err(|e| PyRuntimeError::new_err(format!("failed to create temp file: {e}")))?;
    let output_path = output_file.path().to_string_lossy().to_string();

    {
        let file = std::fs::File::create(&output_path)
            .map_err(|e| PyRuntimeError::new_err(format!("failed to open temp parquet: {e}")))?;
        let mut writer = parquet::arrow::ArrowWriter::try_new(file, batches[0].schema(), None)
            .map_err(|e| PyRuntimeError::new_err(format!("ArrowWriter init failed: {e}")))?;
        for batch in &batches {
            writer.write(batch).map_err(|e| PyRuntimeError::new_err(format!("write batch: {e}")))?;
        }
        writer.close().map_err(|e| PyRuntimeError::new_err(format!("writer close: {e}")))?;
    }

    // Hand the Parquet to pyarrow (fast path) and return a real Table
    let pyarrow = py.import("pyarrow")?;
    let pyarrow_parquet = pyarrow.getattr("parquet")?;
    let table = pyarrow_parquet.call_method1("read_table", (output_path.as_str(),))?;

    // Delete the internal file before returning — caller never sees it
    let _ = fs::remove_file(&output_path);

    Ok(table.into())
}

/// Returns PCA results as a Python dict of pyarrow.Tables:
/// {
///     "eigenvalues": Table with component, eigenvalue, variance_fraction, cumulative,
///     "loadings": Table with individual + PC1..PCk,
///     "n_markers", "n_individuals", "n_components", "total_variance"
/// }
///
/// Uses the real in-memory `run_to_arrow` in rsx-core (identical Jacobi +
/// streaming Gram math). Hidden temp Parquet files are cleaned before return.
#[pyfunction]
#[pyo3(signature = (table_path, min_depth=1, n_components=None))]
fn pca_to_arrow(
    py: Python<'_>,
    table_path: &str,
    min_depth: u16,
    n_components: Option<usize>,
) -> PyResult<PyObject> {
    use std::fs;
    use tempfile::NamedTempFile;

    let params = rsx_core::commands::pca::PcaParams {
        markers_table_path: table_path.to_string(),
        output_dir: String::new(), // unused by Arrow path
        min_depth,
        n_components,
    };

    let res = rsx_core::commands::pca::run_to_arrow(&params)
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

    // Helper to write one batch to a hidden temp Parquet and return the path
    fn write_batch_to_hidden_parquet(batch: &arrow::record_batch::RecordBatch) -> PyResult<String> {
        let tmp = NamedTempFile::new()
            .map_err(|e| PyRuntimeError::new_err(format!("temp file: {e}")))?;
        let path = tmp.path().to_string_lossy().to_string();

        let file = std::fs::File::create(&path)
            .map_err(|e| PyRuntimeError::new_err(format!("open parquet: {e}")))?;
        let mut writer = parquet::arrow::ArrowWriter::try_new(file, batch.schema(), None)
            .map_err(|e| PyRuntimeError::new_err(format!("ArrowWriter: {e}")))?;
        writer.write(batch).map_err(|e| PyRuntimeError::new_err(format!("write: {e}")))?;
        writer.close().map_err(|e| PyRuntimeError::new_err(format!("close: {e}")))?;
        Ok(path)
    }

    let ev_path = write_batch_to_hidden_parquet(&res.eigenvalues)?;
    let ld_path = write_batch_to_hidden_parquet(&res.loadings)?;

    let pyarrow = py.import("pyarrow")?;
    let pq = pyarrow.getattr("parquet")?;

    let ev_table = pq.call_method1("read_table", (ev_path.as_str(),))?;
    let ld_table = pq.call_method1("read_table", (ld_path.as_str(),))?;

    // Cleanup
    let _ = fs::remove_file(&ev_path);
    let _ = fs::remove_file(&ld_path);

    // Build Python dict return value (very ergonomic for the high-level layer)
    let dict = pyo3::types::PyDict::new(py);
    dict.set_item("eigenvalues", ev_table)?;
    dict.set_item("loadings", ld_table)?;
    dict.set_item("n_markers", res.n_markers)?;
    dict.set_item("n_individuals", res.n_individuals)?;
    dict.set_item("n_components", res.n_components)?;
    dict.set_item("total_variance", res.total_variance)?;

    Ok(dict.into())
}
