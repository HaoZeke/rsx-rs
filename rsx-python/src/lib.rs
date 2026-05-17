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
    m.add_function(wrap_pyfunction!(freq, m)?)?;
    m.add_function(wrap_pyfunction!(depth, m)?)?;
    m.add_function(wrap_pyfunction!(merge, m)?)?;
    m.add_function(wrap_pyfunction!(pca, m)?)?;
    Ok(())
}
