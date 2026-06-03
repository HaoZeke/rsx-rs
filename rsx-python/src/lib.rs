// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Python bindings for rsx via PyO3.
//!
//! Two surfaces are exposed:
//!
//! 1. **Path-based wrappers** (`process`, `freq`, `depth`, `distrib`,
//!    `signif`, `triage`, `merge`, `pca`) mirror the rsx CLI 1:1 and
//!    drive the existing streaming TSV-based commands.
//!
//! 2. **Arrow entry points** used by the high-level `MarkerTable` API:
//!    - `triage_to_arrow` / `pca_to_arrow` take a TSV path and return
//!      Arrow directly via `commands::*::run_to_arrow`.
//!    - `*_from_arrow` (markers/popmap as Arrow IPC bytes) decode the
//!      bytes into a [`rsx_core::source::MarkerTableSource`] that
//!      either keeps the data in memory or spills it to a Parquet temp
//!      file (driven by a Beissinger-style working-set estimator). The
//!      analysis commands then run against the source through the
//!      generic [`rsx_core::source::MarkerStream`] trait, so no markers
//!      TSV is ever materialised on the Python side.
//!
//! ## Output handling and backend agnosticism
//!
//! When the high-level API (or direct `*_from_arrow` helpers) need to
//! materialize results from commands that currently write TSV outputs
//! (e.g. `freq`, `depth`, `distrib`, `signif`), the Rust side uses pure
//! `pyarrow.csv` to skip an optional leading "#Number of markers" comment and
//! produce a native `pyarrow.Table`. The Python layer then
//! immediately wraps it via `to_narwhals(table)` (see `_adapters.py` and
//! `_read_core_tsv` in `api/markers.py`).
//!
//! The exposed object is always a **narwhals DataFrame** — fully backend
//! agnostic. The concrete implementation is pyarrow by default (lightweight,
//! no forced pandas/polars dependency for internal I/O), but users can
//! convert on demand with `.to_pandas()`, `.to_polars()`, `to_dataframe(backend=...)`,
//! or by passing the result to any narwhals-compatible library (siuba,
//! plotnine, etc.). There is no "pandas fallback" as an implementation detail;
//! pandas is only ever a user-requested output backend.

use arrow::record_batch::RecordBatch;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;

// Custom exception for better error handling from Python side.
// Registered in the module below so Python users see pyrsx.PyrsxError.
pyo3::create_exception!(pyrsx, PyrsxError, PyRuntimeError);
use rsx_core::popmap::Popmap;
use rsx_core::source::MarkerTableSource;
use tempfile::NamedTempFile;

// --------------------------------------------------------------------------- //
// Path-based low-level wrappers (unchanged surface)
// --------------------------------------------------------------------------- //

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
    Python::with_gil(|py| {
        py.allow_threads(|| {
            rsx_core::commands::process::run(&rsx_core::commands::process::ProcessParams {
                input_dir_path: input_dir.to_string(),
                output_file_path: output_file.to_string(),
                n_threads: threads,
                min_depth,
                kmer_dedup,
            })
            .map_err(|e| PyrsxError::new_err(e.to_string()))
        })
    })
}

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
        .map_err(PyrsxError::new_err)?;
    let tm = rsx_core::test_method::TestMethod::parse_str(test).map_err(PyrsxError::new_err)?;

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
    .map_err(|e| PyrsxError::new_err(e.to_string()))
}

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
        .map_err(PyrsxError::new_err)?;
    let tm = rsx_core::test_method::TestMethod::parse_str(test).map_err(PyrsxError::new_err)?;

    Python::with_gil(|py| {
        py.allow_threads(|| {
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
            .map_err(|e| PyrsxError::new_err(e.to_string()))
        })
    })
}

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
    .map_err(|e| PyrsxError::new_err(e.to_string()))
}

#[pyfunction]
#[pyo3(signature = (table_path, output_file, min_depth=1))]
fn freq(table_path: &str, output_file: &str, min_depth: u16) -> PyResult<()> {
    rsx_core::commands::freq::run(&rsx_core::commands::freq::FreqParams {
        markers_table_path: table_path.to_string(),
        output_file_path: output_file.to_string(),
        min_depth,
    })
    .map_err(|e| PyrsxError::new_err(e.to_string()))
}

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
    .map_err(|e| PyrsxError::new_err(e.to_string()))
}

#[pyfunction]
#[pyo3(signature = (input_files, output_file, buffer_size=2000000, output_parquet=false))]
fn merge(
    input_files: Vec<String>,
    output_file: &str,
    buffer_size: usize,
    output_parquet: bool,
) -> PyResult<()> {
    Python::with_gil(|py| {
        py.allow_threads(|| {
            rsx_core::commands::merge::run(&rsx_core::commands::merge::MergeParams {
                input_files,
                output_file_path: output_file.to_string(),
                buffer_size: Some(buffer_size),
                output_parquet,
            })
            .map_err(|e| PyrsxError::new_err(e.to_string()))
        })
    })
}

#[pyfunction]
#[pyo3(signature = (table_path, output_dir, min_depth=1, n_components=None))]
fn pca(
    table_path: &str,
    output_dir: &str,
    min_depth: u16,
    n_components: Option<usize>,
) -> PyResult<()> {
    Python::with_gil(|py| {
        py.allow_threads(|| {
            rsx_core::commands::pca::run(&rsx_core::commands::pca::PcaParams {
                markers_table_path: table_path.to_string(),
                output_dir: output_dir.to_string(),
                min_depth,
                n_components,
            })
            .map_err(|e| PyrsxError::new_err(e.to_string()))
        })
    })
}

// --------------------------------------------------------------------------- //
// Arrow IPC bridge helpers
// --------------------------------------------------------------------------- //

/// Decode an Arrow IPC popmap payload into a real `Popmap`. The popmap is
/// tiny (one row per individual), so we still flow it through a hidden TSV
/// because `Popmap::from_file` is the only constructor exposed by the core.
fn popmap_from_ipc(ipc_bytes: &[u8]) -> PyResult<(Popmap, NamedTempFile)> {
    use std::io::Write;

    if ipc_bytes.is_empty() {
        return Err(PyrsxError::new_err(
            "popmap Arrow payload is empty — pass a non-empty popmap DataFrame",
        ));
    }

    let cursor = std::io::Cursor::new(ipc_bytes);
    let reader = arrow::ipc::reader::StreamReader::try_new(cursor, None)
        .map_err(|e| PyrsxError::new_err(format!("popmap IPC reader: {e}")))?;
    let schema = reader.schema();
    if schema.fields().len() < 2 {
        return Err(PyrsxError::new_err(format!(
            "popmap Arrow table needs at least 2 columns (individual, group); got {}",
            schema.fields().len()
        )));
    }

    let tmp =
        NamedTempFile::new().map_err(|e| PyrsxError::new_err(format!("popmap tempfile: {e}")))?;
    let file = std::fs::File::create(tmp.path())
        .map_err(|e| PyrsxError::new_err(format!("create popmap tsv: {e}")))?;
    let mut w = std::io::BufWriter::new(file);

    for batch in reader {
        let batch = batch.map_err(|e| PyrsxError::new_err(format!("popmap batch: {e}")))?;
        let ind = batch.column(0);
        let grp = batch.column(1);
        for row in 0..batch.num_rows() {
            let i = scalar_to_string(ind.as_ref(), row);
            let g = scalar_to_string(grp.as_ref(), row);
            writeln!(w, "{i}\t{g}")
                .map_err(|e| PyrsxError::new_err(format!("write popmap row: {e}")))?;
        }
    }
    w.flush()
        .map_err(|e| PyrsxError::new_err(format!("flush popmap tsv: {e}")))?;
    drop(w);

    let popmap = Popmap::from_file(tmp.path())
        .map_err(|e| PyrsxError::new_err(format!("read popmap: {e}")))?;
    Ok((popmap, tmp))
}

fn scalar_to_string(array: &dyn arrow::array::Array, row: usize) -> String {
    use arrow::array::*;
    use arrow::datatypes::DataType;
    if array.is_null(row) {
        return String::new();
    }
    match array.data_type() {
        DataType::Utf8 => array
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap()
            .value(row)
            .to_string(),
        DataType::LargeUtf8 => array
            .as_any()
            .downcast_ref::<LargeStringArray>()
            .unwrap()
            .value(row)
            .to_string(),
        DataType::Int32 => array
            .as_any()
            .downcast_ref::<Int32Array>()
            .unwrap()
            .value(row)
            .to_string(),
        DataType::Int64 => array
            .as_any()
            .downcast_ref::<Int64Array>()
            .unwrap()
            .value(row)
            .to_string(),
        DataType::UInt32 => array
            .as_any()
            .downcast_ref::<UInt32Array>()
            .unwrap()
            .value(row)
            .to_string(),
        DataType::UInt64 => array
            .as_any()
            .downcast_ref::<UInt64Array>()
            .unwrap()
            .value(row)
            .to_string(),
        _ => String::new(),
    }
}

/// Read an output TSV produced by a low-level rsx command and return it
/// to Python as a native `pyarrow.Table`.
///
/// Core command outputs may begin with a single "#Number of markers : N"
/// metadata comment line followed by the real TSV header. We skip that
/// comment when it is present for a clean parse without any pandas dependency
/// inside the extension.
///
/// This is deliberately kept as a pure pyarrow Table (the Arrow interchange
/// format between Rust and Python). The *caller* (in the high-level Python
/// API) is expected to wrap it with `to_narwhals(table)`. The resulting
/// narwhals DataFrame is fully **backend agnostic** — its concrete backend
/// is pyarrow by default, but users can request pandas/polars/etc. on demand
/// via the normal narwhals conversion paths (`.to_pandas()`, `to_dataframe(backend=...)`,
/// etc.). This ensures there is never a "pandas fallback" or hard-coded
/// backend as an implementation detail when materializing command outputs.
///
/// See `_read_core_tsv` in `python/pyrsx/api/markers.py` for the matching
/// pure-Python helper used by the path-backed high-level `MarkerTable`
/// methods, and `_adapters.py` for the to/from narwhals conversions.
fn read_tsv_to_pyarrow_table(py: Python<'_>, path: &str) -> PyResult<PyObject> {
    use std::io::BufRead;

    let file = std::fs::File::open(path)
        .map_err(|e| PyrsxError::new_err(format!("read TSV header: {e}")))?;
    let mut reader = std::io::BufReader::new(file);
    let mut first_line = String::new();
    reader
        .read_line(&mut first_line)
        .map_err(|e| PyrsxError::new_err(format!("read TSV header: {e}")))?;
    let skip_rows = if first_line.starts_with('#') { 1 } else { 0 };

    let pa_csv = py.import("pyarrow.csv")?;

    let read_opts_dict = pyo3::types::PyDict::new(py);
    read_opts_dict.set_item("skip_rows", skip_rows)?;
    read_opts_dict.set_item("use_threads", true)?;
    let read_opts = pa_csv
        .getattr("ReadOptions")?
        .call((), Some(&read_opts_dict))?;

    let parse_opts_dict = pyo3::types::PyDict::new(py);
    parse_opts_dict.set_item("delimiter", "\t")?;
    let parse_opts = pa_csv
        .getattr("ParseOptions")?
        .call((), Some(&parse_opts_dict))?;

    let read_dict = pyo3::types::PyDict::new(py);
    read_dict.set_item("read_options", read_opts)?;
    read_dict.set_item("parse_options", parse_opts)?;

    let table = pa_csv
        .getattr("read_csv")?
        .call((path,), Some(&read_dict))?;

    Ok(table.into())
}

/// Build a `pyarrow.Table` from RecordBatches via Arrow IPC. Pure RAM,
/// never touches disk.
fn batches_to_pyarrow_table(py: Python<'_>, batches: &[RecordBatch]) -> PyResult<PyObject> {
    if batches.is_empty() {
        let pyarrow = py.import("pyarrow")?;
        return Ok(pyarrow
            .getattr("Table")?
            .call_method0("from_batches")?
            .into());
    }
    let bytes = batches_to_ipc_bytes(&batches.iter().collect::<Vec<_>>())?;
    let mut tables = ipc_bytes_to_pyarrow_tables(py, &bytes)?;
    if tables.is_empty() {
        let pyarrow = py.import("pyarrow")?;
        return Ok(pyarrow
            .getattr("Table")?
            .call_method0("from_batches")?
            .into());
    }
    Ok(tables.remove(0))
}

fn batches_to_ipc_bytes(batches: &[&RecordBatch]) -> PyResult<Vec<u8>> {
    if batches.is_empty() {
        return Ok(Vec::new());
    }
    let mut buf = Vec::new();
    {
        let mut writer = arrow::ipc::writer::StreamWriter::try_new(&mut buf, &batches[0].schema())
            .map_err(|e| PyrsxError::new_err(format!("IPC writer: {e}")))?;
        for b in batches {
            writer
                .write(b)
                .map_err(|e| PyrsxError::new_err(format!("IPC write: {e}")))?;
        }
        writer
            .finish()
            .map_err(|e| PyrsxError::new_err(format!("IPC finish: {e}")))?;
    }
    Ok(buf)
}

fn ipc_bytes_to_pyarrow_tables(py: Python<'_>, bytes: &[u8]) -> PyResult<Vec<PyObject>> {
    if bytes.is_empty() {
        return Ok(vec![]);
    }
    let py_bytes = pyo3::types::PyBytes::new(py, bytes);
    let pyarrow = py.import("pyarrow")?;
    let ipc = pyarrow.getattr("ipc")?;
    let reader = ipc.call_method1("RecordBatchStreamReader", (py_bytes,))?;
    let mut tables = Vec::new();
    while let Ok(batch) = reader.call_method0("read_next_batch") {
        let t = pyarrow
            .getattr("Table")?
            .call_method1("from_batches", (vec![batch],))?;
        tables.push(t.into());
    }
    Ok(tables)
}

// --------------------------------------------------------------------------- //
// Arrow entry points keyed on a TSV path (unchanged surface).
// --------------------------------------------------------------------------- //

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
        .map_err(|e| PyrsxError::new_err(e.to_string()))?;
    batches_to_pyarrow_table(py, &batches)
}

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
        .map_err(|e| PyrsxError::new_err(e.to_string()))?;

    let eigenvalues = batches_to_pyarrow_table(py, &[res.eigenvalues])?;
    let loadings = batches_to_pyarrow_table(py, &[res.loadings])?;

    let dict = pyo3::types::PyDict::new(py);
    dict.set_item("eigenvalues", eigenvalues)?;
    dict.set_item("loadings", loadings)?;
    dict.set_item("n_markers", res.n_markers)?;
    dict.set_item("n_individuals", res.n_individuals)?;
    dict.set_item("n_components", res.n_components)?;
    dict.set_item("total_variance", res.total_variance)?;
    Ok(dict.into())
}

// --------------------------------------------------------------------------- //
// from-Arrow entry points (Arrow IPC bytes in → Arrow Table / pyarrow out).
// All of these route the markers payload through MarkerTableSource so the
// command runs against either in-memory RecordBatches or a Parquet spill,
// never through a hidden markers TSV.
// --------------------------------------------------------------------------- //

/// Command-specific working-set multiplier for the spill heuristic.
mod cmd_overhead {
    /// freq is a single streaming pass with no per-marker accumulators.
    pub const FREQ: f64 = 1.3;
    /// depth keeps per-individual buffers; default mode is in-memory.
    pub const DEPTH: f64 = 1.8;
    /// distrib accumulates a 2D group table.
    pub const DISTRIB: f64 = 1.5;
    /// signif may keep full p-value vectors for FDR; widen the prediction.
    pub const SIGNIF: f64 = 2.0;
    /// triage builds Arrow output rows in addition to streaming work.
    pub const TRIAGE: f64 = 2.0;
    /// pca holds the full Gram matrix + intermediate Mean vector.
    pub const PCA: f64 = 1.5;
}

#[pyfunction]
#[pyo3(signature = (markers_ipc, popmap_ipc, min_depth=1, posterior_threshold=0.9, prior_probability=0.01, linked_probability=0.9, group1="", group2=""))]
#[allow(clippy::too_many_arguments)]
fn triage_to_arrow_from_arrow(
    py: Python<'_>,
    markers_ipc: &[u8],
    popmap_ipc: &[u8],
    min_depth: u16,
    posterior_threshold: f64,
    prior_probability: f64,
    linked_probability: f64,
    group1: &str,
    group2: &str,
) -> PyResult<PyObject> {
    let (popmap, _popmap_tmp) = popmap_from_ipc(popmap_ipc)?;
    let source = MarkerTableSource::from_arrow_ipc(
        markers_ipc,
        Some(&popmap),
        min_depth,
        cmd_overhead::TRIAGE,
    )
    .map_err(|e| PyrsxError::new_err(format!("MarkerTableSource: {e}")))?;

    let params = rsx_core::commands::triage::TriageParams {
        markers_table_path: String::new(),
        popmap_file_path: String::new(),
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

    let batches = rsx_core::commands::triage::run_to_arrow_with_source(&source, &popmap, &params)
        .map_err(|e| PyrsxError::new_err(e.to_string()))?;
    batches_to_pyarrow_table(py, &batches)
}

#[pyfunction]
#[pyo3(signature = (markers_ipc, min_depth=1, n_components=None))]
fn pca_to_arrow_from_arrow(
    py: Python<'_>,
    markers_ipc: &[u8],
    min_depth: u16,
    n_components: Option<usize>,
) -> PyResult<PyObject> {
    let source = MarkerTableSource::from_arrow_ipc(markers_ipc, None, min_depth, cmd_overhead::PCA)
        .map_err(|e| PyrsxError::new_err(format!("MarkerTableSource: {e}")))?;

    let params = rsx_core::commands::pca::PcaParams {
        markers_table_path: String::new(),
        output_dir: String::new(),
        min_depth,
        n_components,
    };
    let res = rsx_core::commands::pca::run_to_arrow_with_source(&source, &params)
        .map_err(|e| PyrsxError::new_err(e.to_string()))?;

    let eigenvalues = batches_to_pyarrow_table(py, &[res.eigenvalues])?;
    let loadings = batches_to_pyarrow_table(py, &[res.loadings])?;

    let dict = pyo3::types::PyDict::new(py);
    dict.set_item("eigenvalues", eigenvalues)?;
    dict.set_item("loadings", loadings)?;
    dict.set_item("n_markers", res.n_markers)?;
    dict.set_item("n_individuals", res.n_individuals)?;
    dict.set_item("n_components", res.n_components)?;
    dict.set_item("total_variance", res.total_variance)?;
    Ok(dict.into())
}

/// Low-level: run freq on an Arrow IPC markers payload and return a
/// pyarrow.Table.
///
/// The table is produced by writing a temp TSV from the core (via
/// run_with_source) and reading it back with the pure-pyarrow
/// `read_tsv_to_pyarrow_table`. The high-level API (and direct callers)
/// are expected to immediately do `to_narwhals(...)` on the result so the
/// user sees a backend-agnostic narwhals DataFrame (see Python-side
/// `_read_core_tsv` and adapters for the symmetric path used by
/// path-backed `MarkerTable`).
#[pyfunction]
#[pyo3(signature = (markers_ipc, min_depth=1))]
fn freq_from_arrow(py: Python<'_>, markers_ipc: &[u8], min_depth: u16) -> PyResult<PyObject> {
    let source =
        MarkerTableSource::from_arrow_ipc(markers_ipc, None, min_depth, cmd_overhead::FREQ)
            .map_err(|e| PyrsxError::new_err(format!("MarkerTableSource: {e}")))?;

    let out =
        NamedTempFile::new().map_err(|e| PyrsxError::new_err(format!("freq output temp: {e}")))?;
    let out_path = out.path().to_string_lossy().to_string();

    let params = rsx_core::commands::freq::FreqParams {
        markers_table_path: String::new(),
        output_file_path: out_path.clone(),
        min_depth,
    };

    py.allow_threads(|| {
        rsx_core::commands::freq::run_with_source(&source, &params)
            .map_err(|e| PyrsxError::new_err(e.to_string()))
    })?;

    read_tsv_to_pyarrow_table(py, &out_path)
}

/// Low-level: run depth on Arrow IPC payloads (markers + popmap) and
/// return a pyarrow.Table.
///
/// See `freq_from_arrow` for the general pattern and narwhals
/// backend-agnostic consumption story.
#[pyfunction]
#[pyo3(signature = (markers_ipc, popmap_ipc, min_frequency=0.75))]
fn depth_from_arrow(
    py: Python<'_>,
    markers_ipc: &[u8],
    popmap_ipc: &[u8],
    min_frequency: f32,
) -> PyResult<PyObject> {
    let (popmap, _popmap_tmp) = popmap_from_ipc(popmap_ipc)?;
    // depth reads `individual_depths`; min_depth=1 matches the file-based
    // behaviour where the parser stores raw depths and the command does
    // its own thresholding via `min_frequency`.
    let source =
        MarkerTableSource::from_arrow_ipc(markers_ipc, Some(&popmap), 1, cmd_overhead::DEPTH)
            .map_err(|e| PyrsxError::new_err(format!("MarkerTableSource: {e}")))?;

    let out =
        NamedTempFile::new().map_err(|e| PyrsxError::new_err(format!("depth output temp: {e}")))?;
    let out_path = out.path().to_string_lossy().to_string();

    let params = rsx_core::commands::depth::DepthParams {
        markers_table_path: String::new(),
        popmap_file_path: String::new(),
        output_file_path: out_path.clone(),
        min_frequency,
        // External-sort streaming mode would write LZ4 sort chunks to a temp
        // dir; the in-memory path is fine since the source already fit in
        // RAM (or is paged via Parquet, which is roughly the same cost).
        streaming: false,
    };
    py.allow_threads(|| {
        rsx_core::commands::depth::run_with_source(&source, &popmap, &params)
            .map_err(|e| PyrsxError::new_err(e.to_string()))
    })?;

    read_tsv_to_pyarrow_table(py, &out_path)
}

#[pyfunction]
#[pyo3(signature = (markers_ipc, popmap_ipc, min_depth=1, signif_threshold=0.05, group1="", group2="", correction="bonferroni", test="chisq"))]
#[allow(clippy::too_many_arguments)]
fn distrib_from_arrow(
    py: Python<'_>,
    markers_ipc: &[u8],
    popmap_ipc: &[u8],
    min_depth: u16,
    signif_threshold: f32,
    group1: &str,
    group2: &str,
    correction: &str,
    test: &str,
) -> PyResult<PyObject> {
    let (popmap, _popmap_tmp) = popmap_from_ipc(popmap_ipc)?;
    let source = MarkerTableSource::from_arrow_ipc(
        markers_ipc,
        Some(&popmap),
        min_depth,
        cmd_overhead::DISTRIB,
    )
    .map_err(|e| PyrsxError::new_err(format!("MarkerTableSource: {e}")))?;

    let out = NamedTempFile::new()
        .map_err(|e| PyrsxError::new_err(format!("distrib output temp: {e}")))?;
    let out_path = out.path().to_string_lossy().to_string();

    let corr = rsx_core::test_method::CorrectionMethod::parse_str(correction)
        .map_err(PyrsxError::new_err)?;
    let tm = rsx_core::test_method::TestMethod::parse_str(test).map_err(PyrsxError::new_err)?;

    let params = rsx_core::commands::distrib::DistribParams {
        markers_table_path: String::new(),
        popmap_file_path: String::new(),
        output_file_path: out_path.clone(),
        min_depth,
        signif_threshold,
        correction: corr,
        test_method: tm,
        output_bayes: false,
        group1: group1.to_string(),
        group2: group2.to_string(),
    };
    py.allow_threads(|| {
        rsx_core::commands::distrib::run_with_source(&source, &popmap, &params)
            .map_err(|e| PyrsxError::new_err(e.to_string()))
    })?;

    read_tsv_to_pyarrow_table(py, &out_path)
}

#[pyfunction]
#[pyo3(signature = (markers_ipc, popmap_ipc, min_depth=1, signif_threshold=0.05, group1="", group2="", correction="bonferroni", test="chisq", output_fasta=false, bayes=false))]
#[allow(clippy::too_many_arguments)]
fn signif_from_arrow(
    py: Python<'_>,
    markers_ipc: &[u8],
    popmap_ipc: &[u8],
    min_depth: u16,
    signif_threshold: f32,
    group1: &str,
    group2: &str,
    correction: &str,
    test: &str,
    output_fasta: bool,
    bayes: bool,
) -> PyResult<PyObject> {
    let (popmap, _popmap_tmp) = popmap_from_ipc(popmap_ipc)?;
    let source = MarkerTableSource::from_arrow_ipc(
        markers_ipc,
        Some(&popmap),
        min_depth,
        cmd_overhead::SIGNIF,
    )
    .map_err(|e| PyrsxError::new_err(format!("MarkerTableSource: {e}")))?;

    let out = NamedTempFile::new()
        .map_err(|e| PyrsxError::new_err(format!("signif output temp: {e}")))?;
    let out_path = out.path().to_string_lossy().to_string();

    let corr = rsx_core::test_method::CorrectionMethod::parse_str(correction)
        .map_err(PyrsxError::new_err)?;
    let tm = rsx_core::test_method::TestMethod::parse_str(test).map_err(PyrsxError::new_err)?;

    let params = rsx_core::commands::signif::SignifParams {
        markers_table_path: String::new(),
        popmap_file_path: String::new(),
        output_file_path: out_path.clone(),
        min_depth,
        signif_threshold,
        correction: corr,
        test_method: tm,
        output_fasta,
        output_bayes: bayes,
        group1: group1.to_string(),
        group2: group2.to_string(),
    };
    py.allow_threads(|| {
        rsx_core::commands::signif::run_with_source(&source, &popmap, &params)
            .map_err(|e| PyrsxError::new_err(e.to_string()))
    })?;

    read_tsv_to_pyarrow_table(py, &out_path)
}

// --------------------------------------------------------------------------- //
// PyO3 module registration
// --------------------------------------------------------------------------- //

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

    m.add_function(wrap_pyfunction!(triage_to_arrow, m)?)?;
    m.add_function(wrap_pyfunction!(pca_to_arrow, m)?)?;
    m.add_function(wrap_pyfunction!(triage_to_arrow_from_arrow, m)?)?;
    m.add_function(wrap_pyfunction!(pca_to_arrow_from_arrow, m)?)?;
    m.add_function(wrap_pyfunction!(freq_from_arrow, m)?)?;
    m.add_function(wrap_pyfunction!(depth_from_arrow, m)?)?;
    m.add_function(wrap_pyfunction!(distrib_from_arrow, m)?)?;
    m.add_function(wrap_pyfunction!(signif_from_arrow, m)?)?;

    // Custom exception for idiomatic error handling in Python
    m.add("PyrsxError", m.py().get_type::<PyrsxError>())?;

    Ok(())
}
