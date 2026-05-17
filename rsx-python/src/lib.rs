// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Python bindings for rsx via PyO3.
//!
//! The bindings come in two flavours:
//!
//! 1. Path-based low-level wrappers (`process`, `freq`, `depth`, `distrib`,
//!    `signif`, `triage`, `merge`, `pca`) that mirror the rsx CLI 1:1 and
//!    drive the existing streaming TSV-based commands.
//!
//! 2. Arrow entry points used by the high-level `MarkerTable` API:
//!    - `triage_to_arrow` / `pca_to_arrow` take a TSV path and return Arrow
//!      directly via `commands::*::run_to_arrow`.
//!    - `*_from_arrow` (markers/popmap as Arrow IPC bytes) decode the bytes
//!      back to a real markers/popmap TSV in a hidden temp file, then call
//!      the standard command. This keeps Python from ever seeing a temp file
//!      and lets us reuse the well-tested streaming code paths instead of
//!      re-implementing every command on top of in-memory RecordBatches.

use arrow::array::Array;
use arrow::record_batch::RecordBatch;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use std::io::Write;
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

// --------------------------------------------------------------------------- //
// Arrow ↔ TSV bridge helpers (Rust-only, no Python visible artefacts)
// --------------------------------------------------------------------------- //

/// Format any supported Arrow scalar as a TSV-safe string.
///
/// Null becomes the empty string. Floats are rounded to the nearest integer
/// because the markers TSV parser (`fast_parse_u16`) only understands integer
/// ASCII for depth columns — DataFrames coming in as float (e.g. pandas with
/// NaNs) still round-trip correctly because integer depths survive an f64
/// cast unchanged.
fn array_value_as_string(array: &dyn Array, row: usize) -> PyResult<String> {
    use arrow::array::*;
    use arrow::datatypes::DataType;

    if array.is_null(row) {
        return Ok(String::new());
    }
    Ok(match array.data_type() {
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
        DataType::Int8 => array.as_any().downcast_ref::<Int8Array>().unwrap().value(row).to_string(),
        DataType::Int16 => array.as_any().downcast_ref::<Int16Array>().unwrap().value(row).to_string(),
        DataType::Int32 => array.as_any().downcast_ref::<Int32Array>().unwrap().value(row).to_string(),
        DataType::Int64 => array.as_any().downcast_ref::<Int64Array>().unwrap().value(row).to_string(),
        DataType::UInt8 => array.as_any().downcast_ref::<UInt8Array>().unwrap().value(row).to_string(),
        DataType::UInt16 => array.as_any().downcast_ref::<UInt16Array>().unwrap().value(row).to_string(),
        DataType::UInt32 => array.as_any().downcast_ref::<UInt32Array>().unwrap().value(row).to_string(),
        DataType::UInt64 => array.as_any().downcast_ref::<UInt64Array>().unwrap().value(row).to_string(),
        DataType::Float32 => {
            let v = array.as_any().downcast_ref::<Float32Array>().unwrap().value(row);
            (v.round() as i64).to_string()
        }
        DataType::Float64 => {
            let v = array.as_any().downcast_ref::<Float64Array>().unwrap().value(row);
            (v.round() as i64).to_string()
        }
        other => {
            return Err(PyRuntimeError::new_err(format!(
                "Arrow column type {other:?} not supported in markers/popmap conversion"
            )));
        }
    })
}

/// Decode IPC bytes into the full set of RecordBatches.
fn ipc_bytes_to_batches(ipc_bytes: &[u8]) -> PyResult<(arrow::datatypes::SchemaRef, Vec<RecordBatch>)> {
    if ipc_bytes.is_empty() {
        return Err(PyRuntimeError::new_err(
            "received empty IPC payload — pass a non-empty markers/popmap DataFrame",
        ));
    }
    let cursor = std::io::Cursor::new(ipc_bytes);
    let reader = arrow::ipc::reader::StreamReader::try_new(cursor, None)
        .map_err(|e| PyRuntimeError::new_err(format!("IPC reader: {e}")))?;
    let schema = reader.schema();
    let mut batches: Vec<RecordBatch> = Vec::new();
    for b in reader {
        let b = b.map_err(|e| PyRuntimeError::new_err(format!("IPC batch: {e}")))?;
        batches.push(b);
    }
    Ok((schema, batches))
}

/// Materialise an Arrow IPC markers blob into the exact TSV shape that
/// `MarkersTableStream::open` understands:
///   `#Number of markers : N`
///   `id<TAB>sequence<TAB>ind1<TAB>...<TAB>indK`
///   `<id><TAB><sequence><TAB><depth1>...`
///
/// The returned `NamedTempFile` keeps the file alive until the caller drops it.
fn ipc_to_markers_tsv(ipc_bytes: &[u8]) -> PyResult<NamedTempFile> {
    let (schema, batches) = ipc_bytes_to_batches(ipc_bytes)?;
    let n_cols = schema.fields().len();
    if n_cols < 3 {
        return Err(PyRuntimeError::new_err(format!(
            "Markers Arrow table needs id, sequence, and >=1 individual column; got {n_cols}"
        )));
    }
    let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();

    let tmp = NamedTempFile::new()
        .map_err(|e| PyRuntimeError::new_err(format!("markers tempfile: {e}")))?;
    let file = std::fs::File::create(tmp.path())
        .map_err(|e| PyRuntimeError::new_err(format!("create markers tsv: {e}")))?;
    let mut w = std::io::BufWriter::new(file);

    writeln!(w, "#Number of markers : {total_rows}")
        .map_err(|e| PyRuntimeError::new_err(format!("write header: {e}")))?;

    let col_names: Vec<&str> = schema.fields().iter().map(|f| f.name().as_str()).collect();
    writeln!(w, "{}", col_names.join("\t"))
        .map_err(|e| PyRuntimeError::new_err(format!("write column line: {e}")))?;

    let mut row_buf = String::with_capacity(64 + n_cols * 4);
    for batch in &batches {
        let cols: Vec<&dyn Array> = (0..batch.num_columns())
            .map(|i| batch.column(i).as_ref())
            .collect();
        for row in 0..batch.num_rows() {
            row_buf.clear();
            for (idx, col) in cols.iter().enumerate() {
                if idx > 0 {
                    row_buf.push('\t');
                }
                row_buf.push_str(&array_value_as_string(*col, row)?);
            }
            row_buf.push('\n');
            w.write_all(row_buf.as_bytes())
                .map_err(|e| PyRuntimeError::new_err(format!("write row: {e}")))?;
        }
    }
    w.flush().map_err(|e| PyRuntimeError::new_err(format!("flush markers tsv: {e}")))?;
    Ok(tmp)
}

/// Materialise an Arrow IPC popmap blob into `individual<TAB>group` per line.
/// Uses the first two columns of the schema; any extras are ignored.
fn ipc_to_popmap_tsv(ipc_bytes: &[u8]) -> PyResult<NamedTempFile> {
    let (schema, batches) = ipc_bytes_to_batches(ipc_bytes)?;
    if schema.fields().len() < 2 {
        return Err(PyRuntimeError::new_err(format!(
            "Popmap Arrow table needs at least 2 columns (individual, group); got {}",
            schema.fields().len()
        )));
    }
    let tmp = NamedTempFile::new()
        .map_err(|e| PyRuntimeError::new_err(format!("popmap tempfile: {e}")))?;
    let file = std::fs::File::create(tmp.path())
        .map_err(|e| PyRuntimeError::new_err(format!("create popmap tsv: {e}")))?;
    let mut w = std::io::BufWriter::new(file);

    for batch in &batches {
        let ind = batch.column(0);
        let grp = batch.column(1);
        for row in 0..batch.num_rows() {
            let i = array_value_as_string(ind.as_ref(), row)?;
            let g = array_value_as_string(grp.as_ref(), row)?;
            writeln!(w, "{i}\t{g}")
                .map_err(|e| PyRuntimeError::new_err(format!("write popmap row: {e}")))?;
        }
    }
    w.flush().map_err(|e| PyRuntimeError::new_err(format!("flush popmap tsv: {e}")))?;
    Ok(tmp)
}

/// Read an output TSV produced by a low-level command and hand it back to
/// Python as a `pyarrow.Table`. Comments (`#…`) are stripped; tab-delimited;
/// header row honoured.
fn read_tsv_to_pyarrow_table(py: Python<'_>, path: &str) -> PyResult<PyObject> {
    let pyarrow = py.import("pyarrow")?;
    let csv = pyarrow.getattr("csv")?;
    let read_options = csv.getattr("ReadOptions")?;
    let parse_options = csv.getattr("ParseOptions")?;
    let convert_options = csv.getattr("ConvertOptions")?;

    let ro_kw = pyo3::types::PyDict::new(py);
    ro_kw.set_item("use_threads", false)?;
    let read_opts = read_options.call((), Some(&ro_kw))?;

    let po_kw = pyo3::types::PyDict::new(py);
    po_kw.set_item("delimiter", "\t")?;
    let parse_opts = parse_options.call((), Some(&po_kw))?;

    let co_kw = pyo3::types::PyDict::new(py);
    let convert_opts = convert_options.call((), Some(&co_kw))?;

    let kwargs = pyo3::types::PyDict::new(py);
    kwargs.set_item("read_options", read_opts)?;
    kwargs.set_item("parse_options", parse_opts)?;
    kwargs.set_item("convert_options", convert_opts)?;

    let table = csv.call_method("read_csv", (path,), Some(&kwargs))?;
    Ok(table.into())
}

// --------------------------------------------------------------------------- //
// Arrow output paths from a TSV path (unchanged surface, light tidy-up)
// --------------------------------------------------------------------------- //

/// Build a `pyarrow.Table` from RecordBatches via Arrow IPC. Pure RAM,
/// never touches disk.
fn batches_to_pyarrow_table(py: Python<'_>, batches: &[RecordBatch]) -> PyResult<PyObject> {
    if batches.is_empty() {
        let pyarrow = py.import("pyarrow")?;
        return Ok(pyarrow.getattr("Table")?.call_method0("from_batches")?.into());
    }
    let bytes = batches_to_ipc_bytes(&batches.iter().collect::<Vec<_>>())?;
    let mut tables = ipc_bytes_to_pyarrow_tables(py, &bytes)?;
    if tables.is_empty() {
        let pyarrow = py.import("pyarrow")?;
        return Ok(pyarrow.getattr("Table")?.call_method0("from_batches")?.into());
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
            .map_err(|e| PyRuntimeError::new_err(format!("IPC writer: {e}")))?;
        for b in batches {
            writer
                .write(b)
                .map_err(|e| PyRuntimeError::new_err(format!("IPC write: {e}")))?;
        }
        writer
            .finish()
            .map_err(|e| PyRuntimeError::new_err(format!("IPC finish: {e}")))?;
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

/// Returns the triage results directly as a `pyarrow.Table` for a TSV input.
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
    batches_to_pyarrow_table(py, &batches)
}

/// Returns PCA results as a dict of `pyarrow.Table`s for a TSV input.
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

    let bytes = batches_to_ipc_bytes(&[&res.eigenvalues, &res.loadings])?;
    let tables = ipc_bytes_to_pyarrow_tables(py, &bytes)?;

    let dict = pyo3::types::PyDict::new(py);
    if tables.len() >= 2 {
        dict.set_item("eigenvalues", &tables[0])?;
        dict.set_item("loadings", &tables[1])?;
    }
    dict.set_item("n_markers", res.n_markers)?;
    dict.set_item("n_individuals", res.n_individuals)?;
    dict.set_item("n_components", res.n_components)?;
    dict.set_item("total_variance", res.total_variance)?;
    Ok(dict.into())
}

// --------------------------------------------------------------------------- //
// from-Arrow entry points (Arrow IPC bytes in → Arrow Table / pyarrow out)
//
// All of these convert the IPC payload back to the real markers/popmap TSV
// in hidden temp files and reuse the existing low-level commands. This is
// not the eventual MarkerTableSource design (see
// docs/orgmode/reference/arrow-input-for-high-level-api.org), but it is
// correct, exercises the same code paths the CLI exercises, and removes the
// previously-shipped placeholders that fed binary IPC into pandas.read_csv.
// --------------------------------------------------------------------------- //

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
    let markers_tmp = ipc_to_markers_tsv(markers_ipc)?;
    let popmap_tmp = ipc_to_popmap_tsv(popmap_ipc)?;

    let params = rsx_core::commands::triage::TriageParams {
        markers_table_path: markers_tmp.path().to_string_lossy().to_string(),
        popmap_file_path: popmap_tmp.path().to_string_lossy().to_string(),
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
    batches_to_pyarrow_table(py, &batches)
    // markers_tmp and popmap_tmp are dropped here, which removes both files.
}

#[pyfunction]
#[pyo3(signature = (markers_ipc, min_depth=1, n_components=None))]
fn pca_to_arrow_from_arrow(
    py: Python<'_>,
    markers_ipc: &[u8],
    min_depth: u16,
    n_components: Option<usize>,
) -> PyResult<PyObject> {
    let markers_tmp = ipc_to_markers_tsv(markers_ipc)?;
    let params = rsx_core::commands::pca::PcaParams {
        markers_table_path: markers_tmp.path().to_string_lossy().to_string(),
        output_dir: String::new(),
        min_depth,
        n_components,
    };
    let res = rsx_core::commands::pca::run_to_arrow(&params)
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

    let bytes = batches_to_ipc_bytes(&[&res.eigenvalues, &res.loadings])?;
    let tables = ipc_bytes_to_pyarrow_tables(py, &bytes)?;

    let dict = pyo3::types::PyDict::new(py);
    if tables.len() >= 2 {
        dict.set_item("eigenvalues", &tables[0])?;
        dict.set_item("loadings", &tables[1])?;
    }
    dict.set_item("n_markers", res.n_markers)?;
    dict.set_item("n_individuals", res.n_individuals)?;
    dict.set_item("n_components", res.n_components)?;
    dict.set_item("total_variance", res.total_variance)?;
    Ok(dict.into())
}

#[pyfunction]
#[pyo3(signature = (markers_ipc, min_depth=1))]
fn freq_from_arrow(py: Python<'_>, markers_ipc: &[u8], min_depth: u16) -> PyResult<PyObject> {
    let markers_tmp = ipc_to_markers_tsv(markers_ipc)?;
    let out = NamedTempFile::new()
        .map_err(|e| PyRuntimeError::new_err(format!("freq output temp: {e}")))?;
    let out_path = out.path().to_string_lossy().to_string();

    rsx_core::commands::freq::run(&rsx_core::commands::freq::FreqParams {
        markers_table_path: markers_tmp.path().to_string_lossy().to_string(),
        output_file_path: out_path.clone(),
        min_depth,
    })
    .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

    read_tsv_to_pyarrow_table(py, &out_path)
}

#[pyfunction]
#[pyo3(signature = (markers_ipc, popmap_ipc, min_frequency=0.75))]
fn depth_from_arrow(
    py: Python<'_>,
    markers_ipc: &[u8],
    popmap_ipc: &[u8],
    min_frequency: f32,
) -> PyResult<PyObject> {
    let markers_tmp = ipc_to_markers_tsv(markers_ipc)?;
    let popmap_tmp = ipc_to_popmap_tsv(popmap_ipc)?;
    let out = NamedTempFile::new()
        .map_err(|e| PyRuntimeError::new_err(format!("depth output temp: {e}")))?;
    let out_path = out.path().to_string_lossy().to_string();

    let markers_path = markers_tmp.path().to_string_lossy().to_string();
    let file_size = std::fs::metadata(&markers_path).map(|m| m.len()).unwrap_or(0);

    rsx_core::commands::depth::run(&rsx_core::commands::depth::DepthParams {
        markers_table_path: markers_path,
        popmap_file_path: popmap_tmp.path().to_string_lossy().to_string(),
        output_file_path: out_path.clone(),
        min_frequency,
        streaming: file_size > 2_000_000_000,
    })
    .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

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
    let markers_tmp = ipc_to_markers_tsv(markers_ipc)?;
    let popmap_tmp = ipc_to_popmap_tsv(popmap_ipc)?;
    let out = NamedTempFile::new()
        .map_err(|e| PyRuntimeError::new_err(format!("distrib output temp: {e}")))?;
    let out_path = out.path().to_string_lossy().to_string();

    let corr = rsx_core::test_method::CorrectionMethod::parse_str(correction)
        .map_err(PyRuntimeError::new_err)?;
    let tm = rsx_core::test_method::TestMethod::parse_str(test).map_err(PyRuntimeError::new_err)?;

    rsx_core::commands::distrib::run(&rsx_core::commands::distrib::DistribParams {
        markers_table_path: markers_tmp.path().to_string_lossy().to_string(),
        popmap_file_path: popmap_tmp.path().to_string_lossy().to_string(),
        output_file_path: out_path.clone(),
        min_depth,
        signif_threshold,
        correction: corr,
        test_method: tm,
        output_bayes: false,
        group1: group1.to_string(),
        group2: group2.to_string(),
    })
    .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

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
    let markers_tmp = ipc_to_markers_tsv(markers_ipc)?;
    let popmap_tmp = ipc_to_popmap_tsv(popmap_ipc)?;
    let out = NamedTempFile::new()
        .map_err(|e| PyRuntimeError::new_err(format!("signif output temp: {e}")))?;
    let out_path = out.path().to_string_lossy().to_string();

    let corr = rsx_core::test_method::CorrectionMethod::parse_str(correction)
        .map_err(PyRuntimeError::new_err)?;
    let tm = rsx_core::test_method::TestMethod::parse_str(test).map_err(PyRuntimeError::new_err)?;

    rsx_core::commands::signif::run(&rsx_core::commands::signif::SignifParams {
        markers_table_path: markers_tmp.path().to_string_lossy().to_string(),
        popmap_file_path: popmap_tmp.path().to_string_lossy().to_string(),
        output_file_path: out_path.clone(),
        min_depth,
        signif_threshold,
        correction: corr,
        test_method: tm,
        output_fasta,
        output_bayes: bayes,
        group1: group1.to_string(),
        group2: group2.to_string(),
    })
    .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

    read_tsv_to_pyarrow_table(py, &out_path)
}

// --------------------------------------------------------------------------- //
// PyO3 module registration
// --------------------------------------------------------------------------- //

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

    m.add_function(wrap_pyfunction!(triage_to_arrow, m)?)?;
    m.add_function(wrap_pyfunction!(pca_to_arrow, m)?)?;
    m.add_function(wrap_pyfunction!(triage_to_arrow_from_arrow, m)?)?;
    m.add_function(wrap_pyfunction!(pca_to_arrow_from_arrow, m)?)?;
    m.add_function(wrap_pyfunction!(freq_from_arrow, m)?)?;
    m.add_function(wrap_pyfunction!(depth_from_arrow, m)?)?;
    m.add_function(wrap_pyfunction!(distrib_from_arrow, m)?)?;
    m.add_function(wrap_pyfunction!(signif_from_arrow, m)?)?;

    Ok(())
}
