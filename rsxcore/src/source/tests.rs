//! Source-trait parity tests: the same fixture, iterated through
//! `MarkersTableStream` (file), `ArrowMarkerSource` (in-memory) and
//! `ParquetMarkerSource` (spill) must produce identical Marker views and
//! identical command outputs.

use std::io::Write;

use arrow::array::{Int32Array, StringArray};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;

use crate::commands;
use crate::markers_table::{MarkersTableStream, ParserConfig};
use crate::popmap::Popmap;
use crate::source::{ArrowMarkerSource, MarkerStream, ParquetMarkerSource};

fn fixture_tsv(dir: &std::path::Path) -> std::path::PathBuf {
    let table = dir.join("markers.tsv");
    let mut f = std::fs::File::create(&table).unwrap();
    writeln!(f, "#Number of markers : 5").unwrap();
    writeln!(f, "id\tsequence\tm1\tm2\tm3\tf1\tf2\tf3").unwrap();
    writeln!(f, "0\tALL\t10\t10\t10\t10\t10\t10").unwrap();
    writeln!(f, "1\tMONLY\t10\t10\t10\t0\t0\t0").unwrap();
    writeln!(f, "2\tFONLY\t0\t0\t0\t10\t10\t10").unwrap();
    writeln!(f, "3\tMIX\t10\t0\t10\t0\t10\t0").unwrap();
    writeln!(f, "4\tLOW\t1\t2\t1\t2\t1\t2").unwrap();
    table
}

fn fixture_popmap(dir: &std::path::Path) -> std::path::PathBuf {
    let pop = dir.join("popmap.tsv");
    let mut f = std::fs::File::create(&pop).unwrap();
    for ind in ["m1", "m2", "m3"] {
        writeln!(f, "{ind}\tM").unwrap();
    }
    for ind in ["f1", "f2", "f3"] {
        writeln!(f, "{ind}\tF").unwrap();
    }
    pop
}

fn fixture_batch() -> RecordBatch {
    let schema = Schema::new(vec![
        Field::new("id", DataType::Utf8, false),
        Field::new("sequence", DataType::Utf8, false),
        Field::new("m1", DataType::Int32, false),
        Field::new("m2", DataType::Int32, false),
        Field::new("m3", DataType::Int32, false),
        Field::new("f1", DataType::Int32, false),
        Field::new("f2", DataType::Int32, false),
        Field::new("f3", DataType::Int32, false),
    ]);
    let ids = StringArray::from(vec!["0", "1", "2", "3", "4"]);
    let seqs = StringArray::from(vec!["ALL", "MONLY", "FONLY", "MIX", "LOW"]);
    let m1 = Int32Array::from(vec![10, 10, 0, 10, 1]);
    let m2 = Int32Array::from(vec![10, 10, 0, 0, 2]);
    let m3 = Int32Array::from(vec![10, 10, 0, 10, 1]);
    let f1 = Int32Array::from(vec![10, 0, 10, 0, 2]);
    let f2 = Int32Array::from(vec![10, 0, 10, 10, 1]);
    let f3 = Int32Array::from(vec![10, 0, 10, 0, 2]);
    RecordBatch::try_new(
        std::sync::Arc::new(schema),
        vec![
            std::sync::Arc::new(ids),
            std::sync::Arc::new(seqs),
            std::sync::Arc::new(m1),
            std::sync::Arc::new(m2),
            std::sync::Arc::new(m3),
            std::sync::Arc::new(f1),
            std::sync::Arc::new(f2),
            std::sync::Arc::new(f3),
        ],
    )
    .unwrap()
}

/// Count qualifying markers (presence > 0) for `min_depth=1`.
fn count_qualifying<S: MarkerStream>(source: &S) -> u64 {
    source.count_markers().unwrap()
}

#[test]
fn file_arrow_parquet_marker_counts_agree() {
    let dir = tempfile::tempdir().unwrap();
    let table = fixture_tsv(dir.path());
    let popmap_path = fixture_popmap(dir.path());
    let popmap = Popmap::from_file(&popmap_path).unwrap();

    let config = ParserConfig {
        store_sequence: true,
        store_depths: true,
        compute_groups: true,
        min_depth: 1,
    };
    let file_stream = MarkersTableStream::open(&table, Some(&popmap), config).unwrap();
    let arrow = ArrowMarkerSource::from_batches(vec![fixture_batch()], Some(&popmap), 1).unwrap();
    let spilled = ParquetMarkerSource::spill_from_arrow(&arrow).unwrap();

    let n_file = count_qualifying(&file_stream);
    let n_arrow = count_qualifying(&arrow);
    let n_spilled = count_qualifying(&spilled);

    assert_eq!(n_file, n_arrow, "file vs arrow marker count");
    assert_eq!(n_file, n_spilled, "file vs spilled-parquet marker count");
}

#[test]
fn freq_outputs_match_across_sources() {
    let dir = tempfile::tempdir().unwrap();
    let table = fixture_tsv(dir.path());
    let popmap_path = fixture_popmap(dir.path());
    let popmap = Popmap::from_file(&popmap_path).unwrap();

    let arrow = ArrowMarkerSource::from_batches(vec![fixture_batch()], Some(&popmap), 1).unwrap();
    let spilled = ParquetMarkerSource::spill_from_arrow(&arrow).unwrap();

    let dir = std::env::temp_dir().join("rsx_source_parity_freq");
    std::fs::create_dir_all(&dir).unwrap();

    let file_out = dir.join("freq_file.tsv");
    let arrow_out = dir.join("freq_arrow.tsv");
    let spill_out = dir.join("freq_spill.tsv");

    commands::freq::run(&commands::freq::FreqParams {
        markers_table_path: table.to_str().unwrap().to_string(),
        output_file_path: file_out.to_str().unwrap().to_string(),
        min_depth: 1,
    })
    .unwrap();
    commands::freq::run_with_source(
        &arrow,
        &commands::freq::FreqParams {
            markers_table_path: String::new(),
            output_file_path: arrow_out.to_str().unwrap().to_string(),
            min_depth: 1,
        },
    )
    .unwrap();
    commands::freq::run_with_source(
        &spilled,
        &commands::freq::FreqParams {
            markers_table_path: String::new(),
            output_file_path: spill_out.to_str().unwrap().to_string(),
            min_depth: 1,
        },
    )
    .unwrap();

    let file_body = std::fs::read_to_string(&file_out).unwrap();
    let arrow_body = std::fs::read_to_string(&arrow_out).unwrap();
    let spill_body = std::fs::read_to_string(&spill_out).unwrap();

    assert_eq!(file_body, arrow_body, "arrow vs file freq output");
    assert_eq!(file_body, spill_body, "spill vs file freq output");
}
