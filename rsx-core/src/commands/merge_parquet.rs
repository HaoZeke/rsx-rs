// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Parquet output for the merge command.
//!
//! Writes the merged marker depth table as a Parquet file with
//! ZSTD compression. Schema: id (UInt64), sequence (Utf8),
//! sample_1..sample_N (UInt16, nullable).

#[cfg(feature = "parquet-io")]
use arrow::array::{ArrayRef, StringArray, UInt16Array, UInt64Array};
#[cfg(feature = "parquet-io")]
use arrow::datatypes::{DataType, Field, Schema};
#[cfg(feature = "parquet-io")]
use arrow::record_batch::RecordBatch;
#[cfg(feature = "parquet-io")]
use parquet::arrow::ArrowWriter;
#[cfg(feature = "parquet-io")]
use parquet::basic::Compression;
#[cfg(feature = "parquet-io")]
use parquet::file::properties::WriterProperties;

#[cfg(feature = "parquet-io")]
use crate::io::seq_reader::unpack_2bit;
#[cfg(feature = "parquet-io")]
use std::sync::Arc;

/// Write a batch of merged markers to a Parquet file.
/// Called from `merge::run` when `--output-parquet` is set.
#[cfg(feature = "parquet-io")]
pub fn write_parquet(
    path: &str,
    sample_names: &[String],
    rows: impl Iterator<Item = (u64, Vec<u8>, Vec<u16>)>, // (id, packed_seq, depths)
) -> Result<(), Box<dyn std::error::Error>> {
    let n_samples = sample_names.len();

    // Build schema dynamically
    let mut fields = vec![
        Field::new("id", DataType::UInt64, false),
        Field::new("sequence", DataType::Utf8, false),
    ];
    for name in sample_names {
        fields.push(Field::new(name, DataType::UInt16, true));
    }
    let schema = Arc::new(Schema::new(fields));

    // Writer with ZSTD compression
    let file = std::fs::File::create(path)?;
    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        .set_max_row_group_size(1_000_000)
        .build();
    let mut writer = ArrowWriter::try_new(file, schema.clone(), Some(props))?;

    // Accumulate rows into batches
    const BATCH_SIZE: usize = 100_000;
    let mut ids: Vec<u64> = Vec::with_capacity(BATCH_SIZE);
    let mut seqs: Vec<String> = Vec::with_capacity(BATCH_SIZE);
    let mut depth_cols: Vec<Vec<Option<u16>>> = (0..n_samples)
        .map(|_| Vec::with_capacity(BATCH_SIZE))
        .collect();

    let flush = |writer: &mut ArrowWriter<std::fs::File>,
                 ids: &mut Vec<u64>,
                 seqs: &mut Vec<String>,
                 depth_cols: &mut Vec<Vec<Option<u16>>>,
                 schema: &Arc<Schema>|
     -> Result<(), Box<dyn std::error::Error>> {
        if ids.is_empty() {
            return Ok(());
        }

        let mut arrays: Vec<ArrayRef> = vec![
            Arc::new(UInt64Array::from(ids.clone())),
            Arc::new(StringArray::from(seqs.clone())),
        ];
        for col in depth_cols.iter() {
            arrays.push(Arc::new(UInt16Array::from(col.clone())));
        }

        let batch = RecordBatch::try_new(schema.clone(), arrays)?;
        writer.write(&batch)?;

        ids.clear();
        seqs.clear();
        for col in depth_cols.iter_mut() {
            col.clear();
        }
        Ok(())
    };

    for (id, packed_seq, depths) in rows {
        ids.push(id);
        let unpacked = unpack_2bit(&packed_seq);
        seqs.push(String::from_utf8(unpacked).unwrap_or_else(|_| "?".to_string()));

        for (j, col) in depth_cols.iter_mut().enumerate() {
            let d = depths.get(j).copied().unwrap_or(0);
            col.push(if d == 0 { None } else { Some(d) });
        }

        if ids.len() >= BATCH_SIZE {
            flush(&mut writer, &mut ids, &mut seqs, &mut depth_cols, &schema)?;
        }
    }

    flush(&mut writer, &mut ids, &mut seqs, &mut depth_cols, &schema)?;
    writer.close()?;

    Ok(())
}

/// Stub for non-parquet builds.
#[cfg(not(feature = "parquet-io"))]
pub fn write_parquet(
    _path: &str,
    _sample_names: &[String],
    _rows: impl Iterator<Item = (u64, Vec<u8>, Vec<u16>)>,
) -> Result<(), Box<dyn std::error::Error>> {
    Err("Parquet output requires --features parquet-io".into())
}
