// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! `ParquetMarkerSource`: marker stream backed by a Parquet temp file.
//!
//! Used when the predicted working set would not safely fit in RAM. The
//! Parquet payload has the same schema as the in-memory variant
//! (`id`, `sequence`, `ind1`, ..., `indN`) so the per-row decoding is
//! shared via [`super::arrow_source::ArrowMarkerSource::from_batches`].

use std::fs::File;
use std::io;
use std::path::PathBuf;

use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

use crate::io::table_io::TableHeader;
use crate::marker::Marker;

use super::MarkerStream;
use super::arrow_source::{ArrowMarkerSource, ArrowSourceError};

/// Spilled-to-Parquet marker source. Owns the underlying temp file so it
/// is removed automatically when the source is dropped.
pub struct ParquetMarkerSource {
    /// Held only to keep the temp file alive for as long as the source
    /// is reachable; the actual reads go through the path below.
    _tempfile: tempfile::NamedTempFile,
    path: PathBuf,
    header: TableHeader,
    groups: Vec<String>,
    min_depth: u16,
}

impl ParquetMarkerSource {
    /// Materialise an in-memory [`ArrowMarkerSource`] as a Parquet temp
    /// file and wrap it as a `ParquetMarkerSource`. The original batches
    /// are dropped once the file is written.
    pub fn spill_from_arrow(in_mem: &ArrowMarkerSource) -> Result<Self, ParquetSourceError> {
        use parquet::arrow::ArrowWriter;
        use parquet::file::properties::WriterProperties;

        let tempfile = tempfile::NamedTempFile::new()
            .map_err(|e| ParquetSourceError::TempFile(e.to_string()))?;
        let path = tempfile.path().to_path_buf();

        let schema = if let Some(b) = in_mem.batches().first() {
            b.schema()
        } else {
            return Err(ParquetSourceError::EmptyInput);
        };

        let file = File::create(&path).map_err(|e| ParquetSourceError::Write(e.to_string()))?;
        let props = WriterProperties::builder().build();
        let mut writer = ArrowWriter::try_new(file, schema, Some(props))
            .map_err(|e| ParquetSourceError::Write(e.to_string()))?;
        for batch in in_mem.batches() {
            writer
                .write(batch)
                .map_err(|e| ParquetSourceError::Write(e.to_string()))?;
        }
        writer
            .close()
            .map_err(|e| ParquetSourceError::Write(e.to_string()))?;

        Ok(Self {
            _tempfile: tempfile,
            path,
            header: in_mem.header().clone(),
            groups: in_mem.groups().to_vec(),
            min_depth: in_mem.min_depth(),
        })
    }

    fn open_batches(&self) -> Result<Vec<RecordBatch>, io::Error> {
        let file = File::open(&self.path)?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(file)
            .map_err(|e| io::Error::other(format!("parquet open: {e}")))?;
        let reader = builder
            .build()
            .map_err(|e| io::Error::other(format!("parquet build: {e}")))?;
        let mut batches = Vec::new();
        for b in reader {
            batches.push(b.map_err(|e| io::Error::other(format!("parquet batch: {e}")))?);
        }
        Ok(batches)
    }

    /// Path of the on-disk Parquet payload (e.g. for diagnostics).
    pub fn path(&self) -> &std::path::Path {
        &self.path
    }
}

impl MarkerStream for ParquetMarkerSource {
    fn header(&self) -> &TableHeader {
        &self.header
    }
    fn groups(&self) -> &[String] {
        &self.groups
    }

    fn count_markers(&self) -> io::Result<u64> {
        // The header carries a marker count, but it is an *upper bound*
        // because n_individuals=0 rows are not counted by the CLI. Re-derive
        // by streaming so we honour the same definition as the file source.
        let batches = self.open_batches()?;
        let arrow = ArrowMarkerSource::from_batches(batches, None, self.min_depth)
            .map_err(map_arrow_err)?;
        arrow.count_markers()
    }

    fn for_each<F>(&self, f: F) -> io::Result<()>
    where
        F: FnMut(&Marker),
    {
        let batches = self.open_batches()?;
        let arrow = ArrowMarkerSource::from_batches(batches, None, self.min_depth)
            .map_err(map_arrow_err)?;
        // Group annotations live on this struct, not on the temporary arrow source.
        arrow.for_each(f)
    }

    #[cfg(feature = "parallel")]
    fn par_for_each<F>(&self, f: F) -> io::Result<()>
    where
        F: Fn(&Marker) + Send + Sync,
    {
        let batches = self.open_batches()?;
        let arrow = ArrowMarkerSource::from_batches(batches, None, self.min_depth)
            .map_err(map_arrow_err)?;
        arrow.par_for_each(f)
    }

    #[cfg(feature = "parallel")]
    fn par_fold_reduce<Acc, Fold, Reduce>(
        &self,
        init: Acc,
        fold: Fold,
        reduce: Reduce,
    ) -> io::Result<Acc>
    where
        Acc: Send + Sync + Clone,
        Fold: Fn(&mut Acc, &Marker) + Send + Sync + Clone,
        Reduce: Fn(Acc, Acc) -> Acc + Send + Sync,
    {
        let batches = self.open_batches()?;
        let arrow = ArrowMarkerSource::from_batches(batches, None, self.min_depth)
            .map_err(map_arrow_err)?;
        arrow.par_fold_reduce(init, fold, reduce)
    }
}

fn map_arrow_err(e: ArrowSourceError) -> io::Error {
    io::Error::other(format!("arrow source: {e}"))
}

#[derive(Debug)]
pub enum ParquetSourceError {
    EmptyInput,
    TempFile(String),
    Write(String),
}

impl std::fmt::Display for ParquetSourceError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptyInput => write!(f, "cannot spill empty Arrow source to Parquet"),
            Self::TempFile(s) => write!(f, "tempfile: {s}"),
            Self::Write(s) => write!(f, "parquet write: {s}"),
        }
    }
}

impl std::error::Error for ParquetSourceError {}
