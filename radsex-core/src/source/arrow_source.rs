// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! `ArrowMarkerSource`: drives the analysis commands directly from
//! in-memory Arrow `RecordBatch`es with no TSV round-trip.
//!
//! Built from Arrow IPC bytes, the source synthesises the same
//! `TableHeader` / `groups` / `Marker` shapes the path-based reader
//! produces so the existing command logic does not have to care which
//! source it received.

use std::io;
use std::sync::Arc;

use arrow::array::Array;
use arrow::datatypes::SchemaRef;
use arrow::record_batch::RecordBatch;

use crate::io::table_io::TableHeader;
use crate::marker::Marker;
use crate::popmap::Popmap;

use super::MarkerStream;

/// In-memory Arrow source. Holds Arc'd RecordBatches so cloning is cheap
/// (parallel iteration shares the underlying buffers).
pub struct ArrowMarkerSource {
    header: TableHeader,
    groups: Vec<String>,
    min_depth: u16,
    batches: Arc<[RecordBatch]>,
}

impl ArrowMarkerSource {
    /// Build a source from already-decoded RecordBatches.
    /// All batches must share the same schema (id, sequence, ind1, …, indN).
    pub fn from_batches(
        batches: Vec<RecordBatch>,
        popmap: Option<&Popmap>,
        min_depth: u16,
    ) -> Result<Self, ArrowSourceError> {
        let schema: SchemaRef = if let Some(b) = batches.first() {
            b.schema()
        } else {
            return Err(ArrowSourceError::EmptyInput);
        };
        for b in &batches {
            if b.schema() != schema {
                return Err(ArrowSourceError::HeterogeneousSchemas);
            }
        }

        if schema.fields().len() < 3 {
            return Err(ArrowSourceError::ShortSchema(schema.fields().len()));
        }

        let columns: Vec<String> = schema.fields().iter().map(|f| f.name().clone()).collect();
        let n_individuals = (columns.len() - 2) as u16;
        let n_markers: u64 = batches.iter().map(|b| b.num_rows() as u64).sum();

        let header = TableHeader {
            n_markers,
            n_individuals,
            columns,
        };

        let groups: Vec<String> = if let Some(pm) = popmap {
            let mut g = vec![String::new(), String::new()];
            for col in header.columns.iter().skip(2) {
                g.push(pm.get_group(col).unwrap_or("").to_string());
            }
            g
        } else {
            Vec::new()
        };

        Ok(Self {
            header,
            groups,
            min_depth,
            batches: batches.into(),
        })
    }

    /// Decode Arrow IPC stream bytes into an in-memory source.
    pub fn from_ipc_bytes(
        bytes: &[u8],
        popmap: Option<&Popmap>,
        min_depth: u16,
    ) -> Result<Self, ArrowSourceError> {
        let cursor = io::Cursor::new(bytes);
        let reader = arrow::ipc::reader::StreamReader::try_new(cursor, None)
            .map_err(|e| ArrowSourceError::IpcRead(e.to_string()))?;
        let mut batches = Vec::new();
        for b in reader {
            batches.push(b.map_err(|e| ArrowSourceError::IpcRead(e.to_string()))?);
        }
        Self::from_batches(batches, popmap, min_depth)
    }

    /// Number of markers across all batches (cached on the header).
    pub fn n_markers(&self) -> u64 {
        self.header.n_markers
    }

    /// Number of individual depth columns.
    pub fn n_individuals(&self) -> u16 {
        self.header.n_individuals
    }

    /// Borrow the underlying batches (e.g. for re-spilling to Parquet).
    pub fn batches(&self) -> &[RecordBatch] {
        &self.batches
    }

    /// Minimum depth threshold for filtering presence bits during iteration.
    pub fn min_depth(&self) -> u16 {
        self.min_depth
    }

    /// Materialise one `Marker` per row of the supplied batches.
    /// Shared between serial and parallel paths.
    fn for_each_in_batches<F>(batches: &[RecordBatch], n_individuals: u16, min_depth: u16, mut f: F)
    where
        F: FnMut(&Marker),
    {
        let mut marker = Marker::new(n_individuals);
        for batch in batches {
            let cols: Vec<&dyn Array> = (0..batch.num_columns())
                .map(|i| batch.column(i).as_ref())
                .collect();

            for row in 0..batch.num_rows() {
                marker.reset(false);
                marker.id.push_str(&array_value_as_string(cols[0], row));
                marker
                    .sequence
                    .push_str(&array_value_as_string(cols[1], row));

                for ind_idx in 0..n_individuals as usize {
                    let col = cols[ind_idx + 2];
                    let depth = array_value_as_u16(col, row);
                    marker.individual_depths[ind_idx] = depth;
                    if depth >= min_depth {
                        marker.presence.set(ind_idx);
                        marker.n_individuals += 1;
                    }
                }

                f(&marker);
            }
        }
    }
}

impl MarkerStream for ArrowMarkerSource {
    fn header(&self) -> &TableHeader {
        &self.header
    }

    fn groups(&self) -> &[String] {
        &self.groups
    }

    fn count_markers(&self) -> io::Result<u64> {
        let mut n = 0u64;
        Self::for_each_in_batches(
            &self.batches,
            self.header.n_individuals,
            self.min_depth,
            |m| {
                if m.n_individuals > 0 {
                    n += 1;
                }
            },
        );
        Ok(n)
    }

    fn for_each<F>(&self, f: F) -> io::Result<()>
    where
        F: FnMut(&Marker),
    {
        Self::for_each_in_batches(&self.batches, self.header.n_individuals, self.min_depth, f);
        Ok(())
    }

    #[cfg(feature = "parallel")]
    fn par_for_each<F>(&self, f: F) -> io::Result<()>
    where
        F: Fn(&Marker) + Send + Sync,
    {
        use rayon::prelude::*;
        let n_ind = self.header.n_individuals;
        let min_depth = self.min_depth;
        self.batches.par_iter().for_each(|batch| {
            let single = std::slice::from_ref(batch);
            Self::for_each_in_batches(single, n_ind, min_depth, |m| f(m));
        });
        Ok(())
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
        use rayon::prelude::*;
        let n_ind = self.header.n_individuals;
        let min_depth = self.min_depth;

        if self.batches.is_empty() {
            return Ok(init);
        }

        let init_for_reduce = init.clone();
        let res = self
            .batches
            .par_iter()
            .map(|batch| {
                let mut local = init.clone();
                let single = std::slice::from_ref(batch);
                Self::for_each_in_batches(single, n_ind, min_depth, |m| fold(&mut local, m));
                local
            })
            .reduce(|| init_for_reduce.clone(), reduce);
        Ok(res)
    }
}

#[derive(Debug)]
pub enum ArrowSourceError {
    EmptyInput,
    HeterogeneousSchemas,
    ShortSchema(usize),
    IpcRead(String),
}

impl std::fmt::Display for ArrowSourceError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptyInput => write!(f, "empty Arrow IPC payload"),
            Self::HeterogeneousSchemas => {
                write!(f, "RecordBatches in the payload have non-matching schemas")
            }
            Self::ShortSchema(n) => write!(
                f,
                "markers table needs id, sequence, and >=1 individual; got {n} columns"
            ),
            Self::IpcRead(s) => write!(f, "Arrow IPC read failed: {s}"),
        }
    }
}

impl std::error::Error for ArrowSourceError {}

/// Coerce an Arrow column cell into a u16 depth value.
///
/// - Integer types saturate at `u16::MAX` (matches the TSV parser).
/// - Floats round and clamp into u16.
/// - Strings parse as decimal; junk → 0.
/// - Anything else (including null) → 0.
fn array_value_as_u16(array: &dyn Array, row: usize) -> u16 {
    use arrow::array::*;
    use arrow::datatypes::DataType;

    if array.is_null(row) {
        return 0;
    }

    macro_rules! clamp_int {
        ($arr:ty) => {{
            let v = array.as_any().downcast_ref::<$arr>().unwrap().value(row) as i64;
            if v <= 0 {
                0
            } else if v >= u16::MAX as i64 {
                u16::MAX
            } else {
                v as u16
            }
        }};
    }
    macro_rules! clamp_uint {
        ($arr:ty) => {{
            let v = array.as_any().downcast_ref::<$arr>().unwrap().value(row) as u64;
            if v >= u16::MAX as u64 {
                u16::MAX
            } else {
                v as u16
            }
        }};
    }

    match array.data_type() {
        DataType::Int8 => clamp_int!(Int8Array),
        DataType::Int16 => clamp_int!(Int16Array),
        DataType::Int32 => clamp_int!(Int32Array),
        DataType::Int64 => clamp_int!(Int64Array),
        DataType::UInt8 => clamp_uint!(UInt8Array),
        DataType::UInt16 => clamp_uint!(UInt16Array),
        DataType::UInt32 => clamp_uint!(UInt32Array),
        DataType::UInt64 => clamp_uint!(UInt64Array),
        DataType::Float32 => {
            let v = array
                .as_any()
                .downcast_ref::<Float32Array>()
                .unwrap()
                .value(row);
            if !v.is_finite() || v <= 0.0 {
                0
            } else if v >= u16::MAX as f32 {
                u16::MAX
            } else {
                v.round() as u16
            }
        }
        DataType::Float64 => {
            let v = array
                .as_any()
                .downcast_ref::<Float64Array>()
                .unwrap()
                .value(row);
            if !v.is_finite() || v <= 0.0 {
                0
            } else if v >= u16::MAX as f64 {
                u16::MAX
            } else {
                v.round() as u16
            }
        }
        DataType::Utf8 => array
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap()
            .value(row)
            .parse()
            .unwrap_or(0),
        DataType::LargeUtf8 => array
            .as_any()
            .downcast_ref::<LargeStringArray>()
            .unwrap()
            .value(row)
            .parse()
            .unwrap_or(0),
        _ => 0,
    }
}

/// Coerce an Arrow column cell into a TSV-safe string.
fn array_value_as_string(array: &dyn Array, row: usize) -> String {
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
        DataType::Int8 => array
            .as_any()
            .downcast_ref::<Int8Array>()
            .unwrap()
            .value(row)
            .to_string(),
        DataType::Int16 => array
            .as_any()
            .downcast_ref::<Int16Array>()
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
        DataType::UInt8 => array
            .as_any()
            .downcast_ref::<UInt8Array>()
            .unwrap()
            .value(row)
            .to_string(),
        DataType::UInt16 => array
            .as_any()
            .downcast_ref::<UInt16Array>()
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
        DataType::Float32 => {
            let v = array
                .as_any()
                .downcast_ref::<Float32Array>()
                .unwrap()
                .value(row);
            if v.is_finite() {
                (v.round() as i64).to_string()
            } else {
                String::new()
            }
        }
        DataType::Float64 => {
            let v = array
                .as_any()
                .downcast_ref::<Float64Array>()
                .unwrap()
                .value(row);
            if v.is_finite() {
                (v.round() as i64).to_string()
            } else {
                String::new()
            }
        }
        _ => String::new(),
    }
}
