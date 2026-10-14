// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Marker source abstraction.
//!
//! A `MarkerStream` is anything that can hand a stream of markers to the
//! analysis commands, regardless of where the data physically lives (mmap'd
//! TSV, in-memory Arrow RecordBatches, or a spilled Parquet temp file).
//!
//! The existing TSV reader (`crate::markers_table::MarkersTableStream`) gets
//! a blanket impl so the path-based code keeps working unchanged. The
//! optional `arrow-input` feature pulls in the in-memory Arrow source and
//! the Parquet spill source.

use crate::io::table_io::TableHeader;
use crate::marker::Marker;

#[cfg(feature = "arrow-input")]
pub mod arrow_source;
#[cfg(feature = "arrow-input")]
pub mod estimator;
#[cfg(feature = "arrow-input")]
pub mod parquet_source;

#[cfg(all(test, feature = "arrow-input"))]
mod tests;

/// In-memory Arrow marker source.
#[cfg(feature = "arrow-input")]
pub use arrow_source::ArrowMarkerSource;
/// Working-set estimate and the Arrow or Parquet source selected from it.
#[cfg(feature = "arrow-input")]
pub use estimator::{
    MarkerTableSource, SizeEstimate, estimate_working_set_bytes, spill_threshold_bytes,
};
/// Marker source backed by a temporary Parquet file.
#[cfg(feature = "arrow-input")]
pub use parquet_source::ParquetMarkerSource;

/// Stream of markers consumed by the analysis commands.
///
/// Implementations must surface:
///
/// - `header()` and `groups()`: same shape `MarkersTableStream` already
///   exposes (header columns include id + sequence + per-individual names,
///   groups vec has two placeholder entries followed by the popmap group
///   for each individual column).
/// - `count_markers()`: needed for the Bonferroni denominator.
/// - `for_each` / `par_for_each` / `par_fold_reduce`: the iteration shapes
///   the commands already use.
pub trait MarkerStream {
    /// Column names and individual count. Columns 0 and 1 are `id` and `sequence`.
    fn header(&self) -> &TableHeader;
    /// Popmap group of each column. Entries 0 and 1 are empty placeholders.
    fn groups(&self) -> &[String];
    /// Number of markers present in at least one individual.
    fn count_markers(&self) -> std::io::Result<u64>;

    /// Call `f` once per marker, in table order.
    fn for_each<F>(&self, f: F) -> std::io::Result<()>
    where
        F: FnMut(&Marker);

    /// Call `f` once per marker. The order is not specified.
    #[cfg(feature = "parallel")]
    fn par_for_each<F>(&self, f: F) -> std::io::Result<()>
    where
        F: Fn(&Marker) + Send + Sync;

    /// Fold markers in parallel and combine the partial accumulators.
    #[cfg(feature = "parallel")]
    fn par_fold_reduce<Acc, Fold, Reduce>(
        &self,
        init: Acc,
        fold: Fold,
        reduce: Reduce,
    ) -> std::io::Result<Acc>
    where
        Acc: Send + Sync + Clone,
        Fold: Fn(&mut Acc, &Marker) + Send + Sync + Clone,
        Reduce: Fn(Acc, Acc) -> Acc + Send + Sync;

    /// Optional columnar fast path for the `freq` histogram (count markers by
    /// number of individuals with depth ≥ min_depth).
    ///
    /// Default: `None` → callers use the generic `Marker` stream. Arrow sources
    /// override this to walk typed columns as tiles without materialising a
    /// `Marker` per row.
    fn freq_histogram_columnar(&self, _min_depth: u16) -> Option<std::io::Result<Vec<u32>>> {
        None
    }
}

impl MarkerStream for crate::markers_table::MarkersTableStream {
    fn header(&self) -> &TableHeader {
        &self.header
    }

    // freq_histogram_columnar: default (TSV path stays mmap + bitset stream)

    fn groups(&self) -> &[String] {
        &self.groups
    }

    fn count_markers(&self) -> std::io::Result<u64> {
        self.count_markers()
    }

    fn for_each<F>(&self, f: F) -> std::io::Result<()>
    where
        F: FnMut(&Marker),
    {
        self.for_each(f)
    }

    #[cfg(feature = "parallel")]
    fn par_for_each<F>(&self, f: F) -> std::io::Result<()>
    where
        F: Fn(&Marker) + Send + Sync,
    {
        self.par_for_each(f)
    }

    #[cfg(feature = "parallel")]
    fn par_fold_reduce<Acc, Fold, Reduce>(
        &self,
        init: Acc,
        fold: Fold,
        reduce: Reduce,
    ) -> std::io::Result<Acc>
    where
        Acc: Send + Sync + Clone,
        Fold: Fn(&mut Acc, &Marker) + Send + Sync + Clone,
        Reduce: Fn(Acc, Acc) -> Acc + Send + Sync,
    {
        self.par_fold_reduce(init, fold, reduce)
    }
}
