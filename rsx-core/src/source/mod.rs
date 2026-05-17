// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Marker source abstraction.
//!
//! A `MarkerStream` is anything that can hand a stream of `Marker`s to the
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
pub mod parquet_source;
#[cfg(feature = "arrow-input")]
pub mod estimator;

#[cfg(all(test, feature = "arrow-input"))]
mod tests;

#[cfg(feature = "arrow-input")]
pub use arrow_source::ArrowMarkerSource;
#[cfg(feature = "arrow-input")]
pub use parquet_source::ParquetMarkerSource;
#[cfg(feature = "arrow-input")]
pub use estimator::{
    estimate_working_set_bytes, spill_threshold_bytes, MarkerTableSource, SizeEstimate,
};

/// Stream of markers consumed by the analysis commands.
///
/// Implementations must surface:
/// - `header()` and `groups()`: same shape `MarkersTableStream` already
///   exposes (header columns include id + sequence + per-individual names,
///   groups vec has two placeholder entries followed by the popmap group
///   for each individual column).
/// - `count_markers()`: needed for the Bonferroni denominator.
/// - `for_each` / `par_for_each` / `par_fold_reduce`: the iteration shapes
///   the commands already use.
pub trait MarkerStream {
    fn header(&self) -> &TableHeader;
    fn groups(&self) -> &[String];
    fn count_markers(&self) -> std::io::Result<u64>;

    fn for_each<F>(&self, f: F) -> std::io::Result<()>
    where
        F: FnMut(&Marker);

    #[cfg(feature = "parallel")]
    fn par_for_each<F>(&self, f: F) -> std::io::Result<()>
    where
        F: Fn(&Marker) + Send + Sync;

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
}

impl MarkerStream for crate::markers_table::MarkersTableStream {
    fn header(&self) -> &TableHeader {
        &self.header
    }

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
