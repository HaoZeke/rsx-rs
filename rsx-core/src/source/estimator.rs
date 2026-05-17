// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Working-set size estimator and spill decision.
//!
//! The estimator combines the *observed* marker count from the inbound
//! Arrow IPC payload with a conservative overhead multiplier drawn from
//! the RAD-seq literature (Beissinger 2013, TASSEL-GBS, ipyrad).
//!
//! `MarkerTableSource::from_arrow_ipc` decodes the bytes once, asks the
//! estimator whether the implied working set fits, and either keeps the
//! batches in RAM (`InMemory(ArrowMarkerSource)`) or spills them to a
//! Parquet temp file (`Spilled(ParquetMarkerSource)`).

use std::io;

use crate::io::table_io::TableHeader;
use crate::marker::Marker;
use crate::popmap::Popmap;

use super::arrow_source::{ArrowMarkerSource, ArrowSourceError};
use super::parquet_source::{ParquetMarkerSource, ParquetSourceError};
use super::MarkerStream;

/// Bytes per depth cell. We always store u16 in the marker buffer
/// regardless of the inbound Arrow type, so this is a fixed 2 bytes.
pub const BYTES_PER_CELL: u64 = 2;

/// Default overhead multiplier capturing arrow validity buffers, group
/// masks, per-marker accumulators, intermediate Vecs. 6x is conservative
/// for the largest commands (signif FDR, triage Bayesian, depth exact).
pub const DEFAULT_OVERHEAD: f64 = 6.0;

/// Default fraction of available RAM we are willing to use before
/// switching to the spill path.
pub const DEFAULT_SPILL_FRACTION: f64 = 0.55;

/// Working-set estimate.
#[derive(Debug, Clone, Copy)]
pub struct SizeEstimate {
    pub n_samples: usize,
    pub m_markers: usize,
    pub bytes_per_cell: u64,
    pub overhead_factor: f64,
    pub command_multiplier: f64,
    pub estimated_bytes: u64,
}

/// Compute the predicted working-set size in bytes.
///
/// `command_specific_multiplier` lets the caller widen the prediction for
/// the heavier commands (e.g. 2.0 for triage / signif with FDR, 1.3 for
/// freq / depth which mostly stream).
pub fn estimate_working_set_bytes(
    n_samples: usize,
    m_observed_or_predicted: usize,
    bytes_per_cell: u64,
    overhead_factor: f64,
    command_specific_multiplier: f64,
) -> SizeEstimate {
    let raw = (n_samples as u64) * (m_observed_or_predicted as u64) * bytes_per_cell;
    let estimated = (raw as f64 * overhead_factor * command_specific_multiplier) as u64;
    SizeEstimate {
        n_samples,
        m_markers: m_observed_or_predicted,
        bytes_per_cell,
        overhead_factor,
        command_multiplier: command_specific_multiplier,
        estimated_bytes: estimated,
    }
}

/// Available physical memory in bytes. Falls back to a 4 GiB assumption
/// when the platform query fails so we never panic; the caller can always
/// override via `RSX_SPILL_BYTES` for deterministic tests.
fn available_memory_bytes() -> u64 {
    // /proc/meminfo is fine for the linux-only deployment surface this
    // crate targets; reading it avoids pulling in a new dependency
    // (sysinfo/lockfree-object-pool) for a single integer.
    if let Ok(contents) = std::fs::read_to_string("/proc/meminfo") {
        for line in contents.lines() {
            if let Some(rest) = line.strip_prefix("MemAvailable:") {
                if let Some(kb_str) = rest.split_whitespace().next() {
                    if let Ok(kb) = kb_str.parse::<u64>() {
                        return kb.saturating_mul(1024);
                    }
                }
            }
        }
    }
    4 * 1024 * 1024 * 1024
}

/// Bytes above which we should spill rather than keep the source in RAM.
pub fn spill_threshold_bytes() -> u64 {
    if let Ok(raw) = std::env::var("RSX_SPILL_BYTES") {
        if let Ok(v) = raw.parse::<u64>() {
            return v;
        }
    }
    let avail = available_memory_bytes();
    (avail as f64 * DEFAULT_SPILL_FRACTION) as u64
}

/// Resolved marker source: either in-memory Arrow or a Parquet spill.
///
/// Wraps the underlying source so the analysis commands can stay generic
/// over `MarkerStream` without caring which physical backing they got.
pub enum MarkerTableSource {
    InMemory(ArrowMarkerSource),
    Spilled(ParquetMarkerSource),
}

impl MarkerTableSource {
    /// Decode the inbound IPC bytes, consult the estimator, and produce
    /// either an in-memory or spilled source. The popmap is optional but
    /// required by the multi-group commands (distrib/signif/triage/depth).
    pub fn from_arrow_ipc(
        bytes: &[u8],
        popmap: Option<&Popmap>,
        min_depth: u16,
        command_multiplier: f64,
    ) -> Result<Self, MarkerSourceError> {
        let in_mem = ArrowMarkerSource::from_ipc_bytes(bytes, popmap, min_depth)
            .map_err(MarkerSourceError::Arrow)?;

        let estimate = estimate_working_set_bytes(
            in_mem.n_individuals() as usize,
            in_mem.n_markers() as usize,
            BYTES_PER_CELL,
            DEFAULT_OVERHEAD,
            command_multiplier,
        );

        let threshold = spill_threshold_bytes();
        let force = std::env::var("RSX_FORCE_SPILL").ok();
        let force_spill = matches!(force.as_deref(), Some("1") | Some("true") | Some("yes"));
        let force_inmem = matches!(force.as_deref(), Some("0") | Some("false") | Some("no") | Some("never"));

        if force_inmem {
            log::debug!("MarkerTableSource: in-memory forced via RSX_FORCE_SPILL");
            return Ok(MarkerTableSource::InMemory(in_mem));
        }

        if force_spill || estimate.estimated_bytes > threshold {
            log::info!(
                "MarkerTableSource: spilling to Parquet (estimated {} bytes > threshold {})",
                estimate.estimated_bytes, threshold
            );
            let spilled = ParquetMarkerSource::spill_from_arrow(&in_mem)
                .map_err(MarkerSourceError::Parquet)?;
            Ok(MarkerTableSource::Spilled(spilled))
        } else {
            Ok(MarkerTableSource::InMemory(in_mem))
        }
    }

    /// Convenience: was this source materialised to disk?
    pub fn is_spilled(&self) -> bool {
        matches!(self, MarkerTableSource::Spilled(_))
    }
}

impl MarkerStream for MarkerTableSource {
    fn header(&self) -> &TableHeader {
        match self {
            MarkerTableSource::InMemory(s) => s.header(),
            MarkerTableSource::Spilled(s) => s.header(),
        }
    }
    fn groups(&self) -> &[String] {
        match self {
            MarkerTableSource::InMemory(s) => s.groups(),
            MarkerTableSource::Spilled(s) => s.groups(),
        }
    }
    fn count_markers(&self) -> io::Result<u64> {
        match self {
            MarkerTableSource::InMemory(s) => s.count_markers(),
            MarkerTableSource::Spilled(s) => s.count_markers(),
        }
    }
    fn for_each<F>(&self, f: F) -> io::Result<()>
    where
        F: FnMut(&Marker),
    {
        match self {
            MarkerTableSource::InMemory(s) => s.for_each(f),
            MarkerTableSource::Spilled(s) => s.for_each(f),
        }
    }

    #[cfg(feature = "parallel")]
    fn par_for_each<F>(&self, f: F) -> io::Result<()>
    where
        F: Fn(&Marker) + Send + Sync,
    {
        match self {
            MarkerTableSource::InMemory(s) => s.par_for_each(f),
            MarkerTableSource::Spilled(s) => s.par_for_each(f),
        }
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
        match self {
            MarkerTableSource::InMemory(s) => s.par_fold_reduce(init, fold, reduce),
            MarkerTableSource::Spilled(s) => s.par_fold_reduce(init, fold, reduce),
        }
    }
}

/// Combined error for `from_arrow_ipc`.
#[derive(Debug)]
pub enum MarkerSourceError {
    Arrow(ArrowSourceError),
    Parquet(ParquetSourceError),
}

impl std::fmt::Display for MarkerSourceError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Arrow(e) => write!(f, "{e}"),
            Self::Parquet(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for MarkerSourceError {}
