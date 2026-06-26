// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Markers table parser: mmap + algorithmic optimizations.
//!
//! Key optimizations:
//! - mmap for zero-copy I/O
//! - memchr-based line and field iteration
//! - Specialized parser paths for presence-only, depth-only, and full marker rows
//! - Bitset presence tracking via popcount
//! - Optional parallel chunk processing behind the `parallel` feature

use crate::io::table_io::fast_parse_u16;
use crate::marker::Marker;
use crate::popmap::Popmap;

use memmap2::Mmap;
use std::path::Path;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Heuristic: only use parallel dispatch for tables large enough to amortize overhead.
#[cfg(feature = "parallel")]
#[inline]
fn should_use_parallel(len: usize) -> bool {
    len > 4 * 1024 * 1024
}

/// Configuration for the markers table parser (controls what is materialized
/// from the on-disk / Arrow source during streaming).
///
/// Used by all commands to trade off speed vs. functionality (e.g. fasta output
/// or per-marker sequences for mapping requires storing the sequence).
pub struct ParserConfig {
    pub store_sequence: bool,
    pub store_depths: bool,
    pub compute_groups: bool,
    pub min_depth: u16,
}

impl Default for ParserConfig {
    fn default() -> Self {
        ParserConfig {
            store_sequence: true,
            store_depths: true,
            compute_groups: true,
            min_depth: 1,
        }
    }
}

/// Streaming / bounded-memory view over a marker depth table (mmap or Arrow).
///
/// This is the central abstraction: commands consume rows without ever
/// materializing the full n_markers × n_ind matrix in RAM. See the paper and
/// `docs/orgmode/reference/architecture.org` for the design.
pub struct MarkersTableStream {
    pub header: crate::io::table_io::TableHeader,
    pub groups: Vec<String>,
    config: ParserConfig,
    mmap: Mmap,
    data_start: usize,
    n_individuals: u16,
}

impl MarkersTableStream {
    pub fn open(
        path: &Path,
        popmap: Option<&Popmap>,
        config: ParserConfig,
    ) -> std::io::Result<Self> {
        let header = crate::io::table_io::TableHeader::from_file(path)?;

        let groups: Vec<String> = if config.compute_groups {
            if let Some(pm) = popmap {
                let mut g = vec![String::new(), String::new()];
                for col in header.columns.iter().skip(2) {
                    g.push(pm.get_group(col).unwrap_or("").to_string());
                }
                g
            } else {
                Vec::new()
            }
        } else {
            Vec::new()
        };

        let file = std::fs::File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };

        let mut data_start = 0;
        loop {
            if data_start >= mmap.len() {
                break;
            }
            let line_end = memchr::memchr(b'\n', &mmap[data_start..])
                .map(|p| data_start + p + 1)
                .unwrap_or(mmap.len());
            if mmap[data_start] == b'#' || mmap[data_start..].starts_with(b"id\t") {
                data_start = line_end;
            } else {
                break;
            }
        }

        let n_individuals = header.n_individuals;

        Ok(MarkersTableStream {
            header,
            groups,
            config,
            mmap,
            data_start,
            n_individuals,
        })
    }

    /// Count markers with n_individuals > 0 (for Bonferroni correction).
    /// Streaming: O(1) memory, just counts.
    pub fn count_markers(&self) -> std::io::Result<u64> {
        let mut n = 0u64;
        self.for_each(|m| {
            if m.n_individuals > 0 {
                n += 1;
            }
        })?;
        Ok(n)
    }

    /// Process all markers. Uses fast path when sequence isn't needed.
    pub fn for_each<F>(&self, f: F) -> std::io::Result<()>
    where
        F: FnMut(&Marker),
    {
        let data = &self.mmap[self.data_start..];
        if data.is_empty() {
            return Ok(());
        }
        self.dispatch_on_slice(data, f)
    }

    /// Process all markers with a callback that is compatible with parallel execution.
    ///
    /// With the `parallel` feature enabled, callback order is not specified.
    #[cfg(feature = "parallel")]
    pub fn for_each_parallel<F>(&self, f: F) -> std::io::Result<()>
    where
        F: Fn(&Marker) + Send + Sync,
    {
        self.par_for_each(f)
    }

    #[cfg(not(feature = "parallel"))]
    pub fn for_each_parallel<F>(&self, f: F) -> std::io::Result<()>
    where
        F: Fn(&Marker) + Send + Sync,
    {
        self.for_each(|m| f(m))
    }

    /// Fast path: min_depth <= 1 and we don't need to store depths.
    /// Only tracks presence bits (no integer parsing).
    #[inline]
    fn for_each_fast_d1<F>(&self, data: &[u8], f: &mut F)
    where
        F: FnMut(&Marker),
    {
        let mut marker = Marker::new(self.n_individuals);

        for_each_line(data, |line| {
            let (_, tab2) = match find_first_two_tabs(line) {
                Some(t) => t,
                None => return,
            };

            for_each_depth_field(line, tab2, |col, field| {
                let is_zero = field.len() == 1 && field[0] == b'0';
                if !is_zero && !field.is_empty() && col < self.n_individuals as usize {
                    marker.presence.set(col);
                    marker.n_individuals += 1;
                }
            });

            f(&marker);
            marker.reset(true);
        });
    }

    /// Medium path: skip id+seq, parse depths for min_depth > 1.
    #[inline]
    fn for_each_skip_seq<F>(&self, data: &[u8], f: &mut F)
    where
        F: FnMut(&Marker),
    {
        let min_depth = self.config.min_depth;
        let mut marker = Marker::new(self.n_individuals);

        for_each_line(data, |line| {
            let (_, tab2) = match find_first_two_tabs(line) {
                Some(t) => t,
                None => return,
            };

            for_each_depth_field(line, tab2, |col, field| {
                let depth = fast_parse_u16(field);

                if col < self.n_individuals as usize {
                    if self.config.store_depths {
                        marker.individual_depths[col] = depth;
                    }
                    if depth >= min_depth {
                        marker.presence.set(col);
                        marker.n_individuals += 1;
                    }
                }
            });

            f(&marker);
            marker.reset(true);
        });
    }

    /// Full path: parse everything including id and sequence.
    #[inline]
    fn for_each_full<F>(&self, data: &[u8], f: &mut F)
    where
        F: FnMut(&Marker),
    {
        let min_depth = self.config.min_depth;
        let mut marker = Marker::new(self.n_individuals);

        for_each_line(data, |line| {
            let (tab1, tab2) = match find_first_two_tabs(line) {
                Some(t) => t,
                None => return,
            };

            marker.id.clear();
            marker
                .id
                .push_str(std::str::from_utf8(&line[0..tab1]).unwrap_or(""));

            let seq_start = tab1 + 1;
            let seq_end = tab2;
            marker.sequence.clear();
            marker
                .sequence
                .push_str(std::str::from_utf8(&line[seq_start..seq_end]).unwrap_or(""));

            for_each_depth_field(line, tab2, |col, field| {
                let depth = fast_parse_u16(field);
                if col < self.n_individuals as usize {
                    if self.config.store_depths {
                        marker.individual_depths[col] = depth;
                    }
                    if depth >= min_depth {
                        marker.presence.set(col);
                        marker.n_individuals += 1;
                    }
                }
            });

            f(&marker);
            marker.reset(false);
        });
    }

    pub fn collect(&self) -> std::io::Result<Vec<Marker>> {
        let mut markers = Vec::new();
        self.for_each(|m| markers.push(m.clone()))?;
        Ok(markers)
    }

    pub fn iter(&self) -> impl Iterator<Item = Marker> {
        self.collect().unwrap_or_default().into_iter()
    }

    /// Returns line-aligned chunks of the data for parallel processing.
    /// Target chunk size ~1 MiB for good granularity on large tables.
    #[cfg(feature = "parallel")]
    fn make_line_aligned_chunks<'a>(&self, data: &'a [u8], target: usize) -> Vec<&'a [u8]> {
        let mut chunks = Vec::new();
        let mut start = 0usize;
        while start < data.len() {
            let mut end = (start + target).min(data.len());
            if end < data.len() {
                if let Some(rel) = memchr::memchr(b'\n', &data[end..]) {
                    end = end + rel + 1;
                } else {
                    end = data.len();
                }
            }
            if end > start {
                chunks.push(&data[start..end]);
            }
            start = end;
        }
        chunks
    }

    /// Parallel for_each when the "parallel" feature is enabled.
    /// Splits the mmap into ~1 MiB line-aligned chunks and processes them
    /// concurrently with rayon. The closure must be `Send + Sync`.
    ///
    /// Use this for strong scaling on large marker tables (100k+ rows) on
    /// multi-core machines for commands like distrib, signif (non-FDR), freq, depth.
    ///
    /// For FDR in signif, the caller must use a thread-safe collector (e.g. DashMap
    /// or crossbeam channel + final sort by original order).
    #[cfg(feature = "parallel")]
    pub fn par_for_each<F>(&self, f: F) -> std::io::Result<()>
    where
        F: Fn(&Marker) + Send + Sync,
    {
        let data = &self.mmap[self.data_start..];
        if data.is_empty() {
            return Ok(());
        }

        if !should_use_parallel(data.len()) {
            return self.dispatch_on_slice(data, |m| f(m));
        }

        let target_chunk = 1 << 20;
        let chunks = self.make_line_aligned_chunks(data, target_chunk);

        chunks
            .into_par_iter()
            .try_for_each(|chunk| self.process_slice_serial(chunk, &f))
    }

    /// Collect mapped marker values in table order, using parallel parsing for large inputs.
    ///
    /// The filter/mapping closure runs independently per marker. Results are buffered per
    /// line-aligned chunk, then concatenated in chunk order before returning.
    #[cfg(feature = "parallel")]
    pub fn par_filter_map_collect<T, F>(&self, filter_map: F) -> std::io::Result<Vec<T>>
    where
        T: Send,
        F: Fn(&Marker) -> Option<T> + Send + Sync,
    {
        let data = &self.mmap[self.data_start..];
        if data.is_empty() {
            return Ok(Vec::new());
        }

        if !should_use_parallel(data.len()) {
            let mut collected = Vec::new();
            self.dispatch_on_slice(data, |marker| {
                if let Some(item) = filter_map(marker) {
                    collected.push(item);
                }
            })?;
            return Ok(collected);
        }

        let target_chunk = 1 << 20;
        let chunks = self.make_line_aligned_chunks(data, target_chunk);
        let chunked: Vec<std::io::Result<Vec<T>>> = chunks
            .into_par_iter()
            .map(|chunk| {
                let mut local = Vec::new();
                self.dispatch_on_slice(chunk, |marker| {
                    if let Some(item) = filter_map(marker) {
                        local.push(item);
                    }
                })
                .map(|()| local)
            })
            .collect();

        let mut collected = Vec::new();
        for chunk in chunked {
            let mut chunk = chunk?;
            collected.append(&mut chunk);
        }
        Ok(collected)
    }

    /// Internal dispatcher: run the appropriate specialized parser on a data slice,
    /// feeding every marker to the given visitor.
    fn dispatch_on_slice<F>(&self, slice: &[u8], mut visit: F) -> std::io::Result<()>
    where
        F: FnMut(&Marker),
    {
        if !self.config.store_sequence && !self.config.store_depths && self.config.min_depth <= 1 {
            self.for_each_fast_d1(slice, &mut visit);
        } else if !self.config.store_sequence {
            self.for_each_skip_seq(slice, &mut visit);
        } else {
            self.for_each_full(slice, &mut visit);
        }
        Ok(())
    }

    #[cfg(feature = "parallel")]
    fn process_slice_serial<F>(&self, slice: &[u8], f: &F) -> std::io::Result<()>
    where
        F: Fn(&Marker),
    {
        self.dispatch_on_slice(slice, |m| f(m))
    }

    /// Parallel fold + reduce for accumulation without per-marker locking.
    ///
    /// This is the ergonomic high-level API for strong scaling.
    ///
    /// Each rayon thread processes one or more chunks, maintaining a local `Acc`
    /// and calling `fold(&mut local, &marker)` for every marker. At the end the
    /// per-chunk accumulators are combined with `reduce`.
    ///
    /// This enables lock-free parallel accumulation for:
    /// - distrib 2D tables
    /// - per-individual depth/freq stats
    /// - FDR p-value collection (fold into Vec of (p, metadata))
    #[cfg(feature = "parallel")]
    pub fn par_fold_reduce<Acc, Fold, Reduce>(
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
        let data = &self.mmap[self.data_start..];
        if data.is_empty() {
            return Ok(init);
        }

        if !should_use_parallel(data.len()) {
            let mut local = init;
            self.dispatch_on_slice(data, |marker| {
                fold(&mut local, marker);
            })?;
            return Ok(local);
        }

        let target_chunk = 1 << 20;
        let chunks = self.make_line_aligned_chunks(data, target_chunk);

        let result = chunks
            .into_par_iter()
            .map(|chunk| {
                let mut local = init.clone();
                let _ = self.dispatch_on_slice(chunk, |marker| {
                    fold(&mut local, marker);
                });
                local
            })
            .reduce(|| init.clone(), reduce);

        Ok(result)
    }
}

/// Calls `f(col, field_bytes)` for every depth field after the second tab.
/// Uses `memchr_iter` for a single pass over the depth section.
#[inline]
fn for_each_depth_field<F>(line: &[u8], tab2: usize, mut f: F)
where
    F: FnMut(usize, &[u8]),
{
    let depth_start = tab2 + 1;
    if depth_start >= line.len() {
        return;
    }

    let depth_section = &line[depth_start..];
    let mut col = 0usize;
    let mut field_start = 0usize;

    for tab_pos in memchr::memchr_iter(b'\t', depth_section) {
        let field = &depth_section[field_start..tab_pos];
        f(col, field);
        col += 1;
        field_start = tab_pos + 1;
    }

    // Last field
    if field_start < depth_section.len() {
        let field = &depth_section[field_start..];
        f(col, field);
    }
}

/// Iterate over lines using `memchr_iter(b'\n')`.
/// Always strips trailing \r if present.
#[inline]
fn for_each_line<F>(data: &[u8], mut f: F)
where
    F: FnMut(&[u8]),
{
    let mut pos = 0;

    for end in memchr::memchr_iter(b'\n', data) {
        let raw = &data[pos..end];
        let line = raw.strip_suffix(b"\r").unwrap_or(raw);
        if !line.is_empty() {
            f(line);
        }
        pos = end + 1;
    }

    // Final line without trailing \n
    if pos < data.len() {
        let raw = &data[pos..];
        let line = raw.strip_suffix(b"\r").unwrap_or(raw);
        if !line.is_empty() {
            f(line);
        }
    }
}

/// Returns the byte offsets of the first and second tab in the line.
/// Returns `None` if there are fewer than two tabs.
#[inline]
fn find_first_two_tabs(line: &[u8]) -> Option<(usize, usize)> {
    let mut it = memchr::memchr_iter(b'\t', line);
    let t1 = it.next()?;
    let t2 = it.next()?;
    Some((t1, t2))
}
