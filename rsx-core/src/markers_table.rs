// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Markers table parser: mmap + algorithmic optimizations.
//!
//! Key optimizations:
//! - mmap for zero-copy I/O
//! - Skip id+sequence fields when not needed (memchr to 2nd tab)
//! - Fast min_depth=1 threshold: check "0" vs non-"0" (1 byte, no parse)
//! - Bitset presence tracking via popcount

use crate::io::table_io::fast_parse_u16;
use crate::marker::Marker;
use crate::popmap::Popmap;

use memmap2::Mmap;
use std::path::Path;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Configuration for the markers table parser.
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

/// Markers table backed by mmap.
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

    /// Fast path: min_depth=1, no sequence or depth values needed.
    /// Skip id+seq fields. Check "0" vs non-"0" without integer parse.
    /// Uses one delimiter pass over each line.
    #[inline]
    fn for_each_fast_d1<F>(&self, data: &[u8], f: &mut F)
    where
        F: FnMut(&Marker),
    {
        let mut marker = Marker::new(self.n_individuals);
        let mut pos = 0;

        while pos < data.len() {
            let line_end = memchr::memchr(b'\n', &data[pos..])
                .map(|p| pos + p)
                .unwrap_or(data.len());

            let line = strip_cr(&data[pos..line_end]);
            pos = line_end + 1;

            // The second tab starts the depth columns.
            let mut tab_iter = memchr::memchr_iter(b'\t', line);
            let _tab1 = tab_iter.next();
            let tab2 = match tab_iter.next() {
                Some(p) => p,
                None => continue,
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
        }
    }

    /// Medium path: skip id+seq, parse depths for min_depth > 1.
    #[inline]
    fn for_each_skip_seq<F>(&self, data: &[u8], f: &mut F)
    where
        F: FnMut(&Marker),
    {
        let min_depth = self.config.min_depth;
        let mut marker = Marker::new(self.n_individuals);
        let mut pos = 0;

        while pos < data.len() {
            let line_end = memchr::memchr(b'\n', &data[pos..])
                .map(|p| pos + p)
                .unwrap_or(data.len());

            let line = strip_cr(&data[pos..line_end]);
            pos = line_end + 1;

            // The second tab starts the depth columns.
            let mut tab_iter = memchr::memchr_iter(b'\t', line);
            let _tab1 = tab_iter.next();
            let tab2 = match tab_iter.next() {
                Some(p) => p,
                None => continue,
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
        }
    }

    /// Full path: parse everything including id and sequence.
    #[inline]
    fn for_each_full<F>(&self, data: &[u8], f: &mut F)
    where
        F: FnMut(&Marker),
    {
        let min_depth = self.config.min_depth;
        let mut marker = Marker::new(self.n_individuals);
        let mut pos = 0;

        while pos < data.len() {
            let line_end = memchr::memchr(b'\n', &data[pos..])
                .map(|p| pos + p)
                .unwrap_or(data.len());

            let line = strip_cr(&data[pos..line_end]);
            pos = line_end + 1;

            let tabs: Vec<usize> = memchr::memchr_iter(b'\t', line).collect();
            if tabs.len() < 2 {
                continue;
            }

            marker.id.clear();
            marker
                .id
                .push_str(std::str::from_utf8(&line[0..tabs[0]]).unwrap_or(""));

            let seq_start = tabs[0] + 1;
            let seq_end = tabs[1];
            marker.sequence.clear();
            marker
                .sequence
                .push_str(std::str::from_utf8(&line[seq_start..seq_end]).unwrap_or(""));

            for_each_depth_field(line, tabs[1], |col, field| {
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
        }
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

        let target_chunk = 1 << 20;
        let chunks = self.make_line_aligned_chunks(data, target_chunk);

        chunks
            .into_par_iter()
            .try_for_each(|chunk| self.process_slice_serial(chunk, &f))
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

#[inline(always)]
fn strip_cr(line: &[u8]) -> &[u8] {
    line.strip_suffix(b"\r").unwrap_or(line)
}
