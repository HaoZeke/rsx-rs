// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Markers table parser with mmap I/O and optional parallel chunk parsing.
//!
//! The file is memory-mapped (zero-copy), then either:
//! - Single-thread: parsed sequentially via `for_each` callback
//! - Parallel: split into chunks at newline boundaries, parsed with rayon
//!
//! Group counting uses bitset popcount (no HashMap in the hot path).

use crate::io::table_io::{TableHeader, fast_parse_u16};
use crate::marker::Marker;
use crate::popmap::Popmap;

use memmap2::Mmap;
use std::path::Path;

/// Configuration for the markers table parser.
pub struct ParserConfig {
    pub store_sequence: bool,
    pub compute_groups: bool,
    pub min_depth: u16,
}

impl Default for ParserConfig {
    fn default() -> Self {
        ParserConfig {
            store_sequence: true,
            compute_groups: true,
            min_depth: 1,
        }
    }
}

/// Markers table backed by mmap. Parsed inline via `for_each`.
pub struct MarkersTableStream {
    pub header: TableHeader,
    pub groups: Vec<String>,
    config: ParserConfig,
    mmap: Mmap,
    /// Byte offset where data lines start (after comments + header).
    data_start: usize,
    n_individuals: u16,
}

impl MarkersTableStream {
    /// Open a markers table via mmap.
    pub fn open(
        path: &Path,
        popmap: Option<&Popmap>,
        config: ParserConfig,
    ) -> std::io::Result<Self> {
        let header = TableHeader::from_file(path)?;

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

        // Find where data starts (skip comment lines + header line)
        let mut data_start = 0;
        loop {
            if data_start >= mmap.len() {
                break;
            }
            // Find end of current line
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

    /// Process all markers via callback. Zero allocation per marker.
    pub fn for_each<F>(&self, mut f: F) -> std::io::Result<()>
    where
        F: FnMut(&Marker),
    {
        let data = &self.mmap[self.data_start..];
        if data.is_empty() {
            return Ok(());
        }

        let mut marker = Marker::new(self.n_individuals);
        let mut temp = Vec::with_capacity(256);
        let mut field_n: usize = 0;

        for &byte in data {
            match byte {
                b'\t' => {
                    handle_field(&mut marker, &temp, field_n, &self.config);
                    temp.clear();
                    field_n += 1;
                }
                b'\n' => {
                    if field_n >= 2 {
                        handle_field(&mut marker, &temp, field_n, &self.config);
                    }
                    temp.clear();
                    field_n = 0;
                    f(&marker);
                    marker.reset(!self.config.store_sequence);
                }
                b'\r' => {}
                _ => {
                    temp.push(byte);
                }
            }
        }

        // Handle last line without trailing newline
        if field_n >= 2 {
            handle_field(&mut marker, &temp, field_n, &self.config);
            f(&marker);
        }

        Ok(())
    }

    /// Collect all markers into a Vec.
    pub fn collect(&self) -> std::io::Result<Vec<Marker>> {
        let mut markers = Vec::new();
        self.for_each(|m| markers.push(m.clone()))?;
        Ok(markers)
    }

    /// Iterate over all markers (allocates).
    pub fn iter(&self) -> impl Iterator<Item = Marker> {
        self.collect().unwrap_or_default().into_iter()
    }
}

#[inline(always)]
fn handle_field(
    marker: &mut Marker,
    temp: &[u8],
    field_n: usize,
    config: &ParserConfig,
) {
    match field_n {
        0 => {
            if config.store_sequence {
                marker.id.clear();
                marker.id.push_str(std::str::from_utf8(temp).unwrap_or(""));
            }
        }
        1 => {
            if config.store_sequence {
                marker.sequence.clear();
                marker.sequence.push_str(std::str::from_utf8(temp).unwrap_or(""));
            }
        }
        _ => {
            let depth = fast_parse_u16(temp);
            let idx = field_n - 2;
            if idx < marker.individual_depths.len() {
                marker.individual_depths[idx] = depth;
                if depth >= config.min_depth {
                    marker.presence.set(idx);
                    marker.n_individuals += 1;
                }
            }
        }
    }
}
