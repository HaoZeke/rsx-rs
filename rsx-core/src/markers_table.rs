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

    /// Process all markers. Uses fast path when sequence isn't needed.
    pub fn for_each<F>(&self, mut f: F) -> std::io::Result<()>
    where
        F: FnMut(&Marker),
    {
        let data = &self.mmap[self.data_start..];
        if data.is_empty() {
            return Ok(());
        }

        if !self.config.store_sequence && self.config.min_depth <= 1 {
            // FAST PATH: skip id+seq, threshold=1 means just check != "0"
            self.for_each_fast_d1(data, &mut f);
        } else if !self.config.store_sequence {
            // MEDIUM PATH: skip id+seq, but need integer parse for depth
            self.for_each_skip_seq(data, &mut f);
        } else {
            // SLOW PATH: need sequence (for signif/subset FASTA output)
            self.for_each_full(data, &mut f);
        }

        Ok(())
    }

    /// Fast path: min_depth=1, no sequence needed.
    /// Skip id+seq fields. Check "0" vs non-"0" without integer parse.
    #[inline]
    fn for_each_fast_d1<F>(&self, data: &[u8], f: &mut F)
    where
        F: FnMut(&Marker),
    {
        let mut marker = Marker::new(self.n_individuals);
        let mut pos = 0;

        while pos < data.len() {
            // Find end of line
            let line_end = memchr::memchr(b'\n', &data[pos..])
                .map(|p| pos + p)
                .unwrap_or(data.len());

            let line = &data[pos..line_end];
            pos = line_end + 1;

            // Skip id field (find 1st tab)
            let tab1 = match memchr::memchr(b'\t', line) {
                Some(p) => p,
                None => continue,
            };
            // Skip sequence field (find 2nd tab)
            let tab2 = match memchr::memchr(b'\t', &line[tab1 + 1..]) {
                Some(p) => tab1 + 1 + p,
                None => continue,
            };

            // Parse depth fields: everything after tab2
            let mut col = 0usize;
            let mut field_start = tab2 + 1;

            while field_start < line.len() {
                // Find next tab or end of line
                let field_end = memchr::memchr(b'\t', &line[field_start..])
                    .map(|p| field_start + p)
                    .unwrap_or(line.len());

                let field = &line[field_start..field_end];

                // min_depth=1: present iff field != "0"
                let is_zero = field.len() == 1 && field[0] == b'0';
                if !is_zero && !field.is_empty() && col < self.n_individuals as usize {
                    marker.presence.set(col);
                    marker.n_individuals += 1;
                }

                col += 1;
                field_start = field_end + 1;
            }

            f(&marker);
            marker.reset(true); // keep_sequence=true (no sequence to clear)
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

            let line = &data[pos..line_end];
            pos = line_end + 1;

            // Skip id + sequence
            let tab1 = match memchr::memchr(b'\t', line) {
                Some(p) => p,
                None => continue,
            };
            let tab2 = match memchr::memchr(b'\t', &line[tab1 + 1..]) {
                Some(p) => tab1 + 1 + p,
                None => continue,
            };

            let mut col = 0usize;
            let mut field_start = tab2 + 1;

            while field_start < line.len() {
                let field_end = memchr::memchr(b'\t', &line[field_start..])
                    .map(|p| field_start + p)
                    .unwrap_or(line.len());

                let field = &line[field_start..field_end];
                let depth = fast_parse_u16(field);

                if col < self.n_individuals as usize {
                    marker.individual_depths[col] = depth;
                    if depth >= min_depth {
                        marker.presence.set(col);
                        marker.n_individuals += 1;
                    }
                }

                col += 1;
                field_start = field_end + 1;
            }

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
        let mut temp = Vec::with_capacity(256);
        let mut field_n: usize = 0;

        for &byte in data {
            match byte {
                b'\t' => {
                    handle_field(&mut marker, &temp, field_n, min_depth);
                    temp.clear();
                    field_n += 1;
                }
                b'\n' => {
                    if field_n >= 2 {
                        handle_field(&mut marker, &temp, field_n, min_depth);
                    }
                    temp.clear();
                    field_n = 0;
                    f(&marker);
                    marker.reset(false);
                }
                b'\r' => {}
                _ => {
                    temp.push(byte);
                }
            }
        }

        if field_n >= 2 {
            handle_field(&mut marker, &temp, field_n, min_depth);
            f(&marker);
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
}

#[inline(always)]
fn handle_field(marker: &mut Marker, temp: &[u8], field_n: usize, min_depth: u16) {
    match field_n {
        0 => {
            marker.id.clear();
            marker.id.push_str(std::str::from_utf8(temp).unwrap_or(""));
        }
        1 => {
            marker.sequence.clear();
            marker.sequence.push_str(std::str::from_utf8(temp).unwrap_or(""));
        }
        _ => {
            let depth = fast_parse_u16(temp);
            let idx = field_n - 2;
            if idx < marker.individual_depths.len() {
                marker.individual_depths[idx] = depth;
                if depth >= min_depth {
                    marker.presence.set(idx);
                    marker.n_individuals += 1;
                }
            }
        }
    }
}
