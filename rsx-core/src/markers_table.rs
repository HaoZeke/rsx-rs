// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Streaming markers table parser.
//!
//! Provides both inline (single-thread) and channel-based (two-thread) modes.
//! Inline mode is faster for analysis commands; channel mode available for
//! pipelines that benefit from decoupled parsing.

use crate::io::table_io::{TableHeader, fast_parse_u16};
use crate::marker::Marker;
use crate::popmap::Popmap;

use std::io::Read;
use std::path::Path;

const BUF_SIZE: usize = 65536;

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

/// Inline markers table iterator -- parses directly in the calling thread.
/// No channel overhead, no cloning. The marker is yielded by reference
/// via a callback to avoid allocation.
pub struct MarkersTableStream {
    pub header: TableHeader,
    groups: Vec<String>,
    group_names: Vec<String>,
    config: ParserConfig,
    path: std::path::PathBuf,
}

impl MarkersTableStream {
    /// Open a markers table and prepare for iteration.
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

        let group_names: Vec<String> = if let Some(pm) = popmap {
            pm.group_counts.keys().cloned().collect()
        } else {
            Vec::new()
        };

        Ok(MarkersTableStream {
            header,
            groups,
            group_names,
            config,
            path: path.to_path_buf(),
        })
    }

    /// Process all markers by calling `f` on each one.
    /// The marker is reused across calls (no allocation per marker).
    pub fn for_each<F>(&self, mut f: F) -> std::io::Result<()>
    where
        F: FnMut(&Marker),
    {
        let file = std::fs::File::open(&self.path)?;
        let mut reader = std::io::BufReader::with_capacity(BUF_SIZE, file);

        // Skip comment and header lines
        let mut line_buf = Vec::with_capacity(1024);
        loop {
            line_buf.clear();
            let n = read_line_bytes(&mut reader, &mut line_buf)?;
            if n == 0 {
                return Ok(());
            }
            if line_buf.starts_with(b"id\t") {
                line_buf.clear();
                break;
            }
            if !line_buf.starts_with(b"#") {
                break; // First data line
            }
        }

        let compute_groups = self.config.compute_groups && !self.groups.is_empty();
        let n_individuals = self.header.n_individuals;

        let mut marker = Marker::new(n_individuals);
        for gn in &self.group_names {
            marker.group_counts.insert(gn.clone(), 0);
        }

        let mut temp = Vec::with_capacity(256);
        let mut field_n: usize = 0;

        // Process any data line we already read
        if !line_buf.is_empty() {
            for &byte in &line_buf {
                process_byte(
                    byte, &mut marker, &mut temp, &mut field_n,
                    &self.config, compute_groups, &self.groups, &mut f,
                );
            }
            // Simulate newline if not present
            process_byte(
                b'\n', &mut marker, &mut temp, &mut field_n,
                &self.config, compute_groups, &self.groups, &mut f,
            );
        }

        let mut buffer = [0u8; BUF_SIZE];
        loop {
            let bytes_read = reader.read(&mut buffer)?;
            if bytes_read == 0 {
                break;
            }
            for &byte in &buffer[..bytes_read] {
                process_byte(
                    byte, &mut marker, &mut temp, &mut field_n,
                    &self.config, compute_groups, &self.groups, &mut f,
                );
            }
        }

        Ok(())
    }

    /// Collect all markers into a Vec. Use `for_each` when possible to avoid
    /// allocation.
    pub fn collect(&self) -> std::io::Result<Vec<Marker>> {
        let mut markers = Vec::new();
        self.for_each(|m| markers.push(m.clone()))?;
        Ok(markers)
    }

    /// Iterate over all markers (allocates -- prefer `for_each`).
    pub fn iter(&self) -> impl Iterator<Item = Marker> {
        self.collect().unwrap_or_default().into_iter()
    }
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn process_byte<F>(
    byte: u8,
    marker: &mut Marker,
    temp: &mut Vec<u8>,
    field_n: &mut usize,
    config: &ParserConfig,
    compute_groups: bool,
    groups: &[String],
    f: &mut F,
) where
    F: FnMut(&Marker),
{
    match byte {
        b'\t' => {
            handle_field(marker, temp, *field_n, config, compute_groups, groups);
            temp.clear();
            *field_n += 1;
        }
        b'\n' => {
            if *field_n >= 2 {
                handle_field(marker, temp, *field_n, config, compute_groups, groups);
            }
            temp.clear();
            *field_n = 0;

            f(marker);
            marker.reset(!config.store_sequence);
        }
        b'\r' => {}
        _ => {
            temp.push(byte);
        }
    }
}

#[inline(always)]
fn handle_field(
    marker: &mut Marker,
    temp: &[u8],
    field_n: usize,
    config: &ParserConfig,
    compute_groups: bool,
    groups: &[String],
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
                    if compute_groups
                        && field_n < groups.len()
                        && !groups[field_n].is_empty()
                    {
                        if let Some(count) = marker.group_counts.get_mut(&groups[field_n]) {
                            *count += 1;
                        }
                    }
                    marker.n_individuals += 1;
                }
            }
        }
    }
}

fn read_line_bytes<R: std::io::BufRead>(
    reader: &mut R,
    buf: &mut Vec<u8>,
) -> std::io::Result<usize> {
    let n = reader.read_until(b'\n', buf)?;
    if buf.ends_with(b"\n") {
        buf.pop();
    }
    if buf.ends_with(b"\r") {
        buf.pop();
    }
    Ok(n)
}
