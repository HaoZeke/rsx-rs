// GPL-3.0-or-later
// Copyright 2024--present radsex-rs developers

//! Streaming markers table parser using crossbeam channels.
//!
//! Replaces the C++ mutex+queue+busy-wait pattern with a bounded channel.
//! A producer thread reads the TSV file and sends batches of `Marker` objects
//! through the channel; consumers receive them without polling.

use crate::io::table_io::{TableHeader, fast_parse_u16};
use crate::marker::Marker;
use crate::popmap::Popmap;

use crossbeam::channel::{self, Receiver, Sender};
use std::io::Read;
use std::path::Path;
use std::thread;

const BATCH_SIZE: usize = 100;
const CHANNEL_CAPACITY: usize = 10_000;
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

/// Mutable parsing state passed through the byte-processing loop.
struct ParseState<'a> {
    marker: Marker,
    temp: Vec<u8>,
    field_n: usize,
    batch: Vec<Marker>,
    config: &'a ParserConfig,
    compute_groups: bool,
    groups: &'a [String],
    group_names: &'a [String],
    tx: &'a Sender<Vec<Marker>>,
}

impl<'a> ParseState<'a> {
    fn process_bytes(&mut self, data: &[u8]) {
        for &byte in data {
            match byte {
                b'\t' => {
                    self.handle_field();
                    self.temp.clear();
                    self.field_n += 1;
                }
                b'\n' => {
                    if self.field_n >= 2 {
                        self.handle_field();
                    }
                    self.temp.clear();
                    self.field_n = 0;

                    self.batch.push(self.marker.clone());

                    if self.batch.len() >= BATCH_SIZE {
                        let full_batch =
                            std::mem::replace(&mut self.batch, Vec::with_capacity(BATCH_SIZE));
                        let _ = self.tx.send(full_batch);
                    }

                    self.marker.reset(!self.config.store_sequence);
                    for gn in self.group_names {
                        self.marker.group_counts.insert(gn.clone(), 0);
                    }
                }
                b'\r' => {}
                _ => {
                    self.temp.push(byte);
                }
            }
        }
    }

    fn handle_field(&mut self) {
        match self.field_n {
            0 => {
                if self.config.store_sequence {
                    self.marker.id = String::from_utf8_lossy(&self.temp).into_owned();
                }
            }
            1 => {
                if self.config.store_sequence {
                    self.marker.sequence = String::from_utf8_lossy(&self.temp).into_owned();
                }
            }
            _ => {
                let depth = fast_parse_u16(&self.temp);
                let idx = self.field_n - 2;
                if idx < self.marker.individual_depths.len() {
                    self.marker.individual_depths[idx] = depth;
                    if depth >= self.config.min_depth {
                        if self.compute_groups
                            && self.field_n < self.groups.len()
                            && !self.groups[self.field_n].is_empty()
                        {
                            if let Some(count) =
                                self.marker.group_counts.get_mut(&self.groups[self.field_n])
                            {
                                *count += 1;
                            }
                        }
                        self.marker.n_individuals += 1;
                    }
                }
            }
        }
    }

    fn flush(self) {
        if !self.batch.is_empty() {
            let _ = self.tx.send(self.batch);
        }
    }
}

/// A streaming markers table that produces markers via a channel.
pub struct MarkersTableStream {
    pub header: TableHeader,
    receiver: Receiver<Vec<Marker>>,
    _handle: thread::JoinHandle<()>,
}

impl MarkersTableStream {
    /// Open a markers table file and start the parser thread.
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

        let n_individuals = header.n_individuals;
        let path = path.to_path_buf();

        let (tx, rx): (Sender<Vec<Marker>>, Receiver<Vec<Marker>>) =
            channel::bounded(CHANNEL_CAPACITY / BATCH_SIZE);

        let handle = thread::spawn(move || {
            if let Err(e) = parse_table(&path, n_individuals, &groups, &group_names, &config, &tx)
            {
                log::error!("Table parser error: {e}");
            }
        });

        Ok(MarkersTableStream {
            header,
            receiver: rx,
            _handle: handle,
        })
    }

    /// Receive the next batch of markers. Returns None when parsing is done.
    pub fn next_batch(&self) -> Option<Vec<Marker>> {
        self.receiver.recv().ok()
    }

    /// Iterate over all markers (consuming batches).
    pub fn iter(&self) -> impl Iterator<Item = Marker> + '_ {
        std::iter::from_fn(move || self.next_batch()).flat_map(|batch| batch.into_iter())
    }
}

fn parse_table(
    path: &Path,
    n_individuals: u16,
    groups: &[String],
    group_names: &[String],
    config: &ParserConfig,
    tx: &Sender<Vec<Marker>>,
) -> std::io::Result<()> {
    let file = std::fs::File::open(path)?;
    let mut reader = std::io::BufReader::with_capacity(BUF_SIZE, file);

    // Skip comment and header lines
    let mut line_buf = Vec::with_capacity(1024);
    loop {
        line_buf.clear();
        let n = read_line_bytes(&mut reader, &mut line_buf)?;
        if n == 0 {
            return Ok(());
        }
        if !line_buf.starts_with(b"#") && !line_buf.starts_with(b"id\t") {
            break;
        }
        if line_buf.starts_with(b"id\t") {
            line_buf.clear();
            break;
        }
    }

    let compute_groups = config.compute_groups && !groups.is_empty();

    let mut marker = Marker::new(n_individuals);
    for gn in group_names {
        marker.group_counts.insert(gn.clone(), 0);
    }

    let mut state = ParseState {
        marker,
        temp: Vec::with_capacity(256),
        field_n: 0,
        batch: Vec::with_capacity(BATCH_SIZE),
        config,
        compute_groups,
        groups,
        group_names,
        tx,
    };

    if !line_buf.is_empty() {
        state.process_bytes(&line_buf);
        line_buf.clear();
    }

    let mut buffer = [0u8; BUF_SIZE];
    loop {
        let bytes_read = reader.read(&mut buffer)?;
        if bytes_read == 0 {
            break;
        }
        state.process_bytes(&buffer[..bytes_read]);
    }

    state.flush();
    Ok(())
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
