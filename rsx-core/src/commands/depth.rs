// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! `depth` command: compute retained read statistics per individual.
//!
//! Two modes:
//! - Default: exact median via in-memory accumulation
//! - Streaming (--streaming): exact median via external sort of
//!   (individual, depth) pairs. O(buffer_size) memory.

use crate::markers_table::{MarkersTableStream, ParserConfig};
use crate::popmap::Popmap;
use crate::stats;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::io::{BufWriter, Read, Write};
use std::path::Path;

pub struct DepthParams {
    pub markers_table_path: String,
    pub popmap_file_path: String,
    pub output_file_path: String,
    pub min_frequency: f32,
    pub streaming: bool,
}

pub fn run(params: &DepthParams) -> Result<(), Box<dyn std::error::Error>> {
    if params.streaming {
        run_streaming(params)
    } else {
        run_exact(params)
    }
}

/// Exact mode: accumulates depth vectors for exact median.
/// Works for tables that fit in RAM.
fn run_exact(params: &DepthParams) -> Result<(), Box<dyn std::error::Error>> {
    let table_path = Path::new(&params.markers_table_path);
    let popmap = Popmap::from_file(Path::new(&params.popmap_file_path))?;

    let config = ParserConfig {
        store_sequence: false,
        store_depths: true,
        compute_groups: false,
        min_depth: 1,
    };

    let stream = MarkersTableStream::open(table_path, Some(&popmap), config)?;
    let n_individuals = stream.header.n_individuals as usize;
    let min_individuals = (params.min_frequency * stream.header.n_individuals as f32) as u32;

    let mut depths: Vec<Vec<u16>> = vec![Vec::new(); n_individuals];
    let mut individual_markers_count: Vec<u64> = vec![0; n_individuals];
    let mut individual_reads_count: Vec<u64> = vec![0; n_individuals];

    stream.for_each(|marker| {
        for i in 0..n_individuals {
            let d = marker.individual_depths[i];
            if marker.n_individuals >= min_individuals {
                depths[i].push(d);
            }
            if d > 0 {
                individual_markers_count[i] += 1;
                individual_reads_count[i] += d as u64;
            }
        }
    })?;

    if depths.iter().any(|d| d.is_empty()) {
        return Err(format!(
            "No markers were present in at least {}% of all individuals ({}/{} individuals)",
            (params.min_frequency * 100.0) as u32,
            min_individuals,
            n_individuals
        )
        .into());
    }

    let header_cols = &stream.header.columns;
    let mut output = std::fs::File::create(&params.output_file_path)?;
    writeln!(
        output,
        "Sample\tGroup\tReads\tMarkers\tRetained\tMin_depth\tMax_depth\tMedian_depth\tAverage_depth"
    )?;

    for i in 0..n_individuals {
        let individual_name = &header_cols[i + 2];
        let group = popmap.get_group(individual_name).unwrap_or("");
        depths[i].sort_unstable();
        let size = depths[i].len() as u64;
        let min_d = depths[i][0];
        let max_d = *depths[i].last().unwrap();
        let total: u64 = depths[i].iter().map(|&d| d as u64).sum();
        let median_d = stats::find_median(&mut depths[i]);
        let avg_d = total / size;
        writeln!(
            output,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            individual_name,
            group,
            individual_reads_count[i],
            individual_markers_count[i],
            size,
            min_d,
            max_d,
            median_d,
            avg_d
        )?;
    }
    Ok(())
}

/// Streaming mode: exact median via external sort of (individual_idx, depth) pairs.
/// Memory: O(buffer_size), not O(n_markers * n_individuals).
fn run_streaming(params: &DepthParams) -> Result<(), Box<dyn std::error::Error>> {
    let table_path = Path::new(&params.markers_table_path);
    let popmap = Popmap::from_file(Path::new(&params.popmap_file_path))?;

    let config = ParserConfig {
        store_sequence: false,
        store_depths: true,
        compute_groups: false,
        min_depth: 1,
    };

    let stream = MarkersTableStream::open(table_path, Some(&popmap), config)?;
    let n_individuals = stream.header.n_individuals as usize;
    let min_individuals = (params.min_frequency * stream.header.n_individuals as f32) as u32;

    // Online accumulators: O(n_individuals) memory
    let mut ind_markers: Vec<u64> = vec![0; n_individuals];
    let mut ind_reads: Vec<u64> = vec![0; n_individuals];
    let mut ind_total: Vec<u64> = vec![0; n_individuals]; // total retained (including zeros)
    let mut ind_nonzero: Vec<u64> = vec![0; n_individuals]; // non-zero retained count
    let mut ind_sum: Vec<u64> = vec![0; n_individuals];
    let mut ind_min: Vec<u16> = vec![u16::MAX; n_individuals];
    let mut ind_max: Vec<u16> = vec![0; n_individuals];

    // SPARSE external sort: only buffer NON-ZERO (individual, depth) pairs.
    // Zeros tracked by count (ind_total - ind_nonzero) for exact median.
    // This reduces sort size by ~70% for typical RAD-seq data.
    const BUFFER_ENTRIES: usize = 50_000_000; // 50M entries * 4 bytes = 200MB
    let mut buffer: Vec<(u16, u16)> = Vec::with_capacity(BUFFER_ENTRIES);
    let temp_dir = tempfile::TempDir::new()?;
    let mut chunk_paths: Vec<std::path::PathBuf> = Vec::new();

    log::info!("depth streaming: reading + sparse sorting (zeros skipped)");

    let mut flush_err: Option<Box<dyn std::error::Error>> = None;

    stream.for_each(|marker| {
        if flush_err.is_some() {
            return;
        }
        let retained = marker.n_individuals >= min_individuals;
        for i in 0..n_individuals {
            let d = marker.individual_depths[i];
            if d > 0 {
                ind_markers[i] += 1;
                ind_reads[i] += d as u64;
            }
            if retained {
                ind_total[i] += 1;
                ind_sum[i] += d as u64;
                if d > 0 {
                    // Only sort non-zero depths (sparse optimization)
                    ind_nonzero[i] += 1;
                    if d < ind_min[i] {
                        ind_min[i] = d;
                    }
                    if d > ind_max[i] {
                        ind_max[i] = d;
                    }

                    buffer.push((i as u16, d));
                    if buffer.len() >= BUFFER_ENTRIES {
                        match flush_depth_chunk(&mut buffer, &temp_dir, chunk_paths.len()) {
                            Ok(p) => chunk_paths.push(p),
                            Err(e) => {
                                flush_err = Some(e);
                                return;
                            }
                        }
                    }
                }
                // d == 0: tracked by ind_total - ind_nonzero, not stored
            }
        }
    })?;

    if let Some(e) = flush_err {
        return Err(e);
    }

    if !buffer.is_empty() {
        let p = flush_depth_chunk(&mut buffer, &temp_dir, chunk_paths.len())?;
        chunk_paths.push(p);
    }
    drop(buffer);

    log::info!("depth streaming: {} chunks written", chunk_paths.len());

    if ind_total.contains(&0) {
        return Err(format!(
            "No markers were present in at least {}% of all individuals ({}/{} individuals)",
            (params.min_frequency * 100.0) as u32,
            min_individuals,
            n_individuals
        )
        .into());
    }

    // K-way merge the sorted (individual_idx, depth) pairs
    // Scan: for each individual, depths are contiguous and sorted -> pick median
    log::info!("depth streaming: merging for exact median");

    let mut readers: Vec<DepthChunkReader> = chunk_paths
        .iter()
        .map(|p| DepthChunkReader::open(p))
        .collect::<Result<Vec<_>, _>>()?;

    let mut heap: BinaryHeap<DepthHeapEntry> = BinaryHeap::new();
    for (idx, r) in readers.iter_mut().enumerate() {
        if let Some((ind, dep)) = r.next_pair()? {
            heap.push(DepthHeapEntry {
                ind,
                dep,
                chunk: idx,
            });
        }
    }

    // Compute median accounting for implicit zeros.
    // The full sorted sequence for individual i is:
    //   [0, 0, ..., 0 (n_zeros times), sorted_nonzero_1, sorted_nonzero_2, ...]
    // The median position is ind_total[i] / 2.
    // If median_pos < n_zeros, median = 0 (already initialized).
    // Otherwise, we need position (median_pos - n_zeros) in the sorted non-zeros.
    let mut medians: Vec<u16> = vec![0; n_individuals];
    let mut nonzero_pos: Vec<u64> = vec![0; n_individuals]; // position within sorted non-zeros
    let median_targets: Vec<i64> = (0..n_individuals)
        .map(|i| {
            let n_zeros = ind_total[i] - ind_nonzero[i];
            let median_pos = ind_total[i] / 2;
            if median_pos < n_zeros {
                -1 // median is 0, already set
            } else {
                (median_pos - n_zeros) as i64 // target position in non-zero stream
            }
        })
        .collect();

    while let Some(top) = heap.pop() {
        let i = top.ind as usize;
        if median_targets[i] >= 0 && nonzero_pos[i] == median_targets[i] as u64 {
            medians[i] = top.dep;
        }
        nonzero_pos[i] += 1;

        if let Some((ind, dep)) = readers[top.chunk].next_pair()? {
            heap.push(DepthHeapEntry {
                ind,
                dep,
                chunk: top.chunk,
            });
        }
    }

    // Write output
    let header_cols = &stream.header.columns;
    let mut output = std::fs::File::create(&params.output_file_path)?;
    writeln!(
        output,
        "Sample\tGroup\tReads\tMarkers\tRetained\tMin_depth\tMax_depth\tMedian_depth\tAverage_depth"
    )?;

    for i in 0..n_individuals {
        let individual_name = &header_cols[i + 2];
        let group = popmap.get_group(individual_name).unwrap_or("");
        let avg = ind_sum[i].checked_div(ind_total[i]).unwrap_or(0);
        // All-zero retained depths have zero as their minimum.
        let min_d = if ind_nonzero[i] == 0 { 0 } else { ind_min[i] };
        // If any zeros exist and min_nonzero > 0, the true min is 0
        let min_d = if ind_total[i] > ind_nonzero[i] {
            0
        } else {
            min_d
        };

        writeln!(
            output,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            individual_name,
            group,
            ind_reads[i],
            ind_markers[i],
            ind_total[i],
            min_d,
            ind_max[i],
            medians[i],
            avg
        )?;
    }
    Ok(())
}

// --- Depth external sort helpers ---

fn flush_depth_chunk(
    buffer: &mut Vec<(u16, u16)>,
    temp_dir: &tempfile::TempDir,
    chunk_idx: usize,
) -> Result<std::path::PathBuf, Box<dyn std::error::Error>> {
    // Sort by (individual_idx, depth) for contiguous per-individual sorted runs
    buffer.sort_unstable();
    let path = temp_dir.path().join(format!("depth_{:04}.lz4", chunk_idx));
    let file = std::fs::File::create(&path)?;
    let mut enc = lz4_flex::frame::FrameEncoder::new(BufWriter::new(file));
    for &(ind, dep) in buffer.iter() {
        enc.write_all(&ind.to_le_bytes())?;
        enc.write_all(&dep.to_le_bytes())?;
    }
    enc.finish()?;
    buffer.clear();
    Ok(path)
}

struct DepthChunkReader {
    reader: lz4_flex::frame::FrameDecoder<std::io::BufReader<std::fs::File>>,
}

impl DepthChunkReader {
    fn open(path: &Path) -> std::io::Result<Self> {
        let file = std::fs::File::open(path)?;
        Ok(DepthChunkReader {
            reader: lz4_flex::frame::FrameDecoder::new(std::io::BufReader::new(file)),
        })
    }

    fn next_pair(&mut self) -> std::io::Result<Option<(u16, u16)>> {
        let mut buf = [0u8; 4];
        match self.reader.read_exact(&mut buf) {
            Ok(()) => {
                let ind = u16::from_le_bytes([buf[0], buf[1]]);
                let dep = u16::from_le_bytes([buf[2], buf[3]]);
                Ok(Some((ind, dep)))
            }
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => Ok(None),
            Err(e) => Err(e),
        }
    }
}

struct DepthHeapEntry {
    ind: u16,
    dep: u16,
    chunk: usize,
}

impl PartialEq for DepthHeapEntry {
    fn eq(&self, other: &Self) -> bool {
        self.ind == other.ind && self.dep == other.dep
    }
}
impl Eq for DepthHeapEntry {}
impl PartialOrd for DepthHeapEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for DepthHeapEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        // Min-heap: reversed
        other.ind.cmp(&self.ind).then(other.dep.cmp(&self.dep))
    }
}
