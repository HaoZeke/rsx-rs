// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! `depth` command: compute retained read statistics per individual.
//!
//! Two modes:
//! - Default: exact median via in-memory accumulation (original behavior)
//! - Streaming (--streaming): exact median via external sort of
//!   (individual, depth) pairs. O(buffer_size) memory.

use crate::markers_table::{MarkersTableStream, ParserConfig};
use crate::popmap::Popmap;
use crate::stats;
use std::collections::BinaryHeap;
use std::cmp::Ordering;
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

/// Exact mode: accumulates depth vectors for exact median (original behavior).
/// Works for tables that fit in RAM.
fn run_exact(params: &DepthParams) -> Result<(), Box<dyn std::error::Error>> {
    let table_path = Path::new(&params.markers_table_path);
    let popmap = Popmap::from_file(Path::new(&params.popmap_file_path))?;

    let config = ParserConfig {
        store_sequence: false,
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
            (params.min_frequency * 100.0) as u32, min_individuals, n_individuals
        ).into());
    }

    let header_cols = &stream.header.columns;
    let mut output = std::fs::File::create(&params.output_file_path)?;
    writeln!(output,
        "Sample\tGroup\tReads\tMarkers\tRetained\tMin_depth\tMax_depth\tMedian_depth\tAverage_depth")?;

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
        writeln!(output, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            individual_name, group, individual_reads_count[i],
            individual_markers_count[i], size, min_d, max_d, median_d, avg_d)?;
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
        compute_groups: false,
        min_depth: 1,
    };

    let stream = MarkersTableStream::open(table_path, Some(&popmap), config)?;
    let n_individuals = stream.header.n_individuals as usize;
    let min_individuals = (params.min_frequency * stream.header.n_individuals as f32) as u32;

    // Online accumulators for reads/markers/min/max/sum/count (O(n_ind))
    let mut ind_markers: Vec<u64> = vec![0; n_individuals];
    let mut ind_reads: Vec<u64> = vec![0; n_individuals];
    let mut ind_count: Vec<u64> = vec![0; n_individuals];
    let mut ind_sum: Vec<u64> = vec![0; n_individuals];
    let mut ind_min: Vec<u16> = vec![u16::MAX; n_individuals];
    let mut ind_max: Vec<u16> = vec![0; n_individuals];

    // Buffer (individual_idx: u16, depth: u16) pairs for external sort
    const BUFFER_ENTRIES: usize = 50_000_000; // 50M entries * 4 bytes = 200MB
    let mut buffer: Vec<(u16, u16)> = Vec::with_capacity(BUFFER_ENTRIES);
    let temp_dir = tempfile::TempDir::new()?;
    let mut chunk_paths: Vec<std::path::PathBuf> = Vec::new();

    log::info!("depth streaming: pass 1 - reading + sorting");

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
                ind_count[i] += 1;
                ind_sum[i] += d as u64;
                if d < ind_min[i] { ind_min[i] = d; }
                if d > ind_max[i] { ind_max[i] = d; }

                buffer.push((i as u16, d));
                if buffer.len() >= BUFFER_ENTRIES {
                    match flush_depth_chunk(&mut buffer, &temp_dir, chunk_paths.len()) {
                        Ok(p) => chunk_paths.push(p),
                        Err(e) => { flush_err = Some(e); return; }
                    }
                }
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

    if ind_count.contains(&0) {
        return Err(format!(
            "No markers were present in at least {}% of all individuals ({}/{} individuals)",
            (params.min_frequency * 100.0) as u32, min_individuals, n_individuals
        ).into());
    }

    // K-way merge the sorted (individual_idx, depth) pairs
    // Scan: for each individual, depths are contiguous and sorted -> pick median
    log::info!("depth streaming: merging for exact median");

    let mut readers: Vec<DepthChunkReader> = chunk_paths.iter()
        .map(|p| DepthChunkReader::open(p))
        .collect::<Result<Vec<_>, _>>()?;

    let mut heap: BinaryHeap<DepthHeapEntry> = BinaryHeap::new();
    for (idx, r) in readers.iter_mut().enumerate() {
        if let Some((ind, dep)) = r.next_pair()? {
            heap.push(DepthHeapEntry { ind, dep, chunk: idx });
        }
    }

    let mut medians: Vec<u16> = vec![0; n_individuals];
    let mut current_ind: u16 = 0;
    let mut current_pos: u64 = 0;

    while let Some(top) = heap.pop() {
        if top.ind != current_ind {
            current_ind = top.ind;
            current_pos = 0;
        }

        let target = ind_count[current_ind as usize] / 2;
        if current_pos == target {
            medians[current_ind as usize] = top.dep;
        }
        current_pos += 1;

        if let Some((ind, dep)) = readers[top.chunk].next_pair()? {
            heap.push(DepthHeapEntry { ind, dep, chunk: top.chunk });
        }
    }

    // Write output
    let header_cols = &stream.header.columns;
    let mut output = std::fs::File::create(&params.output_file_path)?;
    writeln!(output,
        "Sample\tGroup\tReads\tMarkers\tRetained\tMin_depth\tMax_depth\tMedian_depth\tAverage_depth")?;

    for i in 0..n_individuals {
        let individual_name = &header_cols[i + 2];
        let group = popmap.get_group(individual_name).unwrap_or("");
        let avg = ind_sum[i].checked_div(ind_count[i]).unwrap_or(0);
        writeln!(output, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            individual_name, group, ind_reads[i], ind_markers[i],
            ind_count[i], ind_min[i], ind_max[i], medians[i], avg)?;
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
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> { Some(self.cmp(other)) }
}
impl Ord for DepthHeapEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        // Min-heap: reversed
        other.ind.cmp(&self.ind).then(other.dep.cmp(&self.dep))
    }
}
