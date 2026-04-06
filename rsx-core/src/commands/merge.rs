// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! External sort-merge for marker depth tables.
//!
//! Bounded-memory merge of multiple marker tables using chunked external sort:
//! 1. Read all input files, buffer entries in memory (configurable limit)
//! 2. When buffer full: sort by packed sequence, write lz4-compressed temp file
//! 3. K-way merge from sorted temp files, coalesce equal sequences
//! 4. Write merged TSV output
//!
//! Memory usage: ~500MB regardless of dataset size (75M+ sequences supported).

use crate::io::seq_reader::{pack_2bit, unpack_2bit};
use crate::io::table_io::fast_parse_u16;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

/// Default: buffer 2M entries before flushing to disk (~400MB)
const DEFAULT_BUFFER_SIZE: usize = 2_000_000;

pub struct MergeParams {
    pub input_files: Vec<String>,
    pub output_file_path: String,
    pub buffer_size: Option<usize>,
}

/// A single merge entry: one input line's contribution.
#[derive(Clone)]
struct MergeEntry {
    packed_seq: Vec<u8>,
    sample_offset: u16,
    depths: Vec<u16>,
}

/// Wrapper for BinaryHeap (min-heap via Reverse ordering).
struct HeapEntry {
    entry: MergeEntry,
    chunk_idx: usize,
}

impl PartialEq for HeapEntry {
    fn eq(&self, other: &Self) -> bool {
        self.entry.packed_seq == other.entry.packed_seq
    }
}

impl Eq for HeapEntry {}

impl PartialOrd for HeapEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for HeapEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        // Reversed for min-heap
        other.entry.packed_seq.cmp(&self.entry.packed_seq)
    }
}

pub fn run(params: &MergeParams) -> Result<(), Box<dyn std::error::Error>> {
    let n_inputs = params.input_files.len();
    let buffer_limit = params.buffer_size.unwrap_or(DEFAULT_BUFFER_SIZE);
    eprintln!(
        "Merging {} marker tables (buffer: {} entries)...",
        n_inputs, buffer_limit
    );

    // Phase 1: read headers
    let mut all_samples: Vec<String> = Vec::new();
    let mut file_samples: Vec<Vec<String>> = Vec::new();
    let mut sample_offsets: Vec<u16> = Vec::new();

    for (i, path) in params.input_files.iter().enumerate() {
        let file = std::fs::File::open(path)?;
        let reader = BufReader::new(file);
        let mut lines = reader.lines();

        let first_line = lines.next().ok_or("Empty file")??;
        let n_markers: u64 = first_line
            .split(':')
            .nth(1)
            .map(|s| s.trim().parse().unwrap_or(0))
            .unwrap_or(0);

        let header_line = lines.next().ok_or("Missing header")??;
        let parts: Vec<&str> = header_line.split('\t').collect();
        let samples: Vec<String> = parts[2..].iter().map(|s| s.to_string()).collect();

        eprintln!(
            "  Input {}: {} ({} samples, {} markers)",
            i + 1,
            Path::new(path)
                .file_name()
                .unwrap_or_default()
                .to_string_lossy(),
            samples.len(),
            n_markers
        );

        sample_offsets.push(all_samples.len() as u16);
        all_samples.extend(samples.clone());
        file_samples.push(samples);
    }

    let total_samples = all_samples.len();
    eprintln!("  Total samples: {}", total_samples);

    // Phase 2: read all files, buffer entries, sort + flush to temp files
    let mut buffer: Vec<MergeEntry> = Vec::with_capacity(buffer_limit);
    let temp_dir = tempfile::TempDir::new()?;
    let mut chunk_paths: Vec<std::path::PathBuf> = Vec::new();
    let mut total_lines = 0u64;

    for (i, path) in params.input_files.iter().enumerate() {
        let n_samples_this = file_samples[i].len();
        let offset = sample_offsets[i];
        eprintln!("  Reading input {}...", i + 1);

        let file = std::fs::File::open(path)?;
        let reader = BufReader::with_capacity(1 << 20, file);

        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') || line.starts_with("id\t") {
                continue;
            }

            let mut parts = line.splitn(3, '\t');
            let _id = parts.next();
            let seq_str = match parts.next() {
                Some(s) => s,
                None => continue,
            };
            let depths_str = parts.next().unwrap_or("");

            let depths: Vec<u16> = depths_str
                .split('\t')
                .take(n_samples_this)
                .map(|s| fast_parse_u16(s.as_bytes()))
                .collect();

            buffer.push(MergeEntry {
                packed_seq: pack_2bit(seq_str.as_bytes()),
                sample_offset: offset,
                depths,
            });

            total_lines += 1;

            if buffer.len() >= buffer_limit {
                let chunk_path = flush_chunk(&mut buffer, &temp_dir, chunk_paths.len())?;
                chunk_paths.push(chunk_path);
                eprintln!(
                    "    Flushed chunk {} ({} total lines so far)",
                    chunk_paths.len(),
                    total_lines
                );
            }
        }
    }

    // Flush remaining buffer
    if !buffer.is_empty() {
        let n = buffer.len();
        let chunk_path = flush_chunk(&mut buffer, &temp_dir, chunk_paths.len())?;
        chunk_paths.push(chunk_path);
        eprintln!(
            "    Flushed final chunk {} ({} entries)",
            chunk_paths.len(),
            n
        );
    }
    drop(buffer);

    eprintln!(
        "  Sort phase done: {} total lines, {} chunks",
        total_lines,
        chunk_paths.len()
    );

    // Phase 3: k-way merge
    let mut chunk_readers: Vec<ChunkReader> = chunk_paths
        .iter()
        .map(|p| ChunkReader::open(p))
        .collect::<Result<Vec<_>, _>>()?;

    let mut heap = BinaryHeap::new();
    for (idx, reader) in chunk_readers.iter_mut().enumerate() {
        if let Some(entry) = reader.next_entry()? {
            heap.push(HeapEntry {
                entry,
                chunk_idx: idx,
            });
        }
    }

    eprintln!("  Writing {}...", params.output_file_path);
    let out_file = std::fs::File::create(&params.output_file_path)?;
    let mut writer = BufWriter::with_capacity(1 << 20, out_file);

    // Write placeholder header (fix count after)
    writeln!(writer, "#Number of markers : {:<20}", 0)?;
    write!(writer, "id\tsequence")?;
    for s in &all_samples {
        write!(writer, "\t{}", s)?;
    }
    writeln!(writer)?;

    let mut current_seq: Option<Vec<u8>> = None;
    let mut current_depths = vec![0u16; total_samples];
    let mut n_merged = 0u64;

    while let Some(top) = heap.pop() {
        let entry = top.entry;
        let chunk_idx = top.chunk_idx;

        let is_new = match &current_seq {
            Some(seq) => *seq != entry.packed_seq,
            None => true,
        };

        if is_new {
            if let Some(ref seq) = current_seq {
                write_merged_row(&mut writer, n_merged, seq, &current_depths)?;
                n_merged += 1;
                current_depths.fill(0);

                if n_merged % 5_000_000 == 0 {
                    eprintln!("    Written {} merged markers", n_merged);
                }
            }
            current_seq = Some(entry.packed_seq.clone());
        }

        // Merge depths
        let offset = entry.sample_offset as usize;
        for (j, &d) in entry.depths.iter().enumerate() {
            let col = offset + j;
            if col < total_samples {
                if current_depths[col] == 0 {
                    current_depths[col] = d;
                } else {
                    current_depths[col] = current_depths[col].saturating_add(d);
                }
            }
        }

        if let Some(next) = chunk_readers[chunk_idx].next_entry()? {
            heap.push(HeapEntry {
                entry: next,
                chunk_idx,
            });
        }
    }

    // Flush final
    if let Some(ref seq) = current_seq {
        write_merged_row(&mut writer, n_merged, seq, &current_depths)?;
        n_merged += 1;
    }

    writer.flush()?;
    drop(writer);

    fix_header_count(&params.output_file_path, n_merged)?;

    eprintln!(
        "  Done: {} merged markers x {} samples -> {}",
        n_merged, total_samples, params.output_file_path
    );
    Ok(())
}

fn write_merged_row(
    writer: &mut impl Write,
    id: u64,
    packed_seq: &[u8],
    depths: &[u16],
) -> std::io::Result<()> {
    let unpacked = unpack_2bit(packed_seq);
    let seq_str = std::str::from_utf8(&unpacked).unwrap_or("?");
    write!(writer, "{}\t{}", id, seq_str)?;
    for &d in depths {
        write!(writer, "\t{d}")?;
    }
    writeln!(writer)
}

fn fix_header_count(path: &str, count: u64) -> std::io::Result<()> {
    use std::io::{Seek, SeekFrom};
    let mut file = std::fs::OpenOptions::new().write(true).open(path)?;
    file.seek(SeekFrom::Start(0))?;
    write!(file, "#Number of markers : {:<20}", count)?;
    Ok(())
}

// === Chunk I/O: sort buffer -> lz4 compressed temp file ===

fn flush_chunk(
    buffer: &mut Vec<MergeEntry>,
    temp_dir: &tempfile::TempDir,
    chunk_idx: usize,
) -> Result<std::path::PathBuf, Box<dyn std::error::Error>> {
    buffer.sort_unstable_by(|a, b| a.packed_seq.cmp(&b.packed_seq));

    let chunk_path = temp_dir.path().join(format!("chunk_{:04}.lz4", chunk_idx));
    let file = std::fs::File::create(&chunk_path)?;
    let mut encoder = lz4_flex::frame::FrameEncoder::new(BufWriter::new(file));

    for entry in buffer.iter() {
        write_entry(&mut encoder, entry)?;
    }

    encoder.finish()?;
    buffer.clear();
    Ok(chunk_path)
}

/// Binary format: [seq_len:u16][packed_seq][offset:u16][n_depths:u16][depths:u16*n]
fn write_entry(w: &mut impl Write, entry: &MergeEntry) -> std::io::Result<()> {
    let seq_len = entry.packed_seq.len() as u16;
    w.write_all(&seq_len.to_le_bytes())?;
    w.write_all(&entry.packed_seq)?;
    w.write_all(&entry.sample_offset.to_le_bytes())?;
    let n_depths = entry.depths.len() as u16;
    w.write_all(&n_depths.to_le_bytes())?;
    for &d in &entry.depths {
        w.write_all(&d.to_le_bytes())?;
    }
    Ok(())
}

fn read_entry(r: &mut impl Read) -> std::io::Result<Option<MergeEntry>> {
    let mut buf2 = [0u8; 2];
    match r.read_exact(&mut buf2) {
        Ok(()) => {}
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(e) => return Err(e),
    }
    let seq_len = u16::from_le_bytes(buf2) as usize;

    let mut packed_seq = vec![0u8; seq_len];
    r.read_exact(&mut packed_seq)?;

    r.read_exact(&mut buf2)?;
    let sample_offset = u16::from_le_bytes(buf2);

    r.read_exact(&mut buf2)?;
    let n_depths = u16::from_le_bytes(buf2) as usize;

    let mut depths = vec![0u16; n_depths];
    for d in &mut depths {
        r.read_exact(&mut buf2)?;
        *d = u16::from_le_bytes(buf2);
    }

    Ok(Some(MergeEntry {
        packed_seq,
        sample_offset,
        depths,
    }))
}

struct ChunkReader {
    reader: lz4_flex::frame::FrameDecoder<BufReader<std::fs::File>>,
}

impl ChunkReader {
    fn open(path: &std::path::Path) -> std::io::Result<Self> {
        let file = std::fs::File::open(path)?;
        let reader = lz4_flex::frame::FrameDecoder::new(BufReader::new(file));
        Ok(ChunkReader { reader })
    }

    fn next_entry(&mut self) -> std::io::Result<Option<MergeEntry>> {
        read_entry(&mut self.reader)
    }
}
