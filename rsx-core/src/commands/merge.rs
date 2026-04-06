// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Merge multiple marker depth tables by sequence identity.
//!
//! Uses 2-bit packed DNA keys (4x memory reduction) and mmap I/O.
//! Sequences are deduplicated by their packed representation.

use crate::io::seq_reader::{pack_2bit, unpack_2bit};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

pub struct MergeParams {
    pub input_files: Vec<String>,
    pub output_file_path: String,
}

pub fn run(params: &MergeParams) -> Result<(), Box<dyn std::error::Error>> {
    let n_inputs = params.input_files.len();
    eprintln!("Merging {} marker tables...", n_inputs);

    // Phase 1: read headers to get sample names
    let mut all_samples: Vec<String> = Vec::new();
    let mut file_samples: Vec<Vec<String>> = Vec::new();

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
            Path::new(path).file_name().unwrap_or_default().to_string_lossy(),
            samples.len(),
            n_markers
        );

        all_samples.extend(samples.clone());
        file_samples.push(samples);
    }

    let total_samples = all_samples.len();
    eprintln!("  Total samples: {}", total_samples);

    // Phase 2: read all files, build packed_seq -> depths map
    // 2-bit packed keys: 100bp -> 26 bytes vs 100 bytes (4x reduction)
    let mut seq_map: ahash::AHashMap<Vec<u8>, Vec<u16>> = ahash::AHashMap::new();
    let mut sample_offset = 0usize;

    for (i, path) in params.input_files.iter().enumerate() {
        let n_samples_this = file_samples[i].len();
        eprintln!("  Reading input {}...", i + 1);

        let file = std::fs::File::open(path)?;
        let reader = BufReader::with_capacity(1 << 20, file);
        let mut line_count = 0u64;

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

            let packed = pack_2bit(seq_str.as_bytes());
            let entry = seq_map
                .entry(packed)
                .or_insert_with(|| vec![0u16; total_samples]);

            for (j, depth) in depths_str.split('\t').enumerate() {
                if j < n_samples_this {
                    entry[sample_offset + j] = depth.parse().unwrap_or(0);
                }
            }

            line_count += 1;
            if line_count % 5_000_000 == 0 {
                eprintln!(
                    "    {} lines, {} unique sequences",
                    line_count,
                    seq_map.len()
                );
            }
        }

        eprintln!(
            "    Done: {} lines, {} unique sequences total",
            line_count,
            seq_map.len()
        );
        sample_offset += n_samples_this;
    }

    let n_merged = seq_map.len();
    eprintln!("  Total unique sequences: {}", n_merged);

    // Phase 3: write output (unpack 2-bit keys)
    eprintln!("  Writing {}...", params.output_file_path);
    let out_file = std::fs::File::create(&params.output_file_path)?;
    let mut writer = BufWriter::with_capacity(1 << 20, out_file);

    writeln!(writer, "#Number of markers : {}", n_merged)?;
    write!(writer, "id\tsequence")?;
    for s in &all_samples {
        write!(writer, "\t{}", s)?;
    }
    writeln!(writer)?;

    for (idx, (packed, depths)) in seq_map.iter().enumerate() {
        let unpacked = unpack_2bit(packed);
        let seq_str = std::str::from_utf8(&unpacked).unwrap_or("?");
        write!(writer, "{}\t{}", idx, seq_str)?;
        for d in depths {
            write!(writer, "\t{d}")?;
        }
        writeln!(writer)?;

        if (idx + 1) % 5_000_000 == 0 {
            eprintln!("    Written {} / {}", idx + 1, n_merged);
        }
    }

    eprintln!(
        "  Done: {} markers x {} samples -> {}",
        n_merged, total_samples, params.output_file_path
    );
    Ok(())
}
