// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! MPI-distributed `process` command.
//!
//! Each MPI rank processes a subset of FASTQ files in parallel (rayon),
//! then results are reduced to rank 0 which writes the output.
//!
//! Usage: `mpirun -np 4 rsx process -i reads/ -o markers.tsv -T 4`
//! Build: `cargo build --release --features mpi`

#[cfg(feature = "mpi")]
use mpi::traits::*;

use crate::commands::process::ProcessParams;

#[cfg(feature = "mpi")]
use crate::io::seq_reader::{count_sequences, get_input_files};
#[cfg(feature = "mpi")]
use std::io::Write;

/// Run process with MPI distribution.
/// Falls back to single-node rayon if MPI is not initialized or size=1.
#[cfg(feature = "mpi")]
pub fn run_mpi(params: &ProcessParams) -> Result<(), Box<dyn std::error::Error>> {
    let universe = mpi::initialize().unwrap();
    let world = universe.world();
    let rank = world.rank() as usize;
    let size = world.size() as usize;

    if size <= 1 {
        // Single rank: fall back to normal process
        return crate::commands::process::run(params);
    }

    log::info!("MPI rank {rank}/{size}: rsx process started");
    let start = std::time::Instant::now();

    // Rank 0 discovers files
    let input_dir = std::path::Path::new(&params.input_dir_path);
    let all_files = if rank == 0 {
        get_input_files(input_dir)?
    } else {
        Vec::new()
    };

    // Broadcast file count
    let mut n_files = if rank == 0 { all_files.len() as u64 } else { 0 };
    world.process_at_rank(0).broadcast_into(&mut n_files);

    // Scatter file paths (serialize as newline-delimited strings)
    let my_files: Vec<InputFile> = if rank == 0 {
        // Distribute files round-robin
        all_files
            .into_iter()
            .enumerate()
            .filter(|(i, _)| i % size == rank)
            .map(|(_, f)| f)
            .collect()
    } else {
        // Other ranks need the full file list to pick their subset
        // Since we can't easily scatter variable-length data, each rank
        // re-discovers files and picks its subset
        let files = get_input_files(input_dir)?;
        files
            .into_iter()
            .enumerate()
            .filter(|(i, _)| i % size == rank)
            .map(|(_, f)| f)
            .collect()
    };

    log::info!(
        "MPI rank {rank}: processing {} files",
        my_files.len()
    );

    // Process local files with rayon
    #[cfg(feature = "parallel")]
    let local_results: Vec<(String, ahash::AHashMap<Vec<u8>, u16>)> = {
        use rayon::prelude::*;
        rayon::ThreadPoolBuilder::new()
            .num_threads(params.n_threads as usize)
            .build_global()
            .ok();

        my_files
            .par_iter()
            .filter_map(|f| {
                match count_sequences(&f.path) {
                    Ok(counts) => {
                        log::info!(
                            "MPI rank {}: finished {}",
                            rank,
                            f.individual_name
                        );
                        Some((f.individual_name.clone(), counts))
                    }
                    Err(e) => {
                        log::error!("MPI rank {}: error {}: {e}", rank, f.path.display());
                        None
                    }
                }
            })
            .collect()
    };

    #[cfg(not(feature = "parallel"))]
    let local_results: Vec<(String, ahash::AHashMap<Vec<u8>, u16>)> = my_files
        .iter()
        .filter_map(|f| {
            match count_sequences(&f.path) {
                Ok(counts) => Some((f.individual_name.clone(), counts)),
                Err(e) => {
                    log::error!("MPI rank {}: error {}: {e}", rank, f.path.display());
                    None
                }
            }
        })
        .collect();

    // Merge local results into a local global map
    let mut local_global: ahash::AHashMap<Vec<u8>, Vec<(String, u16)>> = ahash::AHashMap::new();
    for (name, counts) in &local_results {
        for (seq, &count) in counts {
            local_global
                .entry(seq.clone())
                .or_default()
                .push((name.clone(), count));
        }
    }

    // Serialize local results for MPI gather
    // Format: sequence_len(u32) + sequence + n_entries(u32) + [name_len(u32) + name + count(u16)]
    let local_bytes = serialize_counts(&local_global);
    let local_len = local_bytes.len() as i32;

    // Gather sizes to rank 0
    let all_sizes = if rank == 0 {
        let mut sizes = vec![0i32; size];
        world.process_at_rank(0).gather_into_root(&local_len, &mut sizes);
        sizes
    } else {
        world.process_at_rank(0).gather_into(&local_len);
        vec![]
    };

    // Gatherv data to rank 0
    if rank == 0 {
        let total: usize = all_sizes.iter().map(|&s| s as usize).sum();
        let mut all_bytes = vec![0u8; total];
        let displs: Vec<i32> = all_sizes
            .iter()
            .scan(0i32, |acc, &s| {
                let d = *acc;
                *acc += s;
                Some(d)
            })
            .collect();

        let mut partition =
            mpi::datatype::Partition::new(&mut all_bytes, all_sizes.clone(), displs);
        world
            .process_at_rank(0)
            .gather_varcount_into_root(&local_bytes, &mut partition);

        // Deserialize and merge all results
        let all_files_sorted = get_input_files(input_dir)?;
        let n_individuals = all_files_sorted.len();
        let individual_indices: ahash::AHashMap<String, usize> = all_files_sorted
            .iter()
            .enumerate()
            .map(|(i, f)| (f.individual_name.clone(), i))
            .collect();

        let mut global: ahash::AHashMap<Vec<u8>, Vec<u16>> = ahash::AHashMap::new();

        // Deserialize each rank's contribution
        let mut offset = 0;
        for &sz in &all_sizes {
            let chunk = &all_bytes[offset..offset + sz as usize];
            let counts = deserialize_counts(chunk);
            for (seq, entries) in counts {
                let depths = global
                    .entry(seq)
                    .or_insert_with(|| vec![0u16; n_individuals]);
                for (name, count) in entries {
                    if let Some(&idx) = individual_indices.get(&name) {
                        depths[idx] = count;
                    }
                }
            }
            offset += sz as usize;
        }

        // Write output
        let mut output =
            std::io::BufWriter::new(std::fs::File::create(&params.output_file_path)?);
        writeln!(output, "#Number of markers : {}", global.len())?;
        write!(output, "id\tsequence")?;
        for f in &all_files_sorted {
            write!(output, "\t{}", f.individual_name)?;
        }
        writeln!(output)?;

        let mut id: u64 = 0;
        for (seq, depths) in &global {
            if params.min_depth > 1 && !depths.iter().any(|&d| d >= params.min_depth) {
                continue;
            }
            let seq_str = std::str::from_utf8(seq).unwrap_or("?");
            write!(output, "{}\t{}", id, seq_str)?;
            for &d in depths {
                write!(output, "\t{d}")?;
            }
            writeln!(output)?;
            id += 1;
        }

        let elapsed = start.elapsed();
        log::info!(
            "MPI process done: {} markers, {}h {}m {}s",
            id,
            elapsed.as_secs() / 3600,
            (elapsed.as_secs() % 3600) / 60,
            elapsed.as_secs() % 60
        );
    } else {
        world
            .process_at_rank(0)
            .gather_varcount_into(&local_bytes);
    }

    Ok(())
}

#[cfg(feature = "mpi")]
fn serialize_counts(
    counts: &ahash::AHashMap<Vec<u8>, Vec<(String, u16)>>,
) -> Vec<u8> {
    let mut buf = Vec::new();
    for (seq, entries) in counts {
        buf.extend_from_slice(&(seq.len() as u32).to_le_bytes());
        buf.extend_from_slice(seq);
        buf.extend_from_slice(&(entries.len() as u32).to_le_bytes());
        for (name, count) in entries {
            buf.extend_from_slice(&(name.len() as u32).to_le_bytes());
            buf.extend_from_slice(name.as_bytes());
            buf.extend_from_slice(&count.to_le_bytes());
        }
    }
    buf
}

#[cfg(feature = "mpi")]
fn deserialize_counts(data: &[u8]) -> Vec<(Vec<u8>, Vec<(String, u16)>)> {
    let mut result = Vec::new();
    let mut pos = 0;
    while pos + 4 <= data.len() {
        let seq_len = u32::from_le_bytes(data[pos..pos + 4].try_into().unwrap()) as usize;
        pos += 4;
        if pos + seq_len > data.len() {
            break;
        }
        let seq = data[pos..pos + seq_len].to_vec();
        pos += seq_len;

        if pos + 4 > data.len() {
            break;
        }
        let n_entries = u32::from_le_bytes(data[pos..pos + 4].try_into().unwrap()) as usize;
        pos += 4;

        let mut entries = Vec::with_capacity(n_entries);
        for _ in 0..n_entries {
            if pos + 4 > data.len() {
                break;
            }
            let name_len = u32::from_le_bytes(data[pos..pos + 4].try_into().unwrap()) as usize;
            pos += 4;
            if pos + name_len > data.len() {
                break;
            }
            let name = String::from_utf8_lossy(&data[pos..pos + name_len]).into_owned();
            pos += name_len;
            if pos + 2 > data.len() {
                break;
            }
            let count = u16::from_le_bytes(data[pos..pos + 2].try_into().unwrap());
            pos += 2;
            entries.push((name, count));
        }
        result.push((seq, entries));
    }
    result
}

/// Stub for non-MPI builds.
#[cfg(not(feature = "mpi"))]
pub fn run_mpi(params: &ProcessParams) -> Result<(), Box<dyn std::error::Error>> {
    log::warn!("MPI not enabled. Build with --features mpi. Falling back to single-node.");
    crate::commands::process::run(params)
}
