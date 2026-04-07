// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! `process` command: create marker depth table from demultiplexed reads.
//!
//! Single-phase concurrent merge: rayon threads insert directly into a
//! DashMap during file processing. No sequential merge bottleneck.
//! Sequences stored as 2-bit packed DNA (4x memory reduction).

use crate::io::seq_reader::{count_sequences, get_input_files, unpack_2bit};
use std::io::Write;
use std::path::Path;

/// Parameters for the `process` command.
pub struct ProcessParams {
    pub input_dir_path: String,
    pub output_file_path: String,
    pub n_threads: u32,
    pub min_depth: u16,
    /// If set, group markers by canonical k-mer of this size before output.
    /// Collapses sequencing error variants. Optional (default: disabled).
    pub kmer_dedup: Option<usize>,
}

/// Run the `process` command.
pub fn run(params: &ProcessParams) -> Result<(), Box<dyn std::error::Error>> {
    log::info!("RADSex process started");
    let start = std::time::Instant::now();

    let input_dir = Path::new(&params.input_dir_path);
    let input_files = get_input_files(input_dir)?;
    let n_individuals = input_files.len();

    // Build individual name -> column index mapping
    let individual_indices: ahash::AHashMap<String, usize> = input_files
        .iter()
        .enumerate()
        .map(|(i, f)| (f.individual_name.clone(), i))
        .collect();

    #[cfg(feature = "parallel")]
    let mut global = {
        use rayon::prelude::*;

        rayon::ThreadPoolBuilder::new()
            .num_threads(params.n_threads as usize)
            .build_global()
            .ok();

        // Use DashMap for concurrent merge when many files (>= 8),
        // fall back to collect+merge for fewer files where DashMap
        // sharding overhead dominates.
        if n_individuals >= 8 {
            use dashmap::DashMap;

            let dm: DashMap<Vec<u8>, Vec<u16>, ahash::RandomState> =
                DashMap::with_hasher(ahash::RandomState::new());

            input_files
                .par_iter()
                .for_each(|f| match count_sequences(&f.path) {
                    Ok(counts) => {
                        let idx = individual_indices[&f.individual_name];
                        for (packed_seq, count) in counts {
                            let mut entry = dm
                                .entry(packed_seq)
                                .or_insert_with(|| vec![0u16; n_individuals]);
                            entry[idx] = count;
                        }
                        log::info!("Finished processing individual {}", f.individual_name);
                    }
                    Err(e) => log::error!("Error processing {}: {e}", f.path.display()),
                });

            // Convert to AHashMap for uniform output path
            dm.into_iter().collect::<ahash::AHashMap<_, _>>()
        } else {
            // Few files: collect per-file results, merge sequentially
            let per_file: Vec<_> = input_files
                .par_iter()
                .filter_map(|f| {
                    count_sequences(&f.path).ok().map(|c| {
                        log::info!("Finished processing individual {}", f.individual_name);
                        (f.individual_name.clone(), c)
                    })
                })
                .collect();

            let mut global: ahash::AHashMap<Vec<u8>, Vec<u16>> = ahash::AHashMap::new();
            for (name, counts) in per_file {
                let idx = individual_indices[&name];
                for (packed_seq, count) in counts {
                    let entry = global
                        .entry(packed_seq)
                        .or_insert_with(|| vec![0u16; n_individuals]);
                    entry[idx] = count;
                }
            }
            global
        }
    };

    #[cfg(not(feature = "parallel"))]
    let mut global = {
        let mut global: ahash::AHashMap<Vec<u8>, Vec<u16>> = ahash::AHashMap::new();
        for f in &input_files {
            match count_sequences(&f.path) {
                Ok(counts) => {
                    let idx = individual_indices[&f.individual_name];
                    for (packed_seq, count) in counts {
                        let entry = global
                            .entry(packed_seq)
                            .or_insert_with(|| vec![0u16; n_individuals]);
                        entry[idx] = count;
                    }
                    log::info!("Finished processing individual {}", f.individual_name);
                }
                Err(e) => {
                    log::error!("Error processing {}: {e}", f.path.display());
                }
            }
        }
        global
    };

    // Optional k-mer deduplication: group markers by canonical k-mer
    if let Some(k) = params.kmer_dedup {
        log::info!(
            "K-mer deduplication (k={}): grouping {} markers",
            k,
            global.len()
        );
        let sequences: Vec<Vec<u8>> = global
            .keys()
            .map(|packed| crate::io::seq_reader::unpack_2bit(packed))
            .collect();
        let groups = crate::kmer::group_by_kmer(&sequences, k);
        let n_before = global.len();
        let n_groups = groups.len();
        log::info!(
            "K-mer dedup: {} markers -> {} groups ({:.1}% reduction)",
            n_before,
            n_groups,
            (1.0 - n_groups as f64 / n_before as f64) * 100.0
        );
        // For each group, keep the representative with highest total depth
        let keys: Vec<Vec<u8>> = global.keys().cloned().collect();
        let mut keep: ahash::AHashSet<usize> = ahash::AHashSet::new();
        for (_hash, indices) in &groups {
            let best = indices.iter().copied().max_by_key(|&i| {
                global
                    .get(&keys[i])
                    .map(|d| d.iter().map(|&v| v as u64).sum::<u64>())
                    .unwrap_or(0)
            });
            if let Some(b) = best {
                keep.insert(b);
            }
        }
        let kept_keys: ahash::AHashSet<Vec<u8>> = keep.iter().map(|&i| keys[i].clone()).collect();
        global.retain(|k, _| kept_keys.contains(k));
        log::info!("K-mer dedup: retained {} markers", global.len());
    }

    // Write output (unpack 2-bit keys back to ASCII for TSV)
    let mut output = std::io::BufWriter::new(std::fs::File::create(&params.output_file_path)?);
    writeln!(output, "#Number of markers : {}", global.len())?;

    write!(output, "id\tsequence")?;
    for f in &input_files {
        write!(output, "\t{}", f.individual_name)?;
    }
    writeln!(output)?;

    log::info!("Writing marker depths to output file");
    let mut id: u64 = 0;

    for (packed_seq, depths) in &global {
        if params.min_depth > 1 && !depths.iter().any(|&d| d >= params.min_depth) {
            continue;
        }
        let unpacked = unpack_2bit(packed_seq);
        let seq_str = std::str::from_utf8(&unpacked).unwrap_or("?");
        write!(output, "{}\t{}", id, seq_str)?;
        for &d in depths {
            write!(output, "\t{d}")?;
        }
        writeln!(output)?;
        id += 1;
    }

    let elapsed = start.elapsed();
    let secs = elapsed.as_secs();
    log::info!(
        "RADSex process ended (total runtime: {}h {}m {}s)",
        secs / 3600,
        (secs % 3600) / 60,
        secs % 60
    );

    Ok(())
}
