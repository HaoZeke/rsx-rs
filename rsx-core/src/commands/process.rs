// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! `process` command: create marker depth table from demultiplexed reads.

use crate::io::seq_reader::{InputFile, count_sequences, get_input_files};
use std::io::Write;
use std::path::Path;

/// Parameters for the `process` command.
pub struct ProcessParams {
    pub input_dir_path: String,
    pub output_file_path: String,
    pub n_threads: u32,
    pub min_depth: u16,
}

/// Run the `process` command.
pub fn run(params: &ProcessParams) -> Result<(), Box<dyn std::error::Error>> {
    log::info!("RADSex process started");
    let start = std::time::Instant::now();

    let input_dir = Path::new(&params.input_dir_path);
    let input_files = get_input_files(input_dir)?;
    // Process files in parallel using rayon
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;

        rayon::ThreadPoolBuilder::new()
            .num_threads(params.n_threads as usize)
            .build_global()
            .ok(); // Ignore if already initialized

        let per_file_results: Vec<(String, ahash::AHashMap<Vec<u8>, u16>)> = input_files
            .par_iter()
            .filter_map(|f| {
                match count_sequences(&f.path) {
                    Ok(counts) => {
                        log::info!("Finished processing individual {}", f.individual_name);
                        Some((f.individual_name.clone(), counts))
                    }
                    Err(e) => {
                        log::error!("Error processing {}: {e}", f.path.display());
                        None
                    }
                }
            })
            .collect();

        merge_and_write(params, &input_files, per_file_results)?;
    }

    #[cfg(not(feature = "parallel"))]
    {
        let per_file_results: Vec<(String, ahash::AHashMap<Vec<u8>, u16>)> = input_files
            .iter()
            .filter_map(|f| {
                match count_sequences(&f.path) {
                    Ok(counts) => {
                        log::info!("Finished processing individual {}", f.individual_name);
                        Some((f.individual_name.clone(), counts))
                    }
                    Err(e) => {
                        log::error!("Error processing {}: {e}", f.path.display());
                        None
                    }
                }
            })
            .collect();

        merge_and_write(params, &input_files, per_file_results)?;
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

fn merge_and_write(
    params: &ProcessParams,
    input_files: &[InputFile],
    per_file_results: Vec<(String, ahash::AHashMap<Vec<u8>, u16>)>,
) -> Result<(), Box<dyn std::error::Error>> {
    // Build individual name -> index mapping
    let individual_indices: ahash::AHashMap<String, usize> = input_files
        .iter()
        .enumerate()
        .map(|(i, f)| (f.individual_name.clone(), i))
        .collect();

    let n_individuals = input_files.len();

    // Merge all per-file results into a global map: sequence -> depths[n_individuals]
    let mut global: ahash::AHashMap<Vec<u8>, Vec<u16>> = ahash::AHashMap::new();

    for (individual_name, counts) in &per_file_results {
        let idx = individual_indices[individual_name];
        for (seq, &count) in counts {
            let entry = global
                .entry(seq.clone())
                .or_insert_with(|| vec![0u16; n_individuals]);
            entry[idx] = count;
        }
    }

    // Filter by min_depth and write output
    let mut output = std::io::BufWriter::new(std::fs::File::create(&params.output_file_path)?);
    writeln!(output, "#Number of markers : {}", global.len())?;

    write!(output, "id\tsequence")?;
    for f in input_files {
        write!(output, "\t{}", f.individual_name)?;
    }
    writeln!(output)?;

    log::info!("Writing marker depths to output file");
    let mut id: u64 = 0;

    for (seq, depths) in &global {
        // Check min_depth filter
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

    Ok(())
}
