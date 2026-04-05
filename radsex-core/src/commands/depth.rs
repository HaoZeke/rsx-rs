// GPL-3.0-or-later
// Copyright 2024--present radsex-rs developers

//! `depth` command: compute retained read statistics per individual.

use crate::markers_table::{MarkersTableStream, ParserConfig};
use crate::popmap::Popmap;
use crate::stats;
use std::io::Write;
use std::path::Path;

/// Parameters for the `depth` command.
pub struct DepthParams {
    pub markers_table_path: String,
    pub popmap_file_path: String,
    pub output_file_path: String,
    pub min_frequency: f32,
}

/// Run the `depth` analysis.
pub fn run(params: &DepthParams) -> Result<(), Box<dyn std::error::Error>> {
    let table_path = Path::new(&params.markers_table_path);
    let popmap = Popmap::from_file(Path::new(&params.popmap_file_path))?;

    let config = ParserConfig {
        store_sequence: false,
        compute_groups: false,
        min_depth: 1,
    };

    let stream = MarkersTableStream::open(table_path, Some(&popmap), config)?;
    let n_individuals = stream.header.n_individuals as usize;
    let min_individuals =
        (params.min_frequency * stream.header.n_individuals as f32) as u32;

    // Per-individual accumulators
    let mut depths: Vec<Vec<u16>> = vec![Vec::new(); n_individuals];
    let mut individual_markers_count: Vec<u64> = vec![0; n_individuals];
    let mut individual_reads_count: Vec<u64> = vec![0; n_individuals];

    for marker in stream.iter() {
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
    }

    // Check that we have retained markers
    if depths.iter().any(|d| d.is_empty()) {
        return Err(format!(
            "No markers were present in at least {}% of all individuals ({}/{} individuals)",
            (params.min_frequency * 100.0) as u32,
            min_individuals,
            n_individuals
        )
        .into());
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
        let group = popmap
            .get_group(individual_name)
            .unwrap_or("");

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
