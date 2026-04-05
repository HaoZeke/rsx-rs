// GPL-3.0-or-later
// Copyright 2024--present radsex-rs developers

//! `freq` command: compute marker frequency distribution.

use crate::markers_table::{MarkersTableStream, ParserConfig};
use std::io::Write;
use std::path::Path;

/// Parameters for the `freq` command.
pub struct FreqParams {
    pub markers_table_path: String,
    pub output_file_path: String,
    pub min_depth: u16,
}

/// Run the `freq` analysis.
pub fn run(params: &FreqParams) -> Result<(), Box<dyn std::error::Error>> {
    let table_path = Path::new(&params.markers_table_path);
    let config = ParserConfig {
        store_sequence: false,
        compute_groups: false,
        min_depth: params.min_depth,
    };

    let stream = MarkersTableStream::open(table_path, None, config)?;
    let n_individuals = stream.header.n_individuals as usize;

    // frequency[i] = number of markers present in exactly i individuals
    let mut frequency: Vec<u32> = vec![0; n_individuals + 1];

    for marker in stream.iter() {
        frequency[marker.n_individuals as usize] += 1;
    }

    // Write output
    let mut output = std::fs::File::create(&params.output_file_path)?;
    writeln!(
        output,
        "#source:radsex-freq;min_depth:{}",
        params.min_depth
    )?;
    writeln!(output, "Frequency\tCount")?;

    for i in 1..=n_individuals {
        writeln!(output, "{}\t{}", i, frequency[i])?;
    }

    Ok(())
}
