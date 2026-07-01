// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! `freq` command: compute marker frequency distribution.

use crate::markers_table::{MarkersTableStream, ParserConfig};
use crate::source::MarkerStream;
use std::io::Write;
use std::path::Path;

/// Parameters for the `freq` command.
pub struct FreqParams {
    pub markers_table_path: String,
    pub output_file_path: String,
    pub min_depth: u16,
}

/// Run the `freq` analysis against a TSV path.
pub fn run(params: &FreqParams) -> Result<(), Box<dyn std::error::Error>> {
    let table_path = Path::new(&params.markers_table_path);
    let config = ParserConfig {
        store_sequence: false,
        store_depths: false,
        compute_groups: false,
        min_depth: params.min_depth,
    };
    let stream = MarkersTableStream::open(table_path, None, config)?;
    run_with_source(&stream, params)
}

/// Run the `freq` analysis against any `MarkerStream` (TSV, Arrow, Parquet).
///
/// Prefer a columnar histogram when the source provides one (Arrow batches as
/// tiles); otherwise fall back to the generic row-stream + bitset count path.
pub fn run_with_source<S: MarkerStream>(
    source: &S,
    params: &FreqParams,
) -> Result<(), Box<dyn std::error::Error>> {
    let frequency = compute_frequency(source, params.min_depth)?;
    write_frequency(&frequency, &params.output_file_path, params.min_depth)
}

/// Histogram of marker frequency (index = number of individuals with depth ≥
/// `min_depth`). Exposed for microbenchmarks / Arrow vs stream comparisons.
pub fn compute_frequency<S: MarkerStream>(
    source: &S,
    min_depth: u16,
) -> Result<Vec<u32>, Box<dyn std::error::Error>> {
    if let Some(col) = source.freq_histogram_columnar(min_depth) {
        return Ok(col?);
    }

    let n_individuals = source.header().n_individuals as usize;

    #[cfg(feature = "parallel")]
    {
        let frequency = source.par_fold_reduce(
            vec![0u32; n_individuals + 1],
            |freq, marker| {
                freq[marker.n_individuals as usize] += 1;
            },
            |mut a, b| {
                for (dst, src) in a.iter_mut().zip(b) {
                    *dst += src;
                }
                a
            },
        )?;
        return Ok(frequency);
    }

    #[cfg(not(feature = "parallel"))]
    {
        let mut frequency: Vec<u32> = vec![0; n_individuals + 1];
        source.for_each(|marker| {
            frequency[marker.n_individuals as usize] += 1;
        })?;
        Ok(frequency)
    }
}

fn write_frequency(
    frequency: &[u32],
    output_file_path: &str,
    min_depth: u16,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut output = std::fs::File::create(output_file_path)?;
    writeln!(output, "#source:rsx-freq;min_depth:{min_depth}")?;
    writeln!(output, "Frequency\tCount")?;

    for (i, count) in frequency.iter().enumerate().skip(1) {
        writeln!(output, "{i}\t{count}")?;
    }

    Ok(())
}
