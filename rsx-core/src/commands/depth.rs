// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! `depth` command: compute retained read statistics per individual.
//!
//! Streaming: O(n_individuals) memory, not O(n_markers x n_individuals).
//! Uses online min/max/sum/count accumulators. Median approximated via
//! P-square algorithm (streaming quantile estimation).

use crate::markers_table::{MarkersTableStream, ParserConfig};
use crate::popmap::Popmap;
use std::io::Write;
use std::path::Path;

pub struct DepthParams {
    pub markers_table_path: String,
    pub popmap_file_path: String,
    pub output_file_path: String,
    pub min_frequency: f32,
}

/// Streaming statistics accumulator per individual.
struct StreamingStats {
    count: u64,
    sum: u64,
    min: u16,
    max: u16,
    // Reservoir sample for approximate median (size 1001)
    reservoir: Vec<u16>,
    reservoir_count: u64,
}

impl StreamingStats {
    fn new() -> Self {
        StreamingStats {
            count: 0,
            sum: 0,
            min: u16::MAX,
            max: 0,
            reservoir: Vec::with_capacity(1001),
            reservoir_count: 0,
        }
    }

    fn push(&mut self, d: u16) {
        self.count += 1;
        self.sum += d as u64;
        if d < self.min { self.min = d; }
        if d > self.max { self.max = d; }

        // Reservoir sampling for approximate median
        self.reservoir_count += 1;
        if self.reservoir.len() < 1001 {
            self.reservoir.push(d);
        } else {
            // Replace with probability 1001/reservoir_count
            let j = fastrand_u64(self.reservoir_count);
            if j < 1001 {
                self.reservoir[j as usize] = d;
            }
        }
    }

    fn median(&mut self) -> u16 {
        if self.reservoir.is_empty() {
            return 0;
        }
        self.reservoir.sort_unstable();
        self.reservoir[self.reservoir.len() / 2]
    }

    fn average(&self) -> u64 {
        self.sum.checked_div(self.count).unwrap_or(0)
    }
}

/// Simple LCG-based pseudo-random for reservoir sampling.
/// Not cryptographic, just needs to be uniform enough for sampling.
#[inline]
fn fastrand_u64(max: u64) -> u64 {
    // Use a simple hash of the count as the "random" value
    let mut x = max;
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94d049bb133111eb);
    x ^= x >> 31;
    x % max
}

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
    let min_individuals = (params.min_frequency * stream.header.n_individuals as f32) as u32;

    // Streaming accumulators: O(n_individuals) memory
    let mut retained_stats: Vec<StreamingStats> = (0..n_individuals)
        .map(|_| StreamingStats::new())
        .collect();
    let mut individual_markers_count: Vec<u64> = vec![0; n_individuals];
    let mut individual_reads_count: Vec<u64> = vec![0; n_individuals];

    stream.for_each(|marker| {
        for i in 0..n_individuals {
            let d = marker.individual_depths[i];
            if marker.n_individuals >= min_individuals {
                retained_stats[i].push(d);
            }
            if d > 0 {
                individual_markers_count[i] += 1;
                individual_reads_count[i] += d as u64;
            }
        }
    })?;

    if retained_stats.iter().any(|s| s.count == 0) {
        return Err(format!(
            "No markers were present in at least {}% of all individuals ({}/{} individuals)",
            (params.min_frequency * 100.0) as u32, min_individuals, n_individuals
        )
        .into());
    }

    let header_cols = &stream.header.columns;
    let mut output = std::fs::File::create(&params.output_file_path)?;
    writeln!(output,
        "Sample\tGroup\tReads\tMarkers\tRetained\tMin_depth\tMax_depth\tMedian_depth\tAverage_depth")?;

    for i in 0..n_individuals {
        let individual_name = &header_cols[i + 2];
        let group = popmap.get_group(individual_name).unwrap_or("");
        let median_d = retained_stats[i].median();
        writeln!(
            output,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            individual_name,
            group,
            individual_reads_count[i],
            individual_markers_count[i],
            retained_stats[i].count,
            retained_stats[i].min,
            retained_stats[i].max,
            median_d,
            retained_stats[i].average()
        )?;
    }
    Ok(())
}
