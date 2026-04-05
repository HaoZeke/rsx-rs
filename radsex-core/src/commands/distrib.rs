// GPL-3.0-or-later
// Copyright 2024--present radsex-rs developers

//! `distrib` command: compute marker distribution between two groups.

use crate::markers_table::{MarkersTableStream, ParserConfig};
use crate::popmap::{GroupConfig, Popmap};
use crate::stats;
use std::io::Write;
use std::path::Path;

/// Parameters for the `distrib` command.
pub struct DistribParams {
    pub markers_table_path: String,
    pub popmap_file_path: String,
    pub output_file_path: String,
    pub min_depth: u16,
    pub signif_threshold: f32,
    pub disable_correction: bool,
    pub group1: String,
    pub group2: String,
}

/// Run the `distrib` analysis.
pub fn run(params: &DistribParams) -> Result<(), Box<dyn std::error::Error>> {
    let table_path = Path::new(&params.markers_table_path);
    let popmap = Popmap::from_file(Path::new(&params.popmap_file_path))?;

    let mut groups = GroupConfig {
        group1: params.group1.clone(),
        group2: params.group2.clone(),
    };
    popmap.resolve_groups(&mut groups)?;

    let total_g1 = popmap.get_count(&groups.group1);
    let total_g2 = popmap.get_count(&groups.group2);

    let config = ParserConfig {
        store_sequence: false,
        compute_groups: true,
        min_depth: params.min_depth,
    };

    let stream = MarkersTableStream::open(table_path, Some(&popmap), config)?;

    // distribution[g1_count][g2_count] = marker_count
    let rows = (total_g1 + 1) as usize;
    let cols = (total_g2 + 1) as usize;
    let mut distribution: Vec<Vec<u64>> = vec![vec![0; cols]; rows];
    let mut n_markers: u64 = 0;

    for marker in stream.iter() {
        if marker.n_individuals > 0 {
            let g1 = *marker.group_counts.get(&groups.group1).unwrap_or(&0) as usize;
            let g2 = *marker.group_counts.get(&groups.group2).unwrap_or(&0) as usize;
            distribution[g1][g2] += 1;
            n_markers += 1;
        }
    }

    // Apply Bonferroni correction
    let effective_n_markers = if params.disable_correction {
        1u64
    } else {
        n_markers
    };

    let signif_threshold = if params.disable_correction {
        params.signif_threshold as f64
    } else {
        params.signif_threshold as f64 / n_markers as f64
    };

    // Write output
    let mut output = std::io::BufWriter::new(std::fs::File::create(&params.output_file_path)?);
    writeln!(
        output,
        "#source:radsex-distrib;min_depth:{};signif_threshold:{};bonferroni:{};n_markers:{}",
        params.min_depth,
        signif_threshold,
        !params.disable_correction,
        effective_n_markers
    )?;
    writeln!(
        output,
        "{}\t{}\tMarkers\tP\tCorrectedP\tSignif\tBias",
        groups.group1, groups.group2
    )?;

    for g in 0..=total_g1 {
        for h in 0..=total_g2 {
            if g + h == 0 {
                continue;
            }

            let count = distribution[g as usize][h as usize];
            let p = stats::p_association(g, h, total_g1, total_g2);
            let p_corrected = stats::bonferroni_correct(p, effective_n_markers);
            let signif = p < signif_threshold;
            let bias = stats::group_bias(g, total_g1, h, total_g2);

            writeln!(
                output,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                g,
                h,
                count,
                p,
                p_corrected,
                if signif { "True" } else { "False" },
                bias
            )?;
        }
    }

    Ok(())
}
