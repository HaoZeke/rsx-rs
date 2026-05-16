// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! `distrib` command: compute marker distribution between two groups.

use crate::bitset::GroupMask;
use crate::markers_table::{MarkersTableStream, ParserConfig};
use crate::popmap::{GroupConfig, Popmap};
use crate::stats;
use crate::stats::Cg;
use crate::test_method::{compute_p, CorrectionMethod, TestMethod};
use std::io::Write;
use std::path::Path;

/// Parameters for the `distrib` command.
pub struct DistribParams {
    pub markers_table_path: String,
    pub popmap_file_path: String,
    pub output_file_path: String,
    pub min_depth: u16,
    pub signif_threshold: f32,
    pub correction: CorrectionMethod,
    pub test_method: TestMethod,
    pub output_bayes: bool,
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
        store_depths: false,
        compute_groups: true,
        min_depth: params.min_depth,
    };

    let stream = MarkersTableStream::open(table_path, Some(&popmap), config)?;

    let rows = (total_g1 + 1) as usize;
    let cols = (total_g2 + 1) as usize;
    let mut distribution: Vec<Vec<u64>> = vec![vec![0; cols]; rows];
    let mut n_markers: u64 = 0;

    // Pre-compute group masks for popcount-based counting
    let mask_g1 =
        GroupMask::from_columns(&stream.groups, &groups.group1, stream.header.n_individuals);
    let mask_g2 =
        GroupMask::from_columns(&stream.groups, &groups.group2, stream.header.n_individuals);

    stream.for_each(|marker| {
        if marker.n_individuals > 0 {
            // popcount: O(1) per group instead of HashMap lookup
            let g1 = marker.presence.count_masked(&mask_g1) as usize;
            let g2 = marker.presence.count_masked(&mask_g2) as usize;
            distribution[g1][g2] += 1;
            n_markers += 1;
        }
    })?;

    let is_corrected = !matches!(params.correction, CorrectionMethod::None);
    let effective_n_markers = if is_corrected { n_markers } else { 1u64 };
    let signif_threshold = if is_corrected {
        params.signif_threshold as f64 / n_markers as f64
    } else {
        params.signif_threshold as f64
    };

    let mut output = std::io::BufWriter::new(std::fs::File::create(&params.output_file_path)?);
    writeln!(
        output,
        "#source:rsx-distrib;min_depth:{};signif_threshold:{};bonferroni:{};n_markers:{}",
        params.min_depth,
        Cg(signif_threshold),
        is_corrected,
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
            let p = compute_p(params.test_method, g, h, total_g1, total_g2);
            let p_corrected = stats::bonferroni_correct(p, effective_n_markers);
            let signif = p < signif_threshold;
            let bias = stats::group_bias(g, total_g1, h, total_g2);
            writeln!(
                output,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                g,
                h,
                count,
                Cg(p),
                Cg(p_corrected),
                if signif { "True" } else { "False" },
                Cg(bias)
            )?;
        }
    }
    Ok(())
}
