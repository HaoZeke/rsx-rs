// GPL-3.0-or-later
// Copyright 2024--present radsex-rs developers

//! `signif` command: extract markers significantly associated with a group.

use crate::marker::Marker;
use crate::markers_table::{MarkersTableStream, ParserConfig};
use crate::popmap::{GroupConfig, Popmap};
use crate::stats;
use std::io::Write;
use std::path::Path;

/// Parameters for the `signif` command.
pub struct SignifParams {
    pub markers_table_path: String,
    pub popmap_file_path: String,
    pub output_file_path: String,
    pub min_depth: u16,
    pub signif_threshold: f32,
    pub disable_correction: bool,
    pub output_fasta: bool,
    pub group1: String,
    pub group2: String,
}

/// Run the `signif` analysis.
pub fn run(params: &SignifParams) -> Result<(), Box<dyn std::error::Error>> {
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
        store_sequence: true,
        compute_groups: true,
        min_depth: params.min_depth,
    };

    let stream = MarkersTableStream::open(table_path, Some(&popmap), config)?;
    let header_columns = stream.header.columns.clone();

    // First pass: collect candidate markers (p < threshold before correction)
    let mut candidate_markers: Vec<Marker> = Vec::new();
    let mut n_markers: u64 = 0;

    for mut marker in stream.iter() {
        if marker.n_individuals > 0 {
            n_markers += 1;
            let g1 = *marker.group_counts.get(&groups.group1).unwrap_or(&0);
            let g2 = *marker.group_counts.get(&groups.group2).unwrap_or(&0);
            marker.p = stats::p_association(g1, g2, total_g1, total_g2);

            if (marker.p as f32) < params.signif_threshold {
                candidate_markers.push(marker);
            }
        }
    }

    // Apply Bonferroni correction
    let effective_n_markers = if params.disable_correction {
        1u64
    } else {
        n_markers
    };

    let corrected_threshold = if params.disable_correction {
        params.signif_threshold as f64
    } else {
        params.signif_threshold as f64 / n_markers as f64
    };

    // Write output
    let mut output = std::io::BufWriter::new(std::fs::File::create(&params.output_file_path)?);

    if !params.output_fasta {
        writeln!(
            output,
            "#source:radsex-signif;min_depth:{};signif_threshold:{};bonferroni:{};n_markers:{}",
            params.min_depth,
            params.signif_threshold,
            !params.disable_correction,
            effective_n_markers
        )?;
        writeln!(output, "{}", header_columns.join("\t"))?;
    }

    // Second pass: filter with corrected threshold
    for mut marker in candidate_markers {
        if (marker.p as f32) < corrected_threshold as f32 {
            marker.p_corrected = stats::bonferroni_correct(marker.p, effective_n_markers);

            if params.output_fasta {
                marker.write_as_fasta(&mut output, params.min_depth as u32)?;
            } else {
                marker.write_as_table(&mut output)?;
            }
        }
    }

    Ok(())
}
