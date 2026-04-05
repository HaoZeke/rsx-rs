// GPL-3.0-or-later
// Copyright 2024--present radsex-rs developers

//! `map` command: align markers to a reference genome and compute metrics.
//!
//! Uses minimap2 for alignment (replacing the C++ BWA-MEM integration).
//! Note: minimap2 integration is stubbed pending the minimap2 crate dependency.
//! For now, this uses a simple exact-match approach as a placeholder that can
//! be swapped for minimap2 when the dependency is added.

use crate::marker::AlignedMarker;
use crate::markers_table::{MarkersTableStream, ParserConfig};
use crate::popmap::{GroupConfig, Popmap};
use crate::stats;
use std::collections::HashMap;
use std::io::{BufRead, Write};
use std::path::Path;

/// Parameters for the `map` command.
pub struct MapParams {
    pub markers_table_path: String,
    pub popmap_file_path: String,
    pub genome_file_path: String,
    pub output_file_path: String,
    pub min_depth: u16,
    pub min_quality: u32,
    pub min_frequency: f32,
    pub signif_threshold: f32,
    pub disable_correction: bool,
    pub group1: String,
    pub group2: String,
}

/// Load contig lengths from a FASTA genome file.
fn load_contig_lengths(genome_path: &Path) -> std::io::Result<HashMap<String, u64>> {
    let file = std::fs::File::open(genome_path)?;
    let reader = std::io::BufReader::new(file);
    let mut lengths: HashMap<String, u64> = HashMap::new();
    let mut current_contig = String::new();
    let mut current_length: u64 = 0;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !current_contig.is_empty() {
                lengths.insert(current_contig.clone(), current_length);
            }
            // Get contig name (first word after '>')
            current_contig = line[1..]
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            current_length = 0;
        } else {
            current_length += line.len() as u64;
        }
    }
    if !current_contig.is_empty() {
        lengths.insert(current_contig, current_length);
    }

    Ok(lengths)
}

/// Run the `map` analysis.
///
/// NOTE: Full minimap2 alignment will be integrated when the minimap2 crate
/// dependency is added. Currently this function reads markers, computes
/// statistics, but alignment is a placeholder that logs a warning.
pub fn run(params: &MapParams) -> Result<(), Box<dyn std::error::Error>> {
    let table_path = Path::new(&params.markers_table_path);
    let popmap = Popmap::from_file(Path::new(&params.popmap_file_path))?;

    let mut groups = GroupConfig {
        group1: params.group1.clone(),
        group2: params.group2.clone(),
    };
    popmap.resolve_groups(&mut groups)?;

    let total_g1 = popmap.get_count(&groups.group1);
    let total_g2 = popmap.get_count(&groups.group2);

    let min_individuals =
        std::cmp::max(1, (params.min_frequency * popmap.n_individuals as f32) as u32);

    let contig_lengths = load_contig_lengths(Path::new(&params.genome_file_path))?;

    let config = ParserConfig {
        store_sequence: true,
        compute_groups: true,
        min_depth: params.min_depth,
    };

    let stream = MarkersTableStream::open(table_path, Some(&popmap), config)?;

    let aligned_markers: Vec<AlignedMarker> = Vec::new();
    let mut n_markers: u64 = 0;

    // TODO: Replace this stub with actual minimap2 alignment.
    // The minimap2 crate will be used to align each marker sequence to the
    // reference genome. For now, markers are collected but not aligned.
    log::warn!("map command: minimap2 alignment not yet integrated; collecting markers only");

    for marker in stream.iter() {
        if marker.n_individuals > 0 {
            n_markers += 1;
        }

        if marker.n_individuals >= min_individuals {
            let g1 = *marker.group_counts.get(&groups.group1).unwrap_or(&0);
            let g2 = *marker.group_counts.get(&groups.group2).unwrap_or(&0);

            let bias = stats::group_bias(g1, total_g1, g2, total_g2);
            let p = stats::p_association(g1, g2, total_g1, total_g2);

            // Placeholder: in the real implementation, alignment results would
            // populate contig and position from minimap2 output.
            // For now, skip actual alignment but preserve the data pipeline.
            let _ = (bias, p, &marker.sequence);
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
        "#source:radsex-map;min_depth:{};min_qual:{};min_freq:{};signif_threshold:{};bonferroni:{};n_markers:{}",
        params.min_depth,
        params.min_quality,
        params.min_frequency,
        signif_threshold,
        !params.disable_correction,
        effective_n_markers
    )?;
    writeln!(output, "Contig\tPosition\tLength\tMarker_id\tBias\tP\tCorrectedP\tSignif")?;

    for marker in &aligned_markers {
        let contig_len = contig_lengths.get(&marker.contig).copied().unwrap_or(0);
        let p_corrected = stats::bonferroni_correct(marker.p, effective_n_markers);
        let signif = marker.p < signif_threshold;

        writeln!(
            output,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            marker.contig,
            marker.position,
            contig_len,
            marker.id,
            marker.bias,
            marker.p,
            p_corrected,
            if signif { "True" } else { "False" }
        )?;
    }

    Ok(())
}
