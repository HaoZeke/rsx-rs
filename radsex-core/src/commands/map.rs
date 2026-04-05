// GPL-3.0-or-later
// Copyright 2024--present radsex-rs developers

//! `map` command: align markers to a reference genome and compute metrics.
//!
//! Uses minimap2 for short-read alignment (replacing the C++ BWA-MEM integration).

use crate::marker::AlignedMarker;
use crate::markers_table::{MarkersTableStream, ParserConfig};
use crate::popmap::{GroupConfig, Popmap};
use crate::stats;
use minimap2::Aligner;
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
        if let Some(stripped) = line.strip_prefix('>') {
            if !current_contig.is_empty() {
                lengths.insert(current_contig.clone(), current_length);
            }
            current_contig = stripped
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

/// Run the `map` analysis with minimap2 alignment.
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

    // Build minimap2 index for short-read alignment
    log::info!("Building minimap2 index for {}", params.genome_file_path);
    let aligner = Aligner::builder()
        .sr() // short-read preset (equivalent to BWA-MEM for short reads)
        .with_index(&params.genome_file_path, None)
        .map_err(|e| format!("Failed to build minimap2 index: {e}"))?;

    log::info!("Minimap2 index built successfully");

    let config = ParserConfig {
        store_sequence: true,
        compute_groups: true,
        min_depth: params.min_depth,
    };

    let stream = MarkersTableStream::open(table_path, Some(&popmap), config)?;

    let mut aligned_markers: Vec<AlignedMarker> = Vec::new();
    let mut n_markers: u64 = 0;

    for marker in stream.iter() {
        if marker.n_individuals > 0 {
            n_markers += 1;
        }

        if marker.n_individuals < min_individuals {
            continue;
        }

        // Align marker sequence to the reference genome
        let mappings = aligner
            .map(marker.sequence.as_bytes(), false, false, None, None, None)
            .unwrap_or_default();

        // Find best-scoring unique alignment (matching C++ BWA-MEM logic)
        if mappings.is_empty() {
            continue;
        }

        let mut best_idx = 0usize;
        let mut best_score = 0i32;
        let mut best_count = 0u32;

        for (j, mapping) in mappings.iter().enumerate() {
            let score = mapping.alignment.as_ref().map_or(0, |a| a.alignment_score.unwrap_or(0));
            if score > best_score {
                best_idx = j;
                best_score = score;
                best_count = 1;
            } else if score == best_score {
                best_count += 1;
            }
        }

        let best = &mappings[best_idx];
        let mapq = best.mapq;

        // Retain only unique best alignment with mapq >= min_quality
        if best_count != 1 || mapq < params.min_quality {
            continue;
        }

        let contig = best
            .target_name
            .as_deref()
            .map_or(String::new(), |v| v.to_string());
        let position = best.target_start as i64;

        let g1 = *marker.group_counts.get(&groups.group1).unwrap_or(&0);
        let g2 = *marker.group_counts.get(&groups.group2).unwrap_or(&0);

        aligned_markers.push(AlignedMarker {
            id: marker.id.clone(),
            contig,
            position,
            bias: stats::group_bias(g1, total_g1, g2, total_g2),
            p: stats::p_association(g1, g2, total_g1, total_g2),
        });
    }

    log::info!(
        "Aligned {} markers to the reference genome",
        aligned_markers.len()
    );

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
    writeln!(
        output,
        "Contig\tPosition\tLength\tMarker_id\tBias\tP\tCorrectedP\tSignif"
    )?;

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
