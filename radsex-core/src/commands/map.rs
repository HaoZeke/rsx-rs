// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! `map` command: align markers to a reference genome and compute metrics.
//!
//! Pass 1: count markers for Bonferroni (fast, no alignment).
//! Pass 2: align each candidate marker, compute stats, and write in table order.

use crate::bitset::GroupMask;
use crate::markers_table::{MarkersTableStream, ParserConfig};
use crate::popmap::{GroupConfig, Popmap};
use crate::stats;
use crate::stats::Cg;
use minimap2::Aligner;
use std::collections::HashMap;
use std::io::{BufRead, Write};
use std::path::Path;

pub struct MapParams {
    pub markers_table_path: String,
    pub popmap_file_path: String,
    pub genome_file_path: String,
    pub output_file_path: String,
    pub min_depth: u16,
    pub min_quality: u32,
    pub min_frequency: f32,
    pub signif_threshold: f32,
    pub correction: crate::test_method::CorrectionMethod,
    pub test_method: crate::test_method::TestMethod,
    pub output_bayes: bool,
    pub group1: String,
    pub group2: String,
}

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
            current_contig = stripped.split_whitespace().next().unwrap_or("").to_string();
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

    let min_individuals = std::cmp::max(
        1,
        (params.min_frequency * popmap.n_individuals as f32) as u32,
    );

    let contig_lengths = load_contig_lengths(Path::new(&params.genome_file_path))?;

    // Pass 1: count markers for Bonferroni (no alignment, no sequence, fast)
    log::info!("map pass 1: counting markers");
    let config1 = ParserConfig {
        store_sequence: false,
        store_depths: false,
        compute_groups: true,
        min_depth: params.min_depth,
    };
    let stream1 = MarkersTableStream::open(table_path, Some(&popmap), config1)?;
    let n_markers = stream1.count_markers()?;
    log::info!("map pass 1: {} markers", n_markers);

    let effective_n_markers = if matches!(
        params.correction,
        crate::test_method::CorrectionMethod::None
    ) {
        1u64
    } else {
        n_markers
    };
    let signif_threshold = if matches!(
        params.correction,
        crate::test_method::CorrectionMethod::None
    ) {
        params.signif_threshold as f64
    } else {
        params.signif_threshold as f64 / n_markers as f64
    };

    // Build minimap2 index
    log::info!("Building minimap2 index for {}", params.genome_file_path);
    let aligner = Aligner::builder()
        .sr()
        .with_index(&params.genome_file_path, None)
        .map_err(|e| format!("Failed to build minimap2 index: {e}"))?;
    log::info!("Minimap2 index built");

    // Pass 2: align + write in table order
    log::info!("map pass 2: aligning and writing");
    let config2 = ParserConfig {
        store_sequence: true,
        store_depths: false,
        compute_groups: true,
        min_depth: params.min_depth,
    };
    let stream2 = MarkersTableStream::open(table_path, Some(&popmap), config2)?;

    let mask_g1 = GroupMask::from_columns(
        &stream2.groups,
        &groups.group1,
        stream2.header.n_individuals,
    );
    let mask_g2 = GroupMask::from_columns(
        &stream2.groups,
        &groups.group2,
        stream2.header.n_individuals,
    );

    let mut output = std::io::BufWriter::new(std::fs::File::create(&params.output_file_path)?);
    writeln!(
        output,
        "#source:rsx-map;min_depth:{};min_qual:{};min_freq:{};signif_threshold:{};bonferroni:{};n_markers:{}",
        params.min_depth,
        params.min_quality,
        params.min_frequency,
        Cg(signif_threshold),
        !matches!(
            params.correction,
            crate::test_method::CorrectionMethod::None
        ),
        effective_n_markers
    )?;
    writeln!(
        output,
        "Contig\tPosition\tLength\tMarker_id\tBias\tP\tCorrectedP\tSignif"
    )?;

    let mut n_aligned = 0u64;

    let mut write_alignment = |marker: &crate::marker::Marker| {
        let mappings = aligner
            .map(marker.sequence.as_bytes(), false, false, None, None, None)
            .unwrap_or_default();

        if mappings.is_empty() {
            return;
        }

        let mut best_idx = 0usize;
        let mut best_score = 0i32;
        let mut best_count = 0u32;

        for (j, mapping) in mappings.iter().enumerate() {
            let score = mapping
                .alignment
                .as_ref()
                .map_or(0, |a| a.alignment_score.unwrap_or(0));
            if score > best_score {
                best_idx = j;
                best_score = score;
                best_count = 1;
            } else if score == best_score {
                best_count += 1;
            }
        }

        let best = &mappings[best_idx];
        if best_count != 1 || best.mapq < params.min_quality {
            return;
        }

        let contig = best
            .target_name
            .as_deref()
            .map_or(String::new(), |v| v.to_string());
        let position = best.target_start;

        let g1 = marker.presence.count_masked(&mask_g1);
        let g2 = marker.presence.count_masked(&mask_g2);
        let bias = stats::group_bias(g1, total_g1, g2, total_g2);
        let p = stats::p_association(g1, g2, total_g1, total_g2);
        let p_corrected = stats::bonferroni_correct(p, effective_n_markers);
        let signif = p < signif_threshold;
        let contig_len = contig_lengths.get(&contig).copied().unwrap_or(0);

        let _ = writeln!(
            output,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            contig,
            position,
            contig_len,
            marker.id,
            Cg(bias),
            Cg(p),
            Cg(p_corrected),
            if signif { "True" } else { "False" }
        );

        n_aligned += 1;
    };

    #[cfg(feature = "parallel")]
    {
        let candidates = stream2.par_filter_map_collect(|marker| {
            if marker.n_individuals >= min_individuals {
                Some(marker.clone())
            } else {
                None
            }
        })?;

        for marker in candidates {
            write_alignment(&marker);
        }
    }

    #[cfg(not(feature = "parallel"))]
    stream2.for_each(|marker| {
        if marker.n_individuals >= min_individuals {
            write_alignment(marker);
        }
    })?;

    log::info!("Aligned {} markers to the reference genome", n_aligned);
    Ok(())
}
