// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! `subset` command: extract a filtered subset of markers.

use crate::bitset::GroupMask;
use crate::marker::Marker;
use crate::markers_table::{MarkersTableStream, ParserConfig};
use crate::popmap::{GroupConfig, Popmap};
use crate::stats;
use std::io::Write;
use std::path::Path;

pub struct SubsetParams {
    pub markers_table_path: String,
    pub popmap_file_path: String,
    pub output_file_path: String,
    pub min_depth: u16,
    pub signif_threshold: f32,
    pub disable_correction: bool,
    pub output_fasta: bool,
    pub group1: String,
    pub group2: String,
    pub min_group1: u32,
    pub min_group2: u32,
    pub max_group1: u32,
    pub max_group2: u32,
    pub min_individuals: u32,
    pub max_individuals: u32,
}

pub fn run(params: &SubsetParams) -> Result<(), Box<dyn std::error::Error>> {
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

    let mask_g1 = GroupMask::from_columns(&stream.groups, &groups.group1, stream.header.n_individuals);
    let mask_g2 = GroupMask::from_columns(&stream.groups, &groups.group2, stream.header.n_individuals);

    let mut filtered_markers: Vec<Marker> = Vec::new();
    let mut n_markers: u64 = 0;

    stream.for_each(|marker| {
        if marker.n_individuals > 0 { n_markers += 1; }
        let g1 = marker.presence.count_masked(&mask_g1);
        let g2 = marker.presence.count_masked(&mask_g2);
        if g1 >= params.min_group1 && g1 <= params.max_group1
            && g2 >= params.min_group2 && g2 <= params.max_group2
            && marker.n_individuals >= params.min_individuals
            && marker.n_individuals <= params.max_individuals
        {
            let mut m = marker.clone();
            m.p = stats::p_association(g1, g2, total_g1, total_g2);
            filtered_markers.push(m);
        }
    })?;

    let effective_n_markers = if params.disable_correction { 1u64 } else { n_markers };

    let mut output = std::io::BufWriter::new(std::fs::File::create(&params.output_file_path)?);

    if !params.output_fasta {
        writeln!(output,
            "#source:rsx-subset;min_depth:{};filters:{}=[{},{}],{}=[{},{}],individuals=[{},{}];signif_threshold:{};bonferroni:{};n_markers:{}",
            params.min_depth,
            groups.group1, params.min_group1, params.max_group1,
            groups.group2, params.min_group2, params.max_group2,
            params.min_individuals, params.max_individuals,
            params.signif_threshold,
            !params.disable_correction, effective_n_markers)?;
        writeln!(output, "{}", header_columns.join("\t"))?;
    }

    let fasta_groups = vec![
        (groups.group1.clone(), &mask_g1),
        (groups.group2.clone(), &mask_g2),
    ];

    for mut marker in filtered_markers {
        marker.p_corrected = stats::bonferroni_correct(marker.p, effective_n_markers);
        if params.output_fasta {
            marker.write_as_fasta_bitset(&mut output, params.min_depth as u32, &fasta_groups)?;
        } else {
            marker.write_as_table(&mut output)?;
        }
    }
    Ok(())
}
