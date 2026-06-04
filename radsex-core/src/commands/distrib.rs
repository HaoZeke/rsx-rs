// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! `distrib` command: compute marker distribution between two groups.

use crate::bitset::GroupMask;
use crate::markers_table::{MarkersTableStream, ParserConfig};
use crate::popmap::{GroupConfig, Popmap};
use crate::source::MarkerStream;
use crate::stats;
use crate::stats::Cg;
use crate::test_method::{CorrectionMethod, TestMethod, compute_p};
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

/// Run `distrib` against the on-disk markers TSV + popmap referenced by `params`.
pub fn run(params: &DistribParams) -> Result<(), Box<dyn std::error::Error>> {
    let table_path = Path::new(&params.markers_table_path);
    let popmap = Popmap::from_file(Path::new(&params.popmap_file_path))?;

    let config = ParserConfig {
        store_sequence: false,
        store_depths: false,
        compute_groups: true,
        min_depth: params.min_depth,
    };
    let stream = MarkersTableStream::open(table_path, Some(&popmap), config)?;
    run_with_source(&stream, &popmap, params)
}

/// Run `distrib` against any `MarkerStream`. Caller supplies the `Popmap`.
pub fn run_with_source<S: MarkerStream>(
    source: &S,
    popmap: &Popmap,
    params: &DistribParams,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut groups = GroupConfig {
        group1: params.group1.clone(),
        group2: params.group2.clone(),
    };
    popmap.resolve_groups(&mut groups)?;

    let total_g1 = popmap.get_count(&groups.group1);
    let total_g2 = popmap.get_count(&groups.group2);

    let rows = (total_g1 + 1) as usize;
    let cols = (total_g2 + 1) as usize;

    let n_individuals = source.header().n_individuals;
    let mask_g1 = GroupMask::from_columns(source.groups(), &groups.group1, n_individuals);
    let mask_g2 = GroupMask::from_columns(source.groups(), &groups.group2, n_individuals);

    #[cfg(feature = "parallel")]
    let (distribution, n_markers) = source.par_fold_reduce(
        (vec![vec![0u64; cols]; rows], 0u64),
        |(dist, n), marker| {
            if marker.n_individuals > 0 {
                let g1 = marker.presence.count_masked(&mask_g1) as usize;
                let g2 = marker.presence.count_masked(&mask_g2) as usize;
                dist[g1][g2] += 1;
                *n += 1;
            }
        },
        |(mut a, na), (b, nb)| {
            for (row_a, row_b) in a.iter_mut().zip(b) {
                for (cell_a, cell_b) in row_a.iter_mut().zip(row_b) {
                    *cell_a += cell_b;
                }
            }
            (a, na + nb)
        },
    )?;

    #[cfg(not(feature = "parallel"))]
    let mut distribution: Vec<Vec<u64>> = vec![vec![0; cols]; rows];
    #[cfg(not(feature = "parallel"))]
    let mut n_markers: u64 = 0;

    #[cfg(not(feature = "parallel"))]
    source.for_each(|marker| {
        if marker.n_individuals > 0 {
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
