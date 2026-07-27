// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! `distrib` command: compute marker distribution between two groups.

use crate::bitset::GroupMask;
use crate::markers_table::{MarkersTableStream, ParserConfig};
use crate::popmap::{GroupConfig, Popmap};
use crate::source::MarkerStream;
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

    // Build the ordered list of non-empty cells once so Bonferroni and BH FDR
    // share the same grid iteration.
    let mut cells: Vec<(u32, u32, u64, f64)> = Vec::new();
    for g in 0..=total_g1 {
        for h in 0..=total_g2 {
            if g + h == 0 {
                continue;
            }
            let count = distribution[g as usize][h as usize];
            let p = compute_p(params.test_method, g, h, total_g1, total_g2);
            cells.push((g, h, count, p));
        }
    }

    let n_tests = cells.len().max(1) as u64;
    let (p_corrected_vals, signif_flags, corr_label, threshold_note): (
        Vec<f64>,
        Vec<bool>,
        &str,
        f64,
    ) = match params.correction {
        CorrectionMethod::None => {
            let thr = params.signif_threshold as f64;
            let corr: Vec<f64> = cells.iter().map(|c| c.3).collect();
            let flags: Vec<bool> = cells.iter().map(|c| c.3 < thr).collect();
            (corr, flags, "none", thr)
        }
        CorrectionMethod::Bonferroni => {
            let thr = params.signif_threshold as f64 / n_markers.max(1) as f64;
            let corr: Vec<f64> = cells
                .iter()
                .map(|c| stats::bonferroni_correct(c.3, n_markers.max(1)))
                .collect();
            let flags: Vec<bool> = cells.iter().map(|c| c.3 < thr).collect();
            (corr, flags, "bonferroni", thr)
        }
        CorrectionMethod::Fdr => {
            let thr = params.signif_threshold as f64;
            let weighted: Vec<(f64, u64)> = cells.iter().map(|c| (c.3, c.2)).collect();
            let q = stats::benjamini_hochberg_weighted(&weighted);
            let flags: Vec<bool> = q.iter().map(|&qi| qi < thr).collect();
            (q, flags, "fdr", thr)
        }
    };

    let mut output = std::io::BufWriter::new(std::fs::File::create(&params.output_file_path)?);
    writeln!(
        output,
        "#source:rsx-distrib;min_depth:{};signif_threshold:{};correction:{};n_markers:{};n_tests:{}",
        params.min_depth,
        Cg(threshold_note),
        corr_label,
        n_markers,
        n_tests
    )?;
    if params.output_bayes {
        writeln!(
            output,
            "{}\t{}\tMarkers\tP\tCorrectedP\tSignif\tBias\tBayes_Factor\tPosterior_SexLinked",
            groups.group1, groups.group2
        )?;
    } else {
        writeln!(
            output,
            "{}\t{}\tMarkers\tP\tCorrectedP\tSignif\tBias",
            groups.group1, groups.group2
        )?;
    }

    for (i, &(g, h, count, p)) in cells.iter().enumerate() {
        let bias = stats::group_bias(g, total_g1, h, total_g2);
        write!(
            output,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            g,
            h,
            count,
            Cg(p),
            Cg(p_corrected_vals[i]),
            if signif_flags[i] { "True" } else { "False" },
            Cg(bias)
        )?;
        if params.output_bayes {
            let bf = stats::bayes_factor_2x2(g, h, total_g1, total_g2);
            let post = stats::posterior_sex_linked(g, h, total_g1, total_g2, 0.01, 0.9);
            writeln!(output, "\t{:.4}\t{:.4}", bf, post)?;
        } else {
            writeln!(output)?;
        }
    }
    Ok(())
}
