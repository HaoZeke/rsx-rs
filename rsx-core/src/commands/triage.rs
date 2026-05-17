// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! `triage` command: marker-level biological candidate ranking.
//!
//! The command keeps RADSex-style strict testing and Bayesian marker evidence
//! in one bounded-memory pass over the marker table.

use crate::bitset::GroupMask;
use crate::markers_table::{MarkersTableStream, ParserConfig};
use crate::popmap::{GroupConfig, Popmap};
use crate::stats;
use crate::stats::Cg;
use crate::test_method::{TestMethod, compute_p};
use std::io::Write;
use std::path::Path;

pub struct TriageParams {
    pub markers_table_path: String,
    pub popmap_file_path: String,
    pub output_file_path: String,
    pub min_depth: u16,
    pub signif_threshold: f32,
    pub posterior_threshold: f64,
    pub bayes_factor_threshold: f64,
    pub prior_probability: f64,
    pub linked_probability: f64,
    pub group1: String,
    pub group2: String,
}

fn penetrance(present: u32, total: u32) -> f64 {
    if total == 0 {
        0.0
    } else {
        present as f64 / total as f64
    }
}

fn bias_direction(group1: &str, group2: &str, g1_penetrance: f64, g2_penetrance: f64) -> String {
    if g1_penetrance > g2_penetrance {
        format!("{group1}-biased")
    } else if g2_penetrance > g1_penetrance {
        format!("{group2}-biased")
    } else {
        "balanced".to_string()
    }
}

fn candidate_class(
    strict_call: bool,
    posterior_call: bool,
    bayes_factor_call: bool,
) -> &'static str {
    if strict_call && posterior_call {
        "strict+posterior"
    } else if strict_call {
        "strict_only"
    } else if posterior_call {
        "posterior_only"
    } else if bayes_factor_call {
        "bayes_factor_only"
    } else {
        "exploratory"
    }
}

pub fn run(params: &TriageParams) -> Result<(), Box<dyn std::error::Error>> {
    let table_path = Path::new(&params.markers_table_path);
    let popmap = Popmap::from_file(Path::new(&params.popmap_file_path))?;

    let mut groups = GroupConfig {
        group1: params.group1.clone(),
        group2: params.group2.clone(),
    };
    popmap.resolve_groups(&mut groups)?;

    let total_g1 = popmap.get_count(&groups.group1);
    let total_g2 = popmap.get_count(&groups.group2);

    let config1 = ParserConfig {
        store_sequence: false,
        store_depths: false,
        compute_groups: true,
        min_depth: params.min_depth,
    };
    let stream1 = MarkersTableStream::open(table_path, Some(&popmap), config1)?;
    let n_markers = stream1.count_markers()?;
    let strict_threshold = if n_markers > 0 {
        params.signif_threshold as f64 / n_markers as f64
    } else {
        params.signif_threshold as f64
    };

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
        "#source:rsx-triage;min_depth:{};signif_threshold:{};posterior_threshold:{};bayes_factor_threshold:{};prior_probability:{};linked_probability:{};n_markers:{}",
        params.min_depth,
        params.signif_threshold,
        Cg(params.posterior_threshold),
        Cg(params.bayes_factor_threshold),
        Cg(params.prior_probability),
        Cg(params.linked_probability),
        n_markers
    )?;
    writeln!(
        output,
        "id\tsequence\tGroup1\tGroup1_Present\tGroup1_Total\tGroup1_Penetrance\tGroup2\tGroup2_Present\tGroup2_Total\tGroup2_Penetrance\tBias_Direction\tBias\tP\tCorrectedP\tBayes_Factor\tPosterior_SexLinked\tStrict_Call\tPosterior_Call\tBayes_Factor_Call\tCandidate_Class"
    )?;

    stream2.for_each(|marker| {
        if marker.n_individuals == 0 {
            return;
        }

        let g1 = marker.presence.count_masked(&mask_g1);
        let g2 = marker.presence.count_masked(&mask_g2);
        let p = compute_p(TestMethod::ChiSquared, g1, g2, total_g1, total_g2);
        let p_corrected = stats::bonferroni_correct(p, n_markers);
        let strict_call = p < strict_threshold;
        let bf = stats::bayes_factor_2x2(g1, g2, total_g1, total_g2);
        let posterior = stats::posterior_sex_linked(
            g1,
            g2,
            total_g1,
            total_g2,
            params.prior_probability,
            params.linked_probability,
        );
        let posterior_call = posterior > params.posterior_threshold;
        let bayes_factor_call = bf > params.bayes_factor_threshold;
        if !(strict_call || posterior_call || bayes_factor_call) {
            return;
        }

        let g1_penetrance = penetrance(g1, total_g1);
        let g2_penetrance = penetrance(g2, total_g2);
        let bias = stats::group_bias(g1, total_g1, g2, total_g2);
        let class = candidate_class(strict_call, posterior_call, bayes_factor_call);
        let direction =
            bias_direction(&groups.group1, &groups.group2, g1_penetrance, g2_penetrance);

        let _ = writeln!(
            output,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            marker.id,
            marker.sequence,
            groups.group1,
            g1,
            total_g1,
            Cg(g1_penetrance),
            groups.group2,
            g2,
            total_g2,
            Cg(g2_penetrance),
            direction,
            Cg(bias),
            Cg(p),
            Cg(p_corrected),
            Cg(bf),
            Cg(posterior),
            strict_call,
            posterior_call,
            bayes_factor_call,
            class
        );
    })?;

    Ok(())
}
