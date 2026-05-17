// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! `triage` command: marker-level biological candidate ranking.
//!
//! The command keeps RADSex-style strict testing and Bayesian marker evidence
//! in one bounded-memory pass over the marker table.

use crate::bitset::GroupMask;
use crate::markers_table::{MarkersTableStream, ParserConfig};
use crate::popmap::{GroupConfig, Popmap};
use crate::source::MarkerStream;
use crate::stats;
use crate::stats::Cg;
use crate::test_method::{TestMethod, compute_p};
use std::io::Write;
use std::path::Path;

#[cfg(feature = "arrow-output")]
use arrow::record_batch::RecordBatch;

#[derive(Clone)]
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
    let config = ParserConfig {
        store_sequence: true,
        store_depths: false,
        compute_groups: true,
        min_depth: params.min_depth,
    };
    let stream = MarkersTableStream::open(table_path, Some(&popmap), config)?;
    run_with_source(&stream, &popmap, params)
}

pub fn run_with_source<S: MarkerStream>(
    source: &S,
    popmap: &Popmap,
    params: &TriageParams,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut groups = GroupConfig {
        group1: params.group1.clone(),
        group2: params.group2.clone(),
    };
    popmap.resolve_groups(&mut groups)?;

    let total_g1 = popmap.get_count(&groups.group1);
    let total_g2 = popmap.get_count(&groups.group2);

    let n_markers = source.count_markers()?;
    let strict_threshold = if n_markers > 0 {
        params.signif_threshold as f64 / n_markers as f64
    } else {
        params.signif_threshold as f64
    };

    let n_individuals = source.header().n_individuals;
    let mask_g1 = GroupMask::from_columns(source.groups(), &groups.group1, n_individuals);
    let mask_g2 = GroupMask::from_columns(source.groups(), &groups.group2, n_individuals);

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

    source.for_each(|marker| {
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

#[cfg(feature = "arrow-output")]
pub fn run_to_arrow(params: &TriageParams) -> Result<Vec<RecordBatch>, Box<dyn std::error::Error>> {
    let table_path = Path::new(&params.markers_table_path);
    let popmap = Popmap::from_file(Path::new(&params.popmap_file_path))?;
    let config = ParserConfig {
        store_sequence: true,
        store_depths: false,
        compute_groups: true,
        min_depth: params.min_depth,
    };
    let stream = MarkersTableStream::open(table_path, Some(&popmap), config)?;
    run_to_arrow_with_source(&stream, &popmap, params)
}

#[cfg(feature = "arrow-output")]
pub fn run_to_arrow_with_source<S: MarkerStream>(
    source: &S,
    popmap: &Popmap,
    params: &TriageParams,
) -> Result<Vec<RecordBatch>, Box<dyn std::error::Error>> {
    use arrow::array::builder::{
        BooleanBuilder, Float64Builder, StringBuilder, UInt32Builder,
    };
    use arrow::datatypes::{DataType, Field, Schema};

    let schema = Schema::new(vec![
        Field::new("id", DataType::Utf8, false),
        Field::new("sequence", DataType::Utf8, false),
        Field::new("Group1", DataType::Utf8, false),
        Field::new("Group1_Present", DataType::UInt32, false),
        Field::new("Group1_Total", DataType::UInt32, false),
        Field::new("Group1_Penetrance", DataType::Float64, false),
        Field::new("Group2", DataType::Utf8, false),
        Field::new("Group2_Present", DataType::UInt32, false),
        Field::new("Group2_Total", DataType::UInt32, false),
        Field::new("Group2_Penetrance", DataType::Float64, false),
        Field::new("Bias_Direction", DataType::Utf8, false),
        Field::new("Bias", DataType::Float64, false),
        Field::new("P", DataType::Float64, false),
        Field::new("CorrectedP", DataType::Float64, false),
        Field::new("Bayes_Factor", DataType::Float64, false),
        Field::new("Posterior_SexLinked", DataType::Float64, false),
        Field::new("Strict_Call", DataType::Boolean, false),
        Field::new("Posterior_Call", DataType::Boolean, false),
        Field::new("Bayes_Factor_Call", DataType::Boolean, false),
        Field::new("Candidate_Class", DataType::Utf8, false),
    ]);

    let mut groups = GroupConfig {
        group1: params.group1.clone(),
        group2: params.group2.clone(),
    };
    popmap.resolve_groups(&mut groups)?;

    let g1_name = groups.group1.clone();
    let g2_name = groups.group2.clone();

    let total_g1 = popmap.get_count(&groups.group1);
    let total_g2 = popmap.get_count(&groups.group2);

    let n_individuals = source.header().n_individuals;
    let mask_g1 = GroupMask::from_columns(source.groups(), &groups.group1, n_individuals);
    let mask_g2 = GroupMask::from_columns(source.groups(), &groups.group2, n_individuals);

    let n_markers = source.count_markers()?;
    let corrected_threshold = if n_markers > 0 {
        params.signif_threshold as f64 / n_markers as f64
    } else {
        params.signif_threshold as f64
    };

    let cap = n_markers as usize;
    let mut id_b = StringBuilder::with_capacity(cap, cap * 16);
    let mut seq_b = StringBuilder::with_capacity(cap, cap * 64);
    let mut g1n_b = StringBuilder::with_capacity(cap, cap * 8);
    let mut g1p_b = UInt32Builder::with_capacity(cap);
    let mut g1t_b = UInt32Builder::with_capacity(cap);
    let mut g1pen_b = Float64Builder::with_capacity(cap);
    let mut g2n_b = StringBuilder::with_capacity(cap, cap * 8);
    let mut g2p_b = UInt32Builder::with_capacity(cap);
    let mut g2t_b = UInt32Builder::with_capacity(cap);
    let mut g2pen_b = Float64Builder::with_capacity(cap);
    let mut dir_b = StringBuilder::with_capacity(cap, cap * 16);
    let mut bias_b = Float64Builder::with_capacity(cap);
    let mut p_b = Float64Builder::with_capacity(cap);
    let mut cp_b = Float64Builder::with_capacity(cap);
    let mut bf_b = Float64Builder::with_capacity(cap);
    let mut post_b = Float64Builder::with_capacity(cap);
    let mut strict_b = BooleanBuilder::with_capacity(cap);
    let mut pcall_b = BooleanBuilder::with_capacity(cap);
    let mut bfcall_b = BooleanBuilder::with_capacity(cap);
    let mut class_b = StringBuilder::with_capacity(cap, cap * 16);

    source.for_each(|marker| {
        if marker.n_individuals == 0 {
            return;
        }

        let g1 = marker.presence.count_masked(&mask_g1);
        let g2 = marker.presence.count_masked(&mask_g2);

        let p = compute_p(TestMethod::ChiSquared, g1, g2, total_g1, total_g2);
        let p_corrected = stats::bonferroni_correct(p, n_markers);

        let strict_call = p < corrected_threshold;
        let bf = stats::bayes_factor_2x2(g1, g2, total_g1, total_g2);
        let posterior = stats::posterior_sex_linked(
            g1, g2, total_g1, total_g2,
            params.prior_probability, params.linked_probability,
        );

        let posterior_call = posterior > params.posterior_threshold;
        let bayes_factor_call = bf > params.bayes_factor_threshold;

        if !(strict_call || posterior_call || bayes_factor_call) {
            return;
        }

        let g1_pen = penetrance(g1, total_g1);
        let g2_pen = penetrance(g2, total_g2);
        let bias = stats::group_bias(g1, total_g1, g2, total_g2);
        let dir = bias_direction(&groups.group1, &groups.group2, g1_pen, g2_pen);
        let class = candidate_class(strict_call, posterior_call, bayes_factor_call);

        id_b.append_value(marker.id.clone());
        seq_b.append_value(&marker.sequence);
        g1n_b.append_value(&g1_name);
        g1p_b.append_value(g1);
        g1t_b.append_value(total_g1);
        g1pen_b.append_value(g1_pen);
        g2n_b.append_value(&g2_name);
        g2p_b.append_value(g2);
        g2t_b.append_value(total_g2);
        g2pen_b.append_value(g2_pen);
        dir_b.append_value(&dir);
        bias_b.append_value(bias);
        p_b.append_value(p);
        cp_b.append_value(p_corrected);
        bf_b.append_value(bf);
        post_b.append_value(posterior);
        strict_b.append_value(strict_call);
        pcall_b.append_value(posterior_call);
        bfcall_b.append_value(bayes_factor_call);
        class_b.append_value(class);
    })?;

    let batch = RecordBatch::try_new(
        std::sync::Arc::new(schema),
        vec![
            std::sync::Arc::new(id_b.finish()),
            std::sync::Arc::new(seq_b.finish()),
            std::sync::Arc::new(g1n_b.finish()),
            std::sync::Arc::new(g1p_b.finish()),
            std::sync::Arc::new(g1t_b.finish()),
            std::sync::Arc::new(g1pen_b.finish()),
            std::sync::Arc::new(g2n_b.finish()),
            std::sync::Arc::new(g2p_b.finish()),
            std::sync::Arc::new(g2t_b.finish()),
            std::sync::Arc::new(g2pen_b.finish()),
            std::sync::Arc::new(dir_b.finish()),
            std::sync::Arc::new(bias_b.finish()),
            std::sync::Arc::new(p_b.finish()),
            std::sync::Arc::new(cp_b.finish()),
            std::sync::Arc::new(bf_b.finish()),
            std::sync::Arc::new(post_b.finish()),
            std::sync::Arc::new(strict_b.finish()),
            std::sync::Arc::new(pcall_b.finish()),
            std::sync::Arc::new(bfcall_b.finish()),
            std::sync::Arc::new(class_b.finish()),
        ],
    )?;

    Ok(vec![batch])
}

#[cfg(all(test, feature = "arrow-output"))]
mod tests {
    use super::*;
    use std::io::Write;

    fn make_triage_fixture() -> (std::path::PathBuf, std::path::PathBuf) {
        let dir = std::env::temp_dir().join("rsx_arrow_triage_test");
        std::fs::create_dir_all(&dir).unwrap();

        let table = dir.join("markers.tsv");
        let mut f = std::fs::File::create(&table).unwrap();
        writeln!(f, "#Number of markers : 3").unwrap();
        write!(f, "id\tsequence").unwrap();
        for i in 1..=10 { write!(f, "\tm{i}").unwrap(); }
        for i in 1..=10 { write!(f, "\tf{i}").unwrap(); }
        writeln!(f).unwrap();

        write!(f, "0\tALL").unwrap();
        for _ in 0..20 { write!(f, "\t10").unwrap(); }
        writeln!(f).unwrap();

        write!(f, "1\tMONLY").unwrap();
        for _ in 0..10 { write!(f, "\t10").unwrap(); }
        for _ in 0..10 { write!(f, "\t0").unwrap(); }
        writeln!(f).unwrap();

        write!(f, "2\tFONLY").unwrap();
        for _ in 0..10 { write!(f, "\t0").unwrap(); }
        for _ in 0..10 { write!(f, "\t10").unwrap(); }
        writeln!(f).unwrap();

        let pop = dir.join("popmap.tsv");
        let mut f = std::fs::File::create(&pop).unwrap();
        for i in 1..=10 { writeln!(f, "m{i}\tM").unwrap(); }
        for i in 1..=10 { writeln!(f, "f{i}\tF").unwrap(); }

        (table, pop)
    }

    #[test]
    fn run_to_arrow_matches_file_based_triage() {
        let (table, pop) = make_triage_fixture();

        let params = TriageParams {
            markers_table_path: table.to_str().unwrap().to_string(),
            popmap_file_path: pop.to_str().unwrap().to_string(),
            output_file_path: "/tmp/discard.tsv".to_string(),
            min_depth: 1,
            signif_threshold: 0.05,
            posterior_threshold: 0.9,
            bayes_factor_threshold: 10.0,
            prior_probability: 0.01,
            linked_probability: 0.9,
            group1: "M".to_string(),
            group2: "F".to_string(),
        };

        let tsv_path = std::env::temp_dir().join("arrow_vs_file_triage.tsv");
        let mut file_params = params.clone();
        file_params.output_file_path = tsv_path.to_str().unwrap().to_string();
        run(&file_params).unwrap();

        let tsv = std::fs::read_to_string(&tsv_path).unwrap();
        let tsv_lines: Vec<&str> = tsv.lines().filter(|l| !l.starts_with('#') && !l.is_empty()).collect();
        let n_tsv = tsv_lines.len();

        let batches = run_to_arrow(&params).expect("run_to_arrow must succeed");
        assert!(!batches.is_empty(), "Arrow path produced at least one batch");
        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();

        let n_data_tsv = n_tsv.saturating_sub(1);
        assert_eq!(
            total_rows, n_data_tsv,
            "Arrow and file paths must agree on the number of qualifying markers"
        );
    }
}
