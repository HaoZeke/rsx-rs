// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! `triage` command: marker-level biological candidate ranking.
//!
//! The command keeps RADSex-style strict testing and Bayesian marker evidence
//! in one bounded-memory pass over the marker table.
//!
//! Every marker needs a p-value, a Bayes factor, and a posterior, so the
//! device backend evaluates all three for the whole table before the writing
//! pass. Counts cost eight bytes per marker; marker text is never retained.

use crate::bitset::GroupMask;
use crate::compute_backend::{AssociationCounts, PValueBackend, compute_bayes_evidence_batch};
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
    pub bayes_model: stats::DirectionalModel,
    pub group1: String,
    pub group2: String,
    /// Where the per-marker p-value and Bayesian evidence are evaluated.
    pub backend: PValueBackend,
}

/// Evidence for every marker, in table order, when a device evaluates it.
struct BatchedEvidence {
    p_values: Vec<f64>,
    bayes_factors: Vec<f64>,
    posteriors: Vec<f64>,
}

impl BatchedEvidence {
    fn collect<S: MarkerStream>(
        source: &S,
        mask_g1: &GroupMask,
        mask_g2: &GroupMask,
        total_g1: u32,
        total_g2: u32,
        params: &TriageParams,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        log::info!(
            "triage {} pass 1: collecting marker counts",
            params.backend.name()
        );
        let mut counts = Vec::new();
        source.for_each(|marker| {
            if marker.n_individuals > 0 {
                counts.push(AssociationCounts {
                    group1: marker.presence.count_masked(mask_g1),
                    group2: marker.presence.count_masked(mask_g2),
                });
            }
        })?;

        let p_values = crate::compute_backend::compute_p_batch(
            params.backend,
            TestMethod::ChiSquared,
            &counts,
            total_g1,
            total_g2,
        )?;
        let (bayes_factors, posteriors) = compute_bayes_evidence_batch(
            params.backend,
            &counts,
            total_g1,
            total_g2,
            &params.bayes_model,
        )?;
        Ok(Self {
            p_values,
            bayes_factors,
            posteriors,
        })
    }
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

#[derive(Clone)]
struct TriageRow {
    id: String,
    sequence: String,
    group1: String,
    group1_present: u32,
    group1_total: u32,
    group1_penetrance: f64,
    group2: String,
    group2_present: u32,
    group2_total: u32,
    group2_penetrance: f64,
    bias_direction: String,
    bias: f64,
    p: f64,
    corrected_p: f64,
    bayes_factor: f64,
    posterior: f64,
    strict_call: bool,
    posterior_call: bool,
    bayes_factor_call: bool,
    class: &'static str,
}

fn scores_beat(posterior: f64, bayes_factor: f64, p: f64, other: &TriageRow) -> bool {
    posterior
        .total_cmp(&other.posterior)
        .then_with(|| bayes_factor.total_cmp(&other.bayes_factor))
        .then_with(|| other.p.total_cmp(&p))
        .is_gt()
}

struct TriageLattice {
    stride: usize,
    p: Vec<f64>,
    bayes_factor: Vec<f64>,
    posterior: Vec<f64>,
    total_g1: u32,
    total_g2: u32,
}

impl TriageLattice {
    fn build(model: &stats::DirectionalModel, total_g1: u32, total_g2: u32) -> Self {
        let cells = (total_g1 as usize + 1).saturating_mul(total_g2 as usize + 1);
        if cells > 250_000 {
            return Self {
                stride: 0,
                p: Vec::new(),
                bayes_factor: Vec::new(),
                posterior: Vec::new(),
                total_g1,
                total_g2,
            };
        }
        let mut p = Vec::with_capacity(cells);
        let mut bayes_factor = Vec::with_capacity(cells);
        let mut posterior = Vec::with_capacity(cells);
        for g1 in 0..=total_g1 {
            for g2 in 0..=total_g2 {
                p.push(compute_p(
                    TestMethod::ChiSquared,
                    g1,
                    g2,
                    total_g1,
                    total_g2,
                ));
                bayes_factor.push(stats::bayes_factor_2x2_with_validated_model(
                    g1,
                    g2,
                    total_g1,
                    total_g2,
                    &model.bayes_factor,
                ));
                posterior.push(stats::posterior_sex_linked_with_model(
                    g1, g2, total_g1, total_g2, model,
                ));
            }
        }
        Self {
            stride: total_g2 as usize + 1,
            p,
            bayes_factor,
            posterior,
            total_g1,
            total_g2,
        }
    }

    fn get(&self, g1: u32, g2: u32, model: &stats::DirectionalModel) -> (f64, f64, f64) {
        if self.stride == 0 || g1 > self.total_g1 || g2 > self.total_g2 {
            return (
                compute_p(TestMethod::ChiSquared, g1, g2, self.total_g1, self.total_g2),
                stats::bayes_factor_2x2_with_validated_model(
                    g1,
                    g2,
                    self.total_g1,
                    self.total_g2,
                    &model.bayes_factor,
                ),
                stats::posterior_sex_linked_with_model(g1, g2, self.total_g1, self.total_g2, model),
            );
        }
        let index = g1 as usize * self.stride + g2 as usize;
        (
            self.p[index],
            self.bayes_factor[index],
            self.posterior[index],
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn build_triage_row(
    id: &str,
    sequence: &str,
    group1: &str,
    group2: &str,
    g1: u32,
    g2: u32,
    total_g1: u32,
    total_g2: u32,
    p: f64,
    corrected_p: f64,
    bayes_factor: f64,
    posterior: f64,
    strict_call: bool,
    posterior_call: bool,
    bayes_factor_call: bool,
) -> TriageRow {
    let g1_penetrance = penetrance(g1, total_g1);
    let g2_penetrance = penetrance(g2, total_g2);
    TriageRow {
        id: id.to_string(),
        sequence: sequence.to_string(),
        group1: group1.to_string(),
        group1_present: g1,
        group1_total: total_g1,
        group1_penetrance: g1_penetrance,
        group2: group2.to_string(),
        group2_present: g2,
        group2_total: total_g2,
        group2_penetrance: g2_penetrance,
        bias_direction: bias_direction(group1, group2, g1_penetrance, g2_penetrance),
        bias: stats::group_bias(g1, total_g1, g2, total_g2),
        p,
        corrected_p,
        bayes_factor,
        posterior,
        strict_call,
        posterior_call,
        bayes_factor_call,
        class: candidate_class(strict_call, posterior_call, bayes_factor_call),
    }
}

fn write_row<W: Write>(output: &mut W, row: &TriageRow) -> std::io::Result<()> {
    writeln!(
        output,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        row.id,
        row.sequence,
        row.group1,
        row.group1_present,
        row.group1_total,
        Cg(row.group1_penetrance),
        row.group2,
        row.group2_present,
        row.group2_total,
        Cg(row.group2_penetrance),
        row.bias_direction,
        Cg(row.bias),
        Cg(row.p),
        Cg(row.corrected_p),
        Cg(row.bayes_factor),
        Cg(row.posterior),
        row.strict_call,
        row.posterior_call,
        row.bayes_factor_call,
        row.class
    )
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
    params.bayes_model.bayes_factor.validate()?;
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
        Cg(params.bayes_model.linkage_prior),
        Cg(params.bayes_model.linked_prevalence),
        n_markers
    )?;
    writeln!(
        output,
        "id\tsequence\tGroup1\tGroup1_Present\tGroup1_Total\tGroup1_Penetrance\tGroup2\tGroup2_Present\tGroup2_Total\tGroup2_Penetrance\tBias_Direction\tBias\tP\tCorrectedP\tBayes_Factor\tPosterior_SexLinked\tStrict_Call\tPosterior_Call\tBayes_Factor_Call\tCandidate_Class"
    )?;

    let batched = match params.backend {
        PValueBackend::Cpu => None,
        PValueBackend::Cuda => Some(BatchedEvidence::collect(
            source, &mask_g1, &mask_g2, total_g1, total_g2, params,
        )?),
    };
    if batched.is_some() {
        log::info!(
            "triage {} pass 2: ranking and writing",
            params.backend.name()
        );
    }

    let lattice = TriageLattice::build(&params.bayes_model, total_g1, total_g2);
    let mut write_err: Option<std::io::Error> = None;
    let mut called_rows = 0_u64;
    let mut best_exploratory: Option<TriageRow> = None;
    let mut evidence_index = 0usize;
    source.for_each(|marker| {
        if write_err.is_some() || marker.n_individuals == 0 {
            return;
        }

        let g1 = marker.presence.count_masked(&mask_g1);
        let g2 = marker.presence.count_masked(&mask_g2);
        let (p, bf, posterior) = match &batched {
            Some(evidence) => {
                let index = evidence_index;
                evidence_index += 1;
                (
                    evidence.p_values[index],
                    evidence.bayes_factors[index],
                    evidence.posteriors[index],
                )
            }
            None => lattice.get(g1, g2, &params.bayes_model),
        };
        let p_corrected = stats::bonferroni_correct(p, n_markers);
        let strict_call = p < strict_threshold;
        let posterior_call = posterior > params.posterior_threshold;
        let bayes_factor_call = bf > params.bayes_factor_threshold;
        let called = strict_call || posterior_call || bayes_factor_call;
        if !called {
            let keep = best_exploratory
                .as_ref()
                .is_none_or(|best| scores_beat(posterior, bf, p, best));
            if keep {
                best_exploratory = Some(build_triage_row(
                    &marker.id,
                    &marker.sequence,
                    &groups.group1,
                    &groups.group2,
                    g1,
                    g2,
                    total_g1,
                    total_g2,
                    p,
                    p_corrected,
                    bf,
                    posterior,
                    strict_call,
                    posterior_call,
                    bayes_factor_call,
                ));
            }
            return;
        }

        called_rows += 1;
        let row = build_triage_row(
            &marker.id,
            &marker.sequence,
            &groups.group1,
            &groups.group2,
            g1,
            g2,
            total_g1,
            total_g2,
            p,
            p_corrected,
            bf,
            posterior,
            strict_call,
            posterior_call,
            bayes_factor_call,
        );
        if let Err(e) = write_row(&mut output, &row) {
            write_err = Some(e);
        }
    })?;
    if let Some(e) = write_err {
        return Err(e.into());
    }
    if called_rows == 0
        && let Some(row) = best_exploratory
    {
        write_row(&mut output, &row)?;
    }

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
    use arrow::array::builder::{BooleanBuilder, Float64Builder, StringBuilder, UInt32Builder};
    use arrow::datatypes::{DataType, Field, Schema};

    params.bayes_model.bayes_factor.validate()?;

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

    let mut append_row = |row: &TriageRow| {
        id_b.append_value(&row.id);
        seq_b.append_value(&row.sequence);
        g1n_b.append_value(&row.group1);
        g1p_b.append_value(row.group1_present);
        g1t_b.append_value(row.group1_total);
        g1pen_b.append_value(row.group1_penetrance);
        g2n_b.append_value(&row.group2);
        g2p_b.append_value(row.group2_present);
        g2t_b.append_value(row.group2_total);
        g2pen_b.append_value(row.group2_penetrance);
        dir_b.append_value(&row.bias_direction);
        bias_b.append_value(row.bias);
        p_b.append_value(row.p);
        cp_b.append_value(row.corrected_p);
        bf_b.append_value(row.bayes_factor);
        post_b.append_value(row.posterior);
        strict_b.append_value(row.strict_call);
        pcall_b.append_value(row.posterior_call);
        bfcall_b.append_value(row.bayes_factor_call);
        class_b.append_value(row.class);
    };

    let lattice = TriageLattice::build(&params.bayes_model, total_g1, total_g2);
    let mut called_rows = 0_u64;
    let mut best_exploratory: Option<TriageRow> = None;

    source.for_each(|marker| {
        if marker.n_individuals == 0 {
            return;
        }

        let g1 = marker.presence.count_masked(&mask_g1);
        let g2 = marker.presence.count_masked(&mask_g2);
        let (p, bf, posterior) = lattice.get(g1, g2, &params.bayes_model);
        let p_corrected = stats::bonferroni_correct(p, n_markers);
        let strict_call = p < corrected_threshold;
        let posterior_call = posterior > params.posterior_threshold;
        let bayes_factor_call = bf > params.bayes_factor_threshold;
        let called = strict_call || posterior_call || bayes_factor_call;
        if !called {
            let keep = best_exploratory
                .as_ref()
                .is_none_or(|best| scores_beat(posterior, bf, p, best));
            if keep {
                best_exploratory = Some(build_triage_row(
                    &marker.id,
                    &marker.sequence,
                    &g1_name,
                    &g2_name,
                    g1,
                    g2,
                    total_g1,
                    total_g2,
                    p,
                    p_corrected,
                    bf,
                    posterior,
                    strict_call,
                    posterior_call,
                    bayes_factor_call,
                ));
            }
            return;
        }

        called_rows += 1;
        let row = build_triage_row(
            &marker.id,
            &marker.sequence,
            &g1_name,
            &g2_name,
            g1,
            g2,
            total_g1,
            total_g2,
            p,
            p_corrected,
            bf,
            posterior,
            strict_call,
            posterior_call,
            bayes_factor_call,
        );
        append_row(&row);
    })?;

    if called_rows == 0
        && let Some(row) = best_exploratory
    {
        append_row(&row);
    }

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

#[cfg(test)]
mod tests {
    #[cfg(any(feature = "cuda", feature = "arrow-output"))]
    use super::*;
    #[cfg(any(feature = "cuda", feature = "arrow-output"))]
    use std::io::Write;

    #[cfg(any(feature = "cuda", feature = "arrow-output"))]
    fn make_triage_fixture() -> (std::path::PathBuf, std::path::PathBuf) {
        let dir = tempfile::tempdir().unwrap().keep();
        std::fs::create_dir_all(&dir).unwrap();

        let table = dir.join("markers.tsv");
        let mut f = std::fs::File::create(&table).unwrap();
        writeln!(f, "#Number of markers : 3").unwrap();
        write!(f, "id\tsequence").unwrap();
        for i in 1..=10 {
            write!(f, "\tm{i}").unwrap();
        }
        for i in 1..=10 {
            write!(f, "\tf{i}").unwrap();
        }
        writeln!(f).unwrap();

        write!(f, "0\tALL").unwrap();
        for _ in 0..20 {
            write!(f, "\t10").unwrap();
        }
        writeln!(f).unwrap();

        write!(f, "1\tMONLY").unwrap();
        for _ in 0..10 {
            write!(f, "\t10").unwrap();
        }
        for _ in 0..10 {
            write!(f, "\t0").unwrap();
        }
        writeln!(f).unwrap();

        write!(f, "2\tFONLY").unwrap();
        for _ in 0..10 {
            write!(f, "\t0").unwrap();
        }
        for _ in 0..10 {
            write!(f, "\t10").unwrap();
        }
        writeln!(f).unwrap();

        let pop = dir.join("popmap.tsv");
        let mut f = std::fs::File::create(&pop).unwrap();
        for i in 1..=10 {
            writeln!(f, "m{i}\tM").unwrap();
        }
        for i in 1..=10 {
            writeln!(f, "f{i}\tF").unwrap();
        }

        (table, pop)
    }

    #[cfg(feature = "cuda")]
    #[test]
    fn triage_cuda_matches_cpu_output() {
        let (table, pop) = make_triage_fixture();
        let directory = table.parent().expect("fixture table has parent");
        let cpu_out = directory.join("triage_cpu.tsv");
        let cuda_out = directory.join("triage_cuda.tsv");

        let params = |output: &std::path::Path, backend| TriageParams {
            markers_table_path: table.to_str().unwrap().to_string(),
            popmap_file_path: pop.to_str().unwrap().to_string(),
            output_file_path: output.to_str().unwrap().to_string(),
            min_depth: 1,
            signif_threshold: 0.05,
            posterior_threshold: 0.9,
            bayes_factor_threshold: 10.0,
            bayes_model: stats::DirectionalModel::directional_screening_v1(),
            group1: "M".to_string(),
            group2: "F".to_string(),
            backend,
        };

        run(&params(&cpu_out, PValueBackend::Cpu)).unwrap();
        run(&params(&cuda_out, PValueBackend::Cuda)).unwrap();

        let cpu = std::fs::read_to_string(&cpu_out).unwrap();
        assert!(
            cpu.lines().filter(|line| !line.starts_with('#')).count() > 1,
            "fixture must emit at least one ranked row"
        );
        assert_eq!(cpu, std::fs::read_to_string(&cuda_out).unwrap());
    }

    #[cfg(feature = "arrow-output")]
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
            bayes_model: stats::DirectionalModel::directional_screening_v1(),
            group1: "M".to_string(),
            group2: "F".to_string(),
            backend: PValueBackend::Cpu,
        };

        let tsv_path = table
            .parent()
            .expect("fixture table has parent")
            .join("arrow_vs_file_triage.tsv");
        let mut file_params = params.clone();
        file_params.output_file_path = tsv_path.to_str().unwrap().to_string();
        run(&file_params).unwrap();

        let tsv = std::fs::read_to_string(&tsv_path).unwrap();
        let tsv_lines: Vec<&str> = tsv
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
            .collect();
        let n_tsv = tsv_lines.len();

        let batches = run_to_arrow(&params).expect("run_to_arrow must succeed");
        assert!(
            !batches.is_empty(),
            "Arrow path produced at least one batch"
        );
        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();

        let n_data_tsv = n_tsv.saturating_sub(1);
        assert_eq!(
            total_rows, n_data_tsv,
            "Arrow and file paths must agree on the number of qualifying markers"
        );
    }
}
