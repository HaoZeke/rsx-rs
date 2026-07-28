// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! `signif` command: extract markers significantly associated with a group.

use crate::bitset::GroupMask;
use crate::compute_backend::{
    AssociationCounts, PValueBackend, compute_chi_squared_batch_with_metrics,
};
use crate::markers_table::{MarkersTableStream, ParserConfig};
use crate::popmap::{GroupConfig, Popmap};
use crate::source::MarkerStream;
use crate::stats;
use crate::test_method::{CorrectionMethod, TestMethod, compute_p};
use std::io::Write;
use std::path::Path;

pub struct SignifParams {
    pub markers_table_path: String,
    pub popmap_file_path: String,
    pub output_file_path: String,
    pub min_depth: u16,
    pub signif_threshold: f32,
    pub correction: CorrectionMethod,
    pub test_method: TestMethod,
    pub output_fasta: bool,
    pub output_bayes: bool,
    pub group1: String,
    pub group2: String,
}

pub fn run(params: &SignifParams) -> Result<(), Box<dyn std::error::Error>> {
    run_with_backend(params, PValueBackend::Cpu)
}

pub fn run_with_backend(
    params: &SignifParams,
    backend: PValueBackend,
) -> Result<(), Box<dyn std::error::Error>> {
    let table_path = Path::new(&params.markers_table_path);
    let popmap = Popmap::from_file(Path::new(&params.popmap_file_path))?;
    let config = ParserConfig {
        store_sequence: true,
        store_depths: true,
        compute_groups: true,
        min_depth: params.min_depth,
    };
    let stream = MarkersTableStream::open(table_path, Some(&popmap), config)?;
    run_with_source_and_backend(&stream, &popmap, params, backend)
}

pub fn run_with_source<S: MarkerStream>(
    source: &S,
    popmap: &Popmap,
    params: &SignifParams,
) -> Result<(), Box<dyn std::error::Error>> {
    run_with_source_and_backend(source, popmap, params, PValueBackend::Cpu)
}

pub fn run_with_source_and_backend<S: MarkerStream>(
    source: &S,
    popmap: &Popmap,
    params: &SignifParams,
    backend: PValueBackend,
) -> Result<(), Box<dyn std::error::Error>> {
    // Reject invalid flag combinations *before* creating/truncating the output file.
    if matches!(params.correction, CorrectionMethod::Fdr) && params.output_fasta {
        return Err(
            "signif: --output-fasta is not supported with --correction fdr \
             (FDR needs a full p-value pass then a re-stream table write; \
             use table output or correction=bonferroni/none for FASTA)"
                .into(),
        );
    }

    let mut groups = GroupConfig {
        group1: params.group1.clone(),
        group2: params.group2.clone(),
    };
    popmap.resolve_groups(&mut groups)?;

    let total_g1 = popmap.get_count(&groups.group1);
    let total_g2 = popmap.get_count(&groups.group2);

    let header_columns = source.header().columns.clone();
    let n_individuals = source.header().n_individuals;
    let mask_g1 = GroupMask::from_columns(source.groups(), &groups.group1, n_individuals);
    let mask_g2 = GroupMask::from_columns(source.groups(), &groups.group2, n_individuals);

    let test_name = match params.test_method {
        TestMethod::ChiSquared => "chisq",
        TestMethod::Fisher => "fisher",
        TestMethod::GTest => "gtest",
    };
    let corr_name = match params.correction {
        CorrectionMethod::Bonferroni => "bonferroni",
        CorrectionMethod::Fdr => "fdr",
        CorrectionMethod::None => "none",
    };
    let threshold = params.signif_threshold as f64;

    let fasta_groups = vec![
        (groups.group1.clone(), &mask_g1),
        (groups.group2.clone(), &mask_g2),
    ];

    if matches!(backend, PValueBackend::Cuda) {
        if !matches!(params.test_method, TestMethod::ChiSquared) {
            return Err(format!(
                "signif: --backend cuda supports --test chisq; selected test is {test_name}"
            )
            .into());
        }
        return run_cuda(
            source,
            params,
            &header_columns,
            &mask_g1,
            &mask_g2,
            total_g1,
            total_g2,
            threshold,
            corr_name,
            test_name,
            &fasta_groups,
        );
    }

    // FDR: two full table passes only (p-values, then re-stream write). n_markers
    // for the header is p_values.len() — no separate count_markers pass.
    // Stores O(n_markers) p/q floats only. Documented in commands.org / README.
    if matches!(params.correction, CorrectionMethod::Fdr) {
        log::info!("signif FDR pass 1: collecting p-values");
        let mut p_values: Vec<f64> = Vec::new();
        source.for_each(|marker| {
            if marker.n_individuals > 0 {
                let g1 = marker.presence.count_masked(&mask_g1);
                let g2 = marker.presence.count_masked(&mask_g2);
                let p = compute_p(params.test_method, g1, g2, total_g1, total_g2);
                p_values.push(p);
            }
        })?;
        let n_markers = p_values.len() as u64;
        log::info!("signif FDR pass 1: {} markers", n_markers);

        let q_values = stats::benjamini_hochberg(&p_values);

        log::info!("signif FDR pass 2: filtering and writing");
        let mut output = std::io::BufWriter::new(std::fs::File::create(&params.output_file_path)?);
        writeln!(
            output,
            "#source:rsx-signif;min_depth:{};signif_threshold:{};correction:{};test:{};n_markers:{}",
            params.min_depth, params.signif_threshold, corr_name, test_name, n_markers
        )?;
        if params.output_bayes {
            writeln!(
                output,
                "{}\tBayes_Factor\tPosterior_SexLinked",
                header_columns.join("\t")
            )?;
        } else {
            writeln!(output, "{}", header_columns.join("\t"))?;
        }

        let mut idx = 0usize;
        let mut write_err: Option<std::io::Error> = None;
        source.for_each(|marker| {
            if write_err.is_some() || marker.n_individuals == 0 {
                return;
            }
            let q = q_values[idx];
            idx += 1;
            if q >= threshold {
                return;
            }
            let g1 = marker.presence.count_masked(&mask_g1);
            let g2 = marker.presence.count_masked(&mask_g2);
            let result = if params.output_bayes {
                write_marker_bayes_row(&mut output, marker, g1, g2, total_g1, total_g2)
            } else {
                marker.write_as_table(&mut output)
            };
            if let Err(e) = result {
                write_err = Some(e);
            }
        })?;
        if let Some(e) = write_err {
            return Err(e.into());
        }
        return Ok(());
    }

    // Bonferroni / none: count then single filter/write pass.
    log::info!("signif pass 1: counting markers");
    let n_markers = source.count_markers()?;
    log::info!("signif pass 1: {} markers", n_markers);

    let corrected_threshold = match params.correction {
        CorrectionMethod::Bonferroni => {
            if n_markers == 0 {
                threshold
            } else {
                threshold / n_markers as f64
            }
        }
        CorrectionMethod::None => threshold,
        CorrectionMethod::Fdr => unreachable!("FDR handled above"),
    };
    let effective_n_markers = match params.correction {
        CorrectionMethod::Bonferroni => n_markers.max(1),
        CorrectionMethod::None => 1,
        CorrectionMethod::Fdr => unreachable!("FDR handled above"),
    };

    log::info!("signif pass 2: filtering and writing");
    let mut output = std::io::BufWriter::new(std::fs::File::create(&params.output_file_path)?);

    if !params.output_fasta {
        writeln!(
            output,
            "#source:rsx-signif;min_depth:{};signif_threshold:{};correction:{};test:{};n_markers:{}",
            params.min_depth, params.signif_threshold, corr_name, test_name, n_markers
        )?;

        if params.output_bayes {
            writeln!(
                output,
                "{}\tBayes_Factor\tPosterior_SexLinked",
                header_columns.join("\t")
            )?;
        } else {
            writeln!(output, "{}", header_columns.join("\t"))?;
        }
    }

    let mut write_err: Option<std::io::Error> = None;
    source.for_each(|marker| {
        if write_err.is_some() || marker.n_individuals == 0 {
            return;
        }
        let g1 = marker.presence.count_masked(&mask_g1);
        let g2 = marker.presence.count_masked(&mask_g2);
        let p = compute_p(params.test_method, g1, g2, total_g1, total_g2);

        // Full f64 compare (do not cast to f32 — lossy near thresholds).
        if p >= corrected_threshold {
            return;
        }
        let p_corr = stats::bonferroni_correct(p, effective_n_markers);

        let result = if params.output_fasta {
            let mut m = marker.clone();
            m.p = p;
            m.p_corrected = p_corr;
            m.write_as_fasta_bitset(&mut output, params.min_depth as u32, &fasta_groups)
        } else if params.output_bayes {
            write_marker_bayes_row(&mut output, marker, g1, g2, total_g1, total_g2)
        } else {
            marker.write_as_table(&mut output)
        };
        if let Err(e) = result {
            write_err = Some(e);
        }
    })?;
    if let Some(e) = write_err {
        return Err(e.into());
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn run_cuda<S: MarkerStream>(
    source: &S,
    params: &SignifParams,
    header_columns: &[String],
    mask_g1: &GroupMask,
    mask_g2: &GroupMask,
    total_g1: u32,
    total_g2: u32,
    threshold: f64,
    corr_name: &str,
    test_name: &str,
    fasta_groups: &[(String, &GroupMask)],
) -> Result<(), Box<dyn std::error::Error>> {
    log::info!("signif CUDA pass 1: collecting marker counts");
    let mut counts = Vec::new();
    source.for_each(|marker| {
        if marker.n_individuals > 0 {
            counts.push(AssociationCounts {
                group1: marker.presence.count_masked(mask_g1),
                group2: marker.presence.count_masked(mask_g2),
            });
        }
    })?;
    let n_markers = counts.len() as u64;
    let cuda_result = compute_chi_squared_batch_with_metrics(
        PValueBackend::Cuda,
        &counts,
        total_g1,
        total_g2,
    )?;
    let p_values = cuda_result.p_values.try_as_slice()?;

    let corrected: Option<Vec<f64>> = match params.correction {
        CorrectionMethod::Fdr => Some(stats::benjamini_hochberg(&p_values)),
        CorrectionMethod::Bonferroni | CorrectionMethod::None => None,
    };
    let corrected_threshold = match params.correction {
        CorrectionMethod::Bonferroni if n_markers > 0 => threshold / n_markers as f64,
        CorrectionMethod::Bonferroni | CorrectionMethod::None | CorrectionMethod::Fdr => threshold,
    };
    let effective_n_markers = if matches!(params.correction, CorrectionMethod::Bonferroni) {
        n_markers.max(1)
    } else {
        1
    };

    log::info!("signif CUDA pass 2: filtering and writing");
    let mut output = std::io::BufWriter::new(std::fs::File::create(&params.output_file_path)?);
    if !params.output_fasta {
        writeln!(
            output,
            "#source:rsx-signif;min_depth:{};signif_threshold:{};correction:{};test:{};n_markers:{}",
            params.min_depth, params.signif_threshold, corr_name, test_name, n_markers
        )?;
        if params.output_bayes {
            writeln!(
                output,
                "{}\tBayes_Factor\tPosterior_SexLinked",
                header_columns.join("\t")
            )?;
        } else {
            writeln!(output, "{}", header_columns.join("\t"))?;
        }
    }

    let mut index = 0usize;
    let mut write_err: Option<std::io::Error> = None;
    source.for_each(|marker| {
        if write_err.is_some() || marker.n_individuals == 0 {
            return;
        }
        let p = p_values[index];
        let passes = match &corrected {
            Some(q_values) => q_values[index] < threshold,
            None => p < corrected_threshold,
        };
        index += 1;
        if !passes {
            return;
        }

        let association = counts[index - 1];
        let p_corrected = match &corrected {
            Some(q_values) => q_values[index - 1],
            None => stats::bonferroni_correct(p, effective_n_markers),
        };
        let result = if params.output_fasta {
            let mut marker = marker.clone();
            marker.p = p;
            marker.p_corrected = p_corrected;
            marker.write_as_fasta_bitset(&mut output, params.min_depth as u32, fasta_groups)
        } else if params.output_bayes {
            write_marker_bayes_row(
                &mut output,
                marker,
                association.group1,
                association.group2,
                total_g1,
                total_g2,
            )
        } else {
            marker.write_as_table(&mut output)
        };
        if let Err(error) = result {
            write_err = Some(error);
        }
    })?;
    if let Some(error) = write_err {
        return Err(error.into());
    }
    Ok(())
}

fn write_marker_bayes_row<W: Write>(
    output: &mut W,
    marker: &crate::marker::Marker,
    g1: u32,
    g2: u32,
    total_g1: u32,
    total_g2: u32,
) -> std::io::Result<()> {
    let bf = stats::bayes_factor_2x2(g1, g2, total_g1, total_g2);
    let post = stats::posterior_sex_linked(g1, g2, total_g1, total_g2, 0.01, 0.9);
    write!(output, "{}\t{}", marker.id, marker.sequence)?;
    for &d in &marker.individual_depths {
        write!(output, "\t{d}")?;
    }
    writeln!(output, "\t{:.4}\t{:.4}", bf, post)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_method::{CorrectionMethod, TestMethod};
    use std::io::Write;

    fn write_fixture(dir: &std::path::Path) -> (std::path::PathBuf, std::path::PathBuf) {
        let table = dir.join("markers.tsv");
        let mut f = std::fs::File::create(&table).unwrap();
        // 8M + 8F so Yates χ² has power after Bonferroni (n_markers small).
        write!(f, "#Number of markers : 2\n").unwrap();
        write!(f, "id\tsequence").unwrap();
        for i in 0..8 {
            write!(f, "\tm{i}").unwrap();
        }
        for i in 0..8 {
            write!(f, "\tf{i}").unwrap();
        }
        writeln!(f).unwrap();
        // Strong male-only marker (present in all males, absent in females)
        write!(f, "0\tAAAAAAAA").unwrap();
        for _ in 0..8 {
            write!(f, "\t10").unwrap();
        }
        for _ in 0..8 {
            write!(f, "\t0").unwrap();
        }
        writeln!(f).unwrap();
        // Balanced
        write!(f, "1\tCCCCCCCC").unwrap();
        for _ in 0..16 {
            write!(f, "\t5").unwrap();
        }
        writeln!(f).unwrap();
        let pop = dir.join("popmap.tsv");
        let mut p = std::fs::File::create(&pop).unwrap();
        for i in 0..8 {
            writeln!(p, "m{i}\tM").unwrap();
        }
        for i in 0..8 {
            writeln!(p, "f{i}\tF").unwrap();
        }
        (table, pop)
    }

    #[test]
    fn signif_f64_threshold_emits_strong_marker() {
        let dir = tempfile::tempdir().unwrap();
        let (table, pop) = write_fixture(dir.path());
        let out = dir.path().join("signif.tsv");
        run(&SignifParams {
            markers_table_path: table.to_str().unwrap().to_string(),
            popmap_file_path: pop.to_str().unwrap().to_string(),
            output_file_path: out.to_str().unwrap().to_string(),
            min_depth: 1,
            signif_threshold: 0.05,
            correction: CorrectionMethod::Bonferroni,
            test_method: TestMethod::ChiSquared,
            output_fasta: false,
            output_bayes: false,
            group1: "M".into(),
            group2: "F".into(),
        })
        .unwrap();
        let body = std::fs::read_to_string(&out).unwrap();
        assert!(
            body.contains("AAAAAAAA"),
            "strong association must pass f64 gate: {body}"
        );
        assert!(
            !body
                .lines()
                .any(|l| l.contains("CCCCCCCC") && !l.starts_with('#')),
            "balanced marker must not be emitted: {body}"
        );
    }

    #[test]
    fn signif_fdr_plus_fasta_is_rejected_without_truncating_output() {
        let dir = tempfile::tempdir().unwrap();
        let (table, pop) = write_fixture(dir.path());
        let out = dir.path().join("preexisting.fa");
        std::fs::write(&out, b">keep_me\nACGT\n").unwrap();
        let err = run(&SignifParams {
            markers_table_path: table.to_str().unwrap().to_string(),
            popmap_file_path: pop.to_str().unwrap().to_string(),
            output_file_path: out.to_str().unwrap().to_string(),
            min_depth: 1,
            signif_threshold: 0.05,
            correction: CorrectionMethod::Fdr,
            test_method: TestMethod::ChiSquared,
            output_fasta: true,
            output_bayes: false,
            group1: "M".into(),
            group2: "F".into(),
        });
        assert!(err.is_err(), "FASTA+FDR must error");
        let msg = format!("{}", err.unwrap_err());
        assert!(
            msg.contains("fdr") || msg.contains("FASTA") || msg.contains("fasta"),
            "error should mention fdr/fasta: {msg}"
        );
        // Must not truncate an existing path on the rejected invocation.
        let preserved = std::fs::read_to_string(&out).unwrap();
        assert_eq!(
            preserved, ">keep_me\nACGT\n",
            "rejected FASTA+FDR must not wipe an existing output file"
        );
    }

    #[test]
    fn signif_fdr_two_pass_writes_table() {
        let dir = tempfile::tempdir().unwrap();
        let (table, pop) = write_fixture(dir.path());
        let out = dir.path().join("signif_fdr.tsv");
        run(&SignifParams {
            markers_table_path: table.to_str().unwrap().to_string(),
            popmap_file_path: pop.to_str().unwrap().to_string(),
            output_file_path: out.to_str().unwrap().to_string(),
            min_depth: 1,
            signif_threshold: 0.05,
            correction: CorrectionMethod::Fdr,
            test_method: TestMethod::ChiSquared,
            output_fasta: false,
            output_bayes: false,
            group1: "M".into(),
            group2: "F".into(),
        })
        .unwrap();
        let body = std::fs::read_to_string(&out).unwrap();
        assert!(
            body.contains("#source:rsx-signif"),
            "header present: {body}"
        );
        assert!(
            body.contains("AAAAAAAA"),
            "FDR should keep strong marker: {body}"
        );
    }

    #[cfg(feature = "cuda")]
    #[test]
    fn signif_cuda_matches_cpu_output() {
        let dir = tempfile::tempdir().unwrap();
        let (table, pop) = write_fixture(dir.path());
        let cpu_out = dir.path().join("signif_cpu.tsv");
        let cuda_out = dir.path().join("signif_cuda.tsv");
        let params = |output: &std::path::Path| SignifParams {
            markers_table_path: table.to_str().unwrap().to_string(),
            popmap_file_path: pop.to_str().unwrap().to_string(),
            output_file_path: output.to_str().unwrap().to_string(),
            min_depth: 1,
            signif_threshold: 0.05,
            correction: CorrectionMethod::Fdr,
            test_method: TestMethod::ChiSquared,
            output_fasta: false,
            output_bayes: false,
            group1: "M".into(),
            group2: "F".into(),
        };

        run_with_backend(
            &params(&cpu_out),
            crate::compute_backend::PValueBackend::Cpu,
        )
        .unwrap();
        run_with_backend(
            &params(&cuda_out),
            crate::compute_backend::PValueBackend::Cuda,
        )
        .unwrap();

        assert_eq!(
            std::fs::read_to_string(cpu_out).unwrap(),
            std::fs::read_to_string(cuda_out).unwrap()
        );
    }
}
