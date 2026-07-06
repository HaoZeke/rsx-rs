// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! `signif` command: extract markers significantly associated with a group.

use crate::bitset::GroupMask;
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
    let table_path = Path::new(&params.markers_table_path);
    let popmap = Popmap::from_file(Path::new(&params.popmap_file_path))?;
    let config = ParserConfig {
        store_sequence: true,
        store_depths: true,
        compute_groups: true,
        min_depth: params.min_depth,
    };
    let stream = MarkersTableStream::open(table_path, Some(&popmap), config)?;
    run_with_source(&stream, &popmap, params)
}

pub fn run_with_source<S: MarkerStream>(
    source: &S,
    popmap: &Popmap,
    params: &SignifParams,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut groups = GroupConfig {
        group1: params.group1.clone(),
        group2: params.group2.clone(),
    };
    popmap.resolve_groups(&mut groups)?;

    let total_g1 = popmap.get_count(&groups.group1);
    let total_g2 = popmap.get_count(&groups.group2);

    log::info!("signif pass 1: counting markers");
    let n_markers = source.count_markers()?;
    log::info!("signif pass 1: {} markers", n_markers);

    let threshold = params.signif_threshold as f64;
    let corrected_threshold = match params.correction {
        CorrectionMethod::Bonferroni => threshold / n_markers as f64,
        CorrectionMethod::None => threshold,
        CorrectionMethod::Fdr => threshold,
    };
    let effective_n_markers = match params.correction {
        CorrectionMethod::Bonferroni => n_markers,
        CorrectionMethod::None | CorrectionMethod::Fdr => 1,
    };

    log::info!("signif pass 2: filtering and writing");
    let header_columns = source.header().columns.clone();
    let n_individuals = source.header().n_individuals;
    let mask_g1 = GroupMask::from_columns(source.groups(), &groups.group1, n_individuals);
    let mask_g2 = GroupMask::from_columns(source.groups(), &groups.group2, n_individuals);

    let mut output = std::io::BufWriter::new(std::fs::File::create(&params.output_file_path)?);

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

    let fasta_groups = vec![
        (groups.group1.clone(), &mask_g1),
        (groups.group2.clone(), &mask_g2),
    ];

    // FDR needs all p-values before writing (BH step-up). Pass 2 re-streams the
    // table so we only store O(n_markers) p-values, not full sequences/depths.
    // This is still O(n_markers) memory for p/q, not O(n_individuals). Documented
    // in commands.org / README. FASTA + FDR is rejected (would need full materialize
    // or a third pass with presence metadata in the header).
    if matches!(params.correction, CorrectionMethod::Fdr) {
        if params.output_fasta {
            return Err(
                "signif: --output-fasta is not supported with --correction fdr \
                 (FDR needs a full p-value pass then a re-stream table write; \
                 use table output or correction=bonferroni/none for FASTA)"
                    .into(),
            );
        }

        let mut p_values: Vec<f64> = Vec::new();
        source.for_each(|marker| {
            if marker.n_individuals > 0 {
                let g1 = marker.presence.count_masked(&mask_g1);
                let g2 = marker.presence.count_masked(&mask_g2);
                let p = compute_p(params.test_method, g1, g2, total_g1, total_g2);
                p_values.push(p);
            }
        })?;

        let q_values = stats::benjamini_hochberg(&p_values);
        let mut idx = 0usize;
        source.for_each(|marker| {
            if marker.n_individuals == 0 {
                return;
            }
            let q = q_values[idx];
            idx += 1;
            if q < threshold {
                let g1 = marker.presence.count_masked(&mask_g1);
                let g2 = marker.presence.count_masked(&mask_g2);
                if params.output_bayes {
                    let bf = stats::bayes_factor_2x2(g1, g2, total_g1, total_g2);
                    let post =
                        stats::posterior_sex_linked(g1, g2, total_g1, total_g2, 0.01, 0.9);
                    write!(output, "{}\t{}", marker.id, marker.sequence).ok();
                    for &d in &marker.individual_depths {
                        write!(output, "\t{d}").ok();
                    }
                    let _ = writeln!(output, "\t{:.4}\t{:.4}", bf, post);
                } else {
                    let _ = marker.write_as_table(&mut output);
                }
            }
        })?;
    } else {
        source.for_each(|marker| {
            if marker.n_individuals > 0 {
                let g1 = marker.presence.count_masked(&mask_g1);
                let g2 = marker.presence.count_masked(&mask_g2);
                let p = compute_p(params.test_method, g1, g2, total_g1, total_g2);

                // Full f64 compare (do not cast to f32 — lossy near thresholds).
                if p < corrected_threshold {
                    let p_corr = stats::bonferroni_correct(p, effective_n_markers);

                    if params.output_fasta {
                        let mut m = marker.clone();
                        m.p = p;
                        m.p_corrected = p_corr;
                        let _ = m.write_as_fasta_bitset(
                            &mut output,
                            params.min_depth as u32,
                            &fasta_groups,
                        );
                    } else if params.output_bayes {
                        let bf = stats::bayes_factor_2x2(g1, g2, total_g1, total_g2);
                        let post =
                            stats::posterior_sex_linked(g1, g2, total_g1, total_g2, 0.01, 0.9);
                        write!(output, "{}\t{}", marker.id, marker.sequence).ok();
                        for &d in &marker.individual_depths {
                            write!(output, "\t{d}").ok();
                        }
                        let _ = writeln!(output, "\t{:.4}\t{:.4}", bf, post);
                    } else {
                        let _ = marker.write_as_table(&mut output);
                    }
                }
            }
        })?;
    }

    Ok(())
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
            !body.lines().any(|l| l.contains("CCCCCCCC") && !l.starts_with('#')),
            "balanced marker must not be emitted: {body}"
        );
    }

    #[test]
    fn signif_fdr_plus_fasta_is_rejected() {
        let dir = tempfile::tempdir().unwrap();
        let (table, pop) = write_fixture(dir.path());
        let out = dir.path().join("signif.fa");
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
        assert!(body.contains("#source:rsx-signif"), "header present: {body}");
        assert!(body.contains("AAAAAAAA"), "FDR should keep strong marker: {body}");
    }
}
