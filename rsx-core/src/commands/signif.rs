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

    if matches!(params.correction, CorrectionMethod::Fdr) {
        struct FdrEntry {
            id: String,
            seq: Vec<u8>,
            depths: Vec<u16>,
            g1: u32,
            g2: u32,
        }
        let mut p_values: Vec<f64> = Vec::new();
        let mut marker_data: Vec<FdrEntry> = Vec::new();

        source.for_each(|marker| {
            if marker.n_individuals > 0 {
                let g1 = marker.presence.count_masked(&mask_g1);
                let g2 = marker.presence.count_masked(&mask_g2);
                let p = compute_p(params.test_method, g1, g2, total_g1, total_g2);
                p_values.push(p);
                marker_data.push(FdrEntry {
                    id: marker.id.clone(),
                    seq: marker.sequence.as_bytes().to_vec(),
                    depths: marker.individual_depths.clone(),
                    g1,
                    g2,
                });
            }
        })?;

        let q_values = stats::benjamini_hochberg(&p_values);

        for (i, (_, &q)) in p_values.iter().zip(q_values.iter()).enumerate() {
            if q < threshold {
                let entry = &marker_data[i];
                let seq_str = std::str::from_utf8(&entry.seq).unwrap_or("?");
                write!(output, "{}\t{}", entry.id, seq_str)?;
                for &d in &entry.depths {
                    write!(output, "\t{d}")?;
                }
                if params.output_bayes {
                    let bf = stats::bayes_factor_2x2(entry.g1, entry.g2, total_g1, total_g2);
                    let post = stats::posterior_sex_linked(
                        entry.g1, entry.g2, total_g1, total_g2, 0.01, 0.9,
                    );
                    writeln!(output, "\t{:.4}\t{:.4}", bf, post)?;
                } else {
                    writeln!(output)?;
                }
            }
        }
    } else {
        source.for_each(|marker| {
            if marker.n_individuals > 0 {
                let g1 = marker.presence.count_masked(&mask_g1);
                let g2 = marker.presence.count_masked(&mask_g2);
                let p = compute_p(params.test_method, g1, g2, total_g1, total_g2);

                if (p as f32) < corrected_threshold as f32 {
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
                        writeln!(output, "\t{:.4}\t{:.4}", bf, post).ok();
                    } else {
                        let _ = marker.write_as_table(&mut output);
                    }
                }
            }
        })?;
    }

    Ok(())
}
