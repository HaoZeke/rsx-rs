// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Correctness regressions: parsing, presence, association arithmetic, and identity.

use std::io::Write;
use std::path::Path;

fn temp_dir() -> std::path::PathBuf {
    tempfile::tempdir().expect("tempdir").keep()
}

#[test]
fn fast_parse_u16_rejects_non_digits_and_does_not_panic() {
    assert_eq!(rsx_core::io::table_io::fast_parse_u16(b"null"), 0);
    assert_eq!(rsx_core::io::table_io::fast_parse_u16(b"10a"), 0);
    assert_eq!(rsx_core::io::table_io::fast_parse_u16(b"-1"), 0);
    assert_eq!(rsx_core::io::table_io::fast_parse_u16(b"00"), 0);
    assert_eq!(rsx_core::io::table_io::fast_parse_u16(b"01"), 1);
    assert_eq!(rsx_core::io::table_io::fast_parse_u16(b"65535"), 65535);
    assert_eq!(rsx_core::io::table_io::fast_parse_u16(b"65536"), u16::MAX);
}

#[test]
fn table_header_rejects_more_individuals_than_u16() {
    let dir = temp_dir();
    let path = dir.join("wide.tsv");
    let mut f = std::fs::File::create(&path).unwrap();
    write!(f, "id\tsequence").unwrap();
    for i in 0..=u32::from(u16::MAX) {
        write!(f, "\ts{i}").unwrap();
    }
    writeln!(f).unwrap();
    drop(f);

    let err = rsx_core::io::table_io::TableHeader::from_file(&path).unwrap_err();
    let message = err.to_string();
    assert!(
        message.contains("65535"),
        "oversized cohort must be rejected, got {message}"
    );
}

fn presence(table: &Path, min_depth: u16) -> Vec<u32> {
    let config = rsx_core::markers_table::ParserConfig {
        store_sequence: false,
        store_depths: false,
        compute_groups: false,
        min_depth,
    };
    rsx_core::markers_table::MarkersTableStream::open(table, None, config)
        .unwrap()
        .collect()
        .unwrap()
        .into_iter()
        .map(|marker| marker.n_individuals)
        .collect()
}

#[test]
fn fast_presence_path_matches_parsed_depth_for_zero_spellings() {
    let dir = temp_dir();
    let path = dir.join("markers.tsv");
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, "#Number of markers : 4").unwrap();
    writeln!(f, "id\tsequence\tind1").unwrap();
    writeln!(f, "0\tAAAA\t00").unwrap();
    writeln!(f, "1\tCCCC\t0").unwrap();
    writeln!(f, "2\tGGGG\t000").unwrap();
    writeln!(f, "3\tTTTT\t10").unwrap();
    drop(f);

    assert_eq!(presence(&path, 1), vec![0, 0, 0, 1]);
}

#[test]
fn min_depth_zero_counts_zero_depth_individuals_as_present() {
    let dir = temp_dir();
    let path = dir.join("markers.tsv");
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, "#Number of markers : 2").unwrap();
    writeln!(f, "id\tsequence\tind1").unwrap();
    writeln!(f, "0\tAAAA\t0").unwrap();
    writeln!(f, "1\tCCCC\t5").unwrap();
    drop(f);

    assert_eq!(presence(&path, 0), vec![1, 1]);
}

#[test]
fn full_parser_does_not_panic_on_non_digit_depth_fields() {
    let dir = temp_dir();
    let path = dir.join("markers.tsv");
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, "#Number of markers : 2").unwrap();
    writeln!(f, "id\tsequence\tind1\tind2").unwrap();
    writeln!(f, "0\tAAAA\tnull\t5").unwrap();
    writeln!(f, "1\tCCCC\t10a\t0").unwrap();
    drop(f);

    let config = rsx_core::markers_table::ParserConfig {
        store_sequence: true,
        store_depths: true,
        compute_groups: false,
        min_depth: 1,
    };
    let markers = rsx_core::markers_table::MarkersTableStream::open(&path, None, config)
        .unwrap()
        .collect()
        .unwrap();
    assert_eq!(markers[0].individual_depths, vec![0, 5]);
    assert_eq!(markers[0].n_individuals, 1);
    assert_eq!(markers[1].individual_depths, vec![0, 0]);
    assert_eq!(markers[1].n_individuals, 0);
    assert_eq!(presence(&path, 1), vec![1, 0]);
}

#[test]
fn association_tests_do_not_overflow_u32_group_totals() {
    let p = rsx_core::stats::g_test(0, 0, 3_000_000_000, 2_000_000_000);
    assert!(
        p > 0.5,
        "all-absent table is not an association, g-test p={p}"
    );

    let fisher = rsx_core::stats::fisher_exact(0, 0, 3_000_000_000, 2_000_000_000);
    assert!(
        fisher > 0.5,
        "all-absent table is not an association, fisher p={fisher}"
    );

    // Equal prevalence, but the group totals sum past u32::MAX. The factor is
    // small and must stay finite; a wrapped margin underflows the beta update.
    let bf = rsx_core::stats::bayes_factor_2x2(
        1_500_000_000,
        1_000_000_000,
        3_000_000_000,
        2_000_000_000,
    );
    assert!(
        bf.is_finite() && bf > 0.0,
        "Bayes factor must be finite for valid large counts, got {bf}"
    );

    let post =
        rsx_core::stats::posterior_sex_linked(100, 0, 3_000_000_000, 2_000_000_000, 0.01, 0.9);
    assert!(
        post.is_finite() && (0.0..=1.0).contains(&post),
        "posterior must stay in range for valid large counts, got {post}"
    );
}

#[test]
fn chi_squared_invalid_table_is_not_significant() {
    let p = rsx_core::stats::p_association(2, 0, 1, 1);
    assert_eq!(
        p, 1.0,
        "present count above the group total must not look significant"
    );
}

#[test]
fn empirical_bayes_empty_cohort_is_finite() {
    let (pi, posteriors) = rsx_core::stats::empirical_bayes_em(&[], 10, 10, 0.9, 5);
    assert!(pi.is_finite() && (0.0..=1.0).contains(&pi), "pi={pi}");
    assert!(posteriors.is_empty());
}

#[test]
fn canonical_kmer_hash_accepts_zero_and_does_not_truncate_past_32() {
    let hashed = rsx_core::kmer::canonical_kmer_hash(b"ACGT", 0);
    assert_eq!(hashed, rsx_core::kmer::canonical_kmer_hash(b"ACGT", 0));
    let groups = rsx_core::kmer::group_by_kmer(&[b"ACGT".to_vec(), b"TTTT".to_vec()], 0);
    assert_eq!(groups.len(), 2);

    let seq = b"ACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
    let other = b"TCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
    assert_eq!(seq.len(), 33);
    assert_ne!(
        rsx_core::kmer::canonical_kmer_hash(seq, 33),
        rsx_core::kmer::canonical_kmer_hash(other, 33),
        "a 33-mer hash must keep the first base"
    );

    let rev = b"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGT";
    assert_eq!(
        rsx_core::kmer::canonical_kmer_hash(seq, 33),
        rsx_core::kmer::canonical_kmer_hash(rev, 33),
        "a sequence and its reverse complement must share a canonical hash"
    );
}

#[test]
fn popmap_rejects_duplicate_individuals_and_unknown_groups() {
    let dir = temp_dir();
    let dup = dir.join("dup.tsv");
    std::fs::write(&dup, "ind1\tM\nind1\tF\nind2\tF\n").unwrap();
    let err = rsx_core::popmap::Popmap::from_file(&dup).unwrap_err();
    assert!(
        err.to_string().contains("ind1"),
        "duplicate id must be reported, got {err}"
    );

    let ok = dir.join("ok.tsv");
    std::fs::write(&ok, "ind1\tM\nind2\tF\n").unwrap();
    let popmap = rsx_core::popmap::Popmap::from_file(&ok).unwrap();
    let mut config = rsx_core::popmap::GroupConfig {
        group1: "male".to_string(),
        group2: "female".to_string(),
    };
    let err = popmap.resolve_groups(&mut config).unwrap_err();
    assert!(
        err.contains("male") || err.contains("not found"),
        "unknown groups must be rejected when the popmap has two groups, got {err}"
    );

    let mut same = rsx_core::popmap::GroupConfig {
        group1: "M".to_string(),
        group2: "M".to_string(),
    };
    assert!(popmap.resolve_groups(&mut same).is_err());
}

#[test]
fn process_rejects_files_that_collapse_to_one_individual_name() {
    let dir = temp_dir();
    let input = dir.join("reads");
    std::fs::create_dir_all(&input).unwrap();
    std::fs::write(input.join("ind1.fa"), ">r\nAAAAAAAA\n").unwrap();
    std::fs::write(input.join("ind1.fasta"), ">r\nCCCCCCCC\n").unwrap();
    let err = rsx_core::commands::process::run(&rsx_core::commands::process::ProcessParams {
        input_dir_path: input.to_str().unwrap().to_string(),
        output_file_path: dir.join("markers.tsv").to_str().unwrap().to_string(),
        n_threads: 1,
        min_depth: 1,
        kmer_dedup: None,
    });
    assert!(
        err.is_err(),
        "ind1.fa and ind1.fasta must not share a depth column"
    );
}

#[test]
fn distrib_does_not_panic_when_duplicate_columns_exceed_the_popmap() {
    let dir = temp_dir();
    let table = dir.join("markers.tsv");
    let mut f = std::fs::File::create(&table).unwrap();
    writeln!(f, "#Number of markers : 1").unwrap();
    writeln!(f, "id\tsequence\tind1\tind1").unwrap();
    writeln!(f, "0\tAAAA\t5\t5").unwrap();
    drop(f);
    std::fs::write(dir.join("popmap.tsv"), "ind1\tM\nind2\tF\n").unwrap();

    let err = rsx_core::commands::distrib::run(&rsx_core::commands::distrib::DistribParams {
        markers_table_path: table.to_str().unwrap().to_string(),
        popmap_file_path: dir.join("popmap.tsv").to_str().unwrap().to_string(),
        output_file_path: dir.join("distrib.tsv").to_str().unwrap().to_string(),
        min_depth: 1,
        signif_threshold: 0.05,
        correction: rsx_core::test_method::CorrectionMethod::Bonferroni,
        test_method: rsx_core::test_method::TestMethod::ChiSquared,
        output_bayes: false,
        bayes_model: rsx_core::stats::DirectionalModel::directional_screening_v1(),
        group1: "M".to_string(),
        group2: "F".to_string(),
    });
    assert!(
        err.is_err(),
        "duplicate sample columns must be reported instead of panicking"
    );
}

#[test]
fn merge_reads_a_table_without_a_marker_count_comment() {
    let dir = temp_dir();
    let input = dir.join("markers.tsv");
    std::fs::write(&input, "id\tsequence\ta\tb\n0\tAAAA\t1\t2\n1\tCCCC\t3\t4\n").unwrap();
    let output = dir.join("merged.tsv");
    rsx_core::commands::merge::run(&rsx_core::commands::merge::MergeParams {
        input_files: vec![input.to_str().unwrap().to_string()],
        output_file_path: output.to_str().unwrap().to_string(),
        buffer_size: Some(8),
        output_parquet: false,
    })
    .unwrap();

    let body = std::fs::read_to_string(&output).unwrap();
    let mut lines = body.lines();
    let _comment = lines.next().unwrap();
    let header = lines.next().unwrap();
    assert_eq!(header, "id\tsequence\ta\tb");
    let rows: Vec<&str> = lines.filter(|line| !line.is_empty()).collect();
    assert_eq!(rows.len(), 2, "both data rows must survive, got {body}");
    assert!(rows.iter().any(|row| row.contains("AAAA\t1\t2")), "{body}");
    assert!(rows.iter().any(|row| row.contains("CCCC\t3\t4")), "{body}");
}

#[test]
fn merge_keeps_sample_offsets_past_u16() {
    let dir = temp_dir();
    let wide = dir.join("wide.tsv");
    let extra = dir.join("extra.tsv");
    let n_wide = u32::from(u16::MAX) + 1;

    let mut f = std::fs::File::create(&wide).unwrap();
    writeln!(f, "#Number of markers : 1").unwrap();
    write!(f, "id\tsequence").unwrap();
    for i in 0..n_wide {
        write!(f, "\tw{i}").unwrap();
    }
    writeln!(f).unwrap();
    write!(f, "0\tAAAA").unwrap();
    for _ in 0..n_wide {
        write!(f, "\t0").unwrap();
    }
    writeln!(f).unwrap();
    drop(f);

    std::fs::write(
        &extra,
        "#Number of markers : 1\nid\tsequence\textra\n0\tAAAA\t7\n",
    )
    .unwrap();
    let output = dir.join("merged.tsv");
    rsx_core::commands::merge::run(&rsx_core::commands::merge::MergeParams {
        input_files: vec![
            wide.to_str().unwrap().to_string(),
            extra.to_str().unwrap().to_string(),
        ],
        output_file_path: output.to_str().unwrap().to_string(),
        buffer_size: Some(8),
        output_parquet: false,
    })
    .unwrap();

    let body = std::fs::read_to_string(&output).unwrap();
    let data = body
        .lines()
        .find(|line| line.contains("AAAA"))
        .expect("merged sequence");
    let fields: Vec<&str> = data.split('\t').collect();
    assert_eq!(fields.len(), 2 + n_wide as usize + 1);
    assert_eq!(fields[2], "0", "first wide sample was overwritten: {data}");
    assert_eq!(fields.last().copied(), Some("7"));
}

#[test]
fn canonical_kmer_hash_matches_revcomp_for_the_reported_33mer() {
    let seq = b"TTAGTTGTGCCGCAGCGAAGTAGTGCTTGAAATATG";
    assert_eq!(seq.len(), 36);
    let reverse = seq
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect::<Vec<_>>();
    assert_eq!(
        rsx_core::kmer::canonical_kmer_hash(seq, 33),
        rsx_core::kmer::canonical_kmer_hash(&reverse, 33)
    );
}

fn subset_output(output_bayes: bool) -> String {
    let dir = temp_dir();
    let table = dir.join(if output_bayes {
        "subset_bayes.tsv"
    } else {
        "subset_plain.tsv"
    });
    let popmap = table.with_extension("popmap.tsv");
    let output = table.with_extension("out.tsv");
    let mut markers = std::fs::File::create(&table).unwrap();
    write!(
        markers,
        "#Number of markers : 1\n\
         id\tsequence\tm1\tf1\n\
         0\tACGT\t5\t0\n"
    )
    .unwrap();
    drop(markers);
    std::fs::write(&popmap, "m1\tM\nf1\tF\n").unwrap();
    rsx_core::commands::subset::run(&rsx_core::commands::subset::SubsetParams {
        markers_table_path: table.to_str().unwrap().to_string(),
        popmap_file_path: popmap.to_str().unwrap().to_string(),
        output_file_path: output.to_str().unwrap().to_string(),
        min_depth: 1,
        signif_threshold: 0.05,
        correction: rsx_core::test_method::CorrectionMethod::None,
        test_method: rsx_core::test_method::TestMethod::ChiSquared,
        output_bayes,
        output_fasta: false,
        group1: "M".into(),
        group2: "F".into(),
        min_group1: 0,
        min_group2: 0,
        max_group1: 9,
        max_group2: 9,
        min_individuals: 0,
        max_individuals: 9,
    })
    .unwrap();
    std::fs::read_to_string(output).unwrap()
}

#[test]
fn subset_bayes_flag_changes_output() {
    let plain = subset_output(false);
    let bayes = subset_output(true);
    assert_ne!(plain, bayes);
    assert!(
        bayes.contains("Bayes_Factor") && bayes.contains("Posterior_SexLinked"),
        "{bayes}"
    );
    let row = bayes
        .lines()
        .find(|line| line.contains("ACGT"))
        .expect("subset keeps the marker");
    assert_eq!(row.split('\t').count(), 6, "{row}");
}

#[cfg(feature = "map")]
#[test]
fn map_bayes_flag_changes_output() {
    let dir = temp_dir();
    let table = dir.join("markers.tsv");
    let popmap = dir.join("popmap.tsv");
    let genome = dir.join("genome.fa");
    let mut markers = std::fs::File::create(&table).unwrap();
    write!(
        markers,
        "#Number of markers : 1\n\
         id\tsequence\tind1\tind2\n\
         0\tATCGATCGATCG\t5\t0\n"
    )
    .unwrap();
    drop(markers);
    std::fs::write(&popmap, "ind1\tM\nind2\tF\n").unwrap();
    std::fs::write(&genome, ">chr1\nATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n").unwrap();
    let run = |bayes: bool| {
        let output = dir.join(format!("map_{bayes}.tsv"));
        rsx_core::commands::map::run(&rsx_core::commands::map::MapParams {
            markers_table_path: table.to_str().unwrap().to_string(),
            popmap_file_path: popmap.to_str().unwrap().to_string(),
            genome_file_path: genome.to_str().unwrap().to_string(),
            output_file_path: output.to_str().unwrap().to_string(),
            min_depth: 1,
            min_quality: 0,
            min_frequency: 0.1,
            signif_threshold: 0.05,
            correction: rsx_core::test_method::CorrectionMethod::None,
            test_method: rsx_core::test_method::TestMethod::ChiSquared,
            output_bayes: bayes,
            group1: "M".into(),
            group2: "F".into(),
        })
        .unwrap();
        std::fs::read_to_string(output).unwrap()
    };
    let plain = run(false);
    let bayes = run(true);
    assert_ne!(plain, bayes);
    assert!(
        bayes.contains("Bayes") || bayes.contains("Posterior"),
        "{bayes}"
    );
}
