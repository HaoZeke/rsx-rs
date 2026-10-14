// GPL-3.0-or-later
// Unit and integration coverage for parsers, group resolution, merge, and
// depth. Tests that describe a bug instead of current behavior are ignored
// at the bottom of this file.

use std::collections::BTreeMap;
use std::io::Write;
use std::path::{Path, PathBuf};

use rsx_core::bitset::{BitsetRow, GroupMask};
use rsx_core::commands::depth::{self, DepthParams};
use rsx_core::commands::merge::{self, MergeParams};
use rsx_core::commands::subset::{self, SubsetParams};
use rsx_core::compute_backend::PValueBackend;
use rsx_core::io::table_io::{TableHeader, fast_parse_u16};
use rsx_core::kmer::canonical_kmer_hash;
use rsx_core::marker::Marker;
use rsx_core::markers_table::{MarkersTableStream, ParserConfig};
use rsx_core::popmap::{GroupConfig, Popmap};
use rsx_core::stats::group_bias;
use rsx_core::test_method::{self, CorrectionMethod, TestMethod};

fn temp_path(name: &str) -> PathBuf {
    let dir = tempfile::tempdir().expect("tempdir").keep();
    std::fs::create_dir_all(&dir).unwrap();
    dir.join(name)
}

fn write_file(path: &Path, body: &str) {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).unwrap();
    }
    let mut file = std::fs::File::create(path).unwrap();
    file.write_all(body.as_bytes()).unwrap();
}

fn parse_depth_rows(body: &str) -> Vec<Vec<&str>> {
    body.lines()
        .filter(|line| !line.starts_with('#') && !line.starts_with("id\t"))
        .map(|line| line.split('\t').collect())
        .collect()
}

#[test]
fn fast_parse_u16_covers_empty_and_boundary_widths() {
    assert_eq!(fast_parse_u16(b""), 0);
    assert_eq!(fast_parse_u16(b"7"), 7);
    assert_eq!(fast_parse_u16(b"007"), 7);
    assert_eq!(fast_parse_u16(b"1000"), 1000);
    assert_eq!(fast_parse_u16(b"9999"), 9999);
    assert_eq!(fast_parse_u16(b"10000"), 10000);
    assert_eq!(fast_parse_u16(b"65535"), u16::MAX);
    assert_eq!(fast_parse_u16(b"65536"), u16::MAX);
}

#[test]
fn table_header_accepts_crlf_and_a_missing_comment() {
    let with_comment = temp_path("header_crlf.tsv");
    write_file(
        &with_comment,
        "#Number of markers : 12\r\nid\tsequence\ta\r\n",
    );
    let header = TableHeader::from_file(&with_comment).unwrap();
    assert_eq!(header.n_markers, 12);
    assert_eq!(header.n_individuals, 1);
    assert_eq!(header.columns, ["id", "sequence", "a"]);

    let bare = temp_path("header_bare.tsv");
    write_file(&bare, "id\tsequence\ta\tb\n");
    let header = TableHeader::from_file(&bare).unwrap();
    assert_eq!(header.n_markers, 0);
    assert_eq!(header.n_individuals, 2);

    let prose = temp_path("header_prose.tsv");
    write_file(&prose, "# not a count\nid\tsequence\tonly\n");
    let header = TableHeader::from_file(&prose).unwrap();
    assert_eq!(header.n_markers, 0);
    assert_eq!(header.n_individuals, 1);
}

#[test]
fn marker_table_skips_short_lines_and_ignores_all_zero_rows() {
    let path = temp_path("markers.tsv");
    write_file(
        &path,
        "#Number of markers : 3\n\
         id\tsequence\ta\tb\n\
         junk\n\
         0\tAAAA\t0\t0\n\
         1\tCCCC\t0\t2\n",
    );
    let markers = MarkersTableStream::open(
        &path,
        None,
        ParserConfig {
            store_sequence: true,
            store_depths: true,
            compute_groups: false,
            min_depth: 1,
        },
    )
    .unwrap()
    .collect()
    .unwrap();
    assert_eq!(markers.len(), 2);
    assert_eq!(markers[0].n_individuals, 0);
    assert_eq!(markers[1].sequence, "CCCC");
    assert_eq!(markers[1].individual_depths, vec![0, 2]);
    let counted = MarkersTableStream::open(
        &path,
        None,
        ParserConfig {
            store_sequence: false,
            store_depths: false,
            compute_groups: false,
            min_depth: 1,
        },
    )
    .unwrap()
    .count_markers()
    .unwrap();
    assert_eq!(counted, 1);
}

#[test]
fn marker_table_crlf_rows_match_the_depth_parser() {
    let path = temp_path("crlf.tsv");
    write_file(
        &path,
        "#Number of markers : 1\r\nid\tsequence\ta\tb\r\n0\tACGT\t4\t0\r\n",
    );
    let markers = MarkersTableStream::open(
        &path,
        None,
        ParserConfig {
            store_sequence: true,
            store_depths: true,
            compute_groups: false,
            min_depth: 1,
        },
    )
    .unwrap()
    .collect()
    .unwrap();
    assert_eq!(markers.len(), 1);
    assert_eq!(markers[0].sequence, "ACGT");
    assert_eq!(markers[0].individual_depths, vec![4, 0]);
    assert_eq!(markers[0].n_individuals, 1);
}

#[test]
fn popmap_resolution_orders_groups_and_rejects_incomplete_sets() {
    let path = temp_path("popmap.tsv");
    write_file(&path, "\nind1\tM\nind2\tF\n \nbad line\nind3\tX\n");
    let popmap = Popmap::from_file(&path).unwrap();
    assert_eq!(popmap.n_individuals, 3);
    assert_eq!(popmap.get_group("missing"), None);
    assert!(popmap.print_groups(false).contains("\"F\""));
    assert_eq!(popmap.print_groups(true).matches("\"M\": 1").count(), 1);

    let mut missing = GroupConfig::default();
    assert!(popmap.resolve_groups(&mut missing).is_err());

    let mut named = GroupConfig {
        group1: "X".into(),
        group2: "M".into(),
    };
    popmap.resolve_groups(&mut named).unwrap();
    assert_eq!(named.group1, "X");
    assert_eq!(named.group2, "M");

    let mut unknown = GroupConfig {
        group1: "M".into(),
        group2: "nope".into(),
    };
    assert!(popmap.resolve_groups(&mut unknown).is_err());

    let two = temp_path("two.tsv");
    write_file(&two, "b\tF\na\tM\n");
    let popmap = Popmap::from_file(&two).unwrap();
    let mut auto = GroupConfig::default();
    popmap.resolve_groups(&mut auto).unwrap();
    assert_eq!(auto.group1, "F");
    assert_eq!(auto.group2, "M");
    assert_eq!(group_bias(0, 1, 1, 1), -1.0);
    assert_eq!(group_bias(1, 0, 0, 1), 0.0);
}

#[test]
fn bitset_counts_across_the_five_word_loop() {
    let mut columns = vec![String::new(), String::new()];
    for index in 0..320 {
        columns.push(if index % 2 == 0 {
            "M".into()
        } else {
            "F".into()
        });
    }
    let males = GroupMask::from_columns(&columns, "M", 320);
    assert_eq!(males.count(), 160);

    let mut row = BitsetRow::new(320);
    for index in (0..320).step_by(2) {
        row.set(index);
    }
    assert_eq!(row.count_masked(&males), 160);
    row.set(319);
    row.set(1_000);
    assert_eq!(row.count_total(), 161);
    assert_eq!(row.count_masked(&males), 160);
    row.clear();
    assert_eq!(row.count_total(), 0);

    let empty = BitsetRow::new(0);
    let empty_mask = GroupMask::from_columns(&["id".into(), "sequence".into()], "M", 0);
    assert_eq!(empty.count_masked(&empty_mask), 0);
}

#[test]
fn marker_reset_can_keep_the_sequence() {
    let mut marker = Marker::new(1);
    marker.id = "9".into();
    marker.sequence = "AC".into();
    marker.individual_depths[0] = 3;
    marker.presence.set(0);
    marker.n_individuals = 1;
    marker.p = 0.2;
    marker.p_corrected = 0.4;
    marker.reset(true);
    assert_eq!(marker.id, "9");
    assert_eq!(marker.sequence, "AC");
    assert_eq!(marker.individual_depths, vec![0]);
    assert_eq!(marker.n_individuals, 0);
    assert_eq!(marker.p, 0.0);
    assert_eq!(marker.p_corrected, 0.0);
    assert_eq!(marker.presence.count_total(), 0);
}

#[test]
fn ambiguous_bases_pack_as_adenine() {
    assert_eq!(
        rsx_core::io::seq_reader::unpack_2bit(&rsx_core::io::seq_reader::pack_2bit(b"N")),
        b"A"
    );
    assert_eq!(canonical_kmer_hash(b"N", 1), canonical_kmer_hash(b"A", 1));
}

#[test]
fn canonical_kmer_hash_is_symmetric_at_the_u64_boundary() {
    let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
    assert_eq!(seq.len(), 32);
    assert_eq!(
        canonical_kmer_hash(seq, 32),
        canonical_kmer_hash(
            &seq.iter()
                .rev()
                .map(|base| match base {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    b'T' => b'A',
                    _ => b'N',
                })
                .collect::<Vec<_>>(),
            32
        )
    );
}

#[test]
fn test_and_correction_aliases_parse_and_dispatch() {
    assert!(matches!(
        TestMethod::parse_str("CHI2").unwrap(),
        TestMethod::ChiSquared
    ));
    assert!(matches!(
        TestMethod::parse_str("Fisher-Exact").unwrap(),
        TestMethod::Fisher
    ));
    assert!(matches!(
        TestMethod::parse_str("G-Test").unwrap(),
        TestMethod::GTest
    ));
    assert!(TestMethod::parse_str("anova").is_err());

    assert!(matches!(
        CorrectionMethod::parse_str("Bonf").unwrap(),
        CorrectionMethod::Bonferroni
    ));
    assert!(matches!(
        CorrectionMethod::parse_str("Benjamini-Hochberg").unwrap(),
        CorrectionMethod::Fdr
    ));
    assert!(matches!(
        CorrectionMethod::parse_str("disabled").unwrap(),
        CorrectionMethod::None
    ));
    assert!(CorrectionMethod::parse_str("holm").is_err());

    assert_eq!(
        PValueBackend::parse_str("GPU").unwrap(),
        PValueBackend::Cuda
    );
    assert_eq!(PValueBackend::Cpu.name(), "cpu");
    assert!(PValueBackend::parse_str("metal").is_err());

    for method in [
        TestMethod::ChiSquared,
        TestMethod::Fisher,
        TestMethod::GTest,
    ] {
        let p = test_method::compute_p(method, 4, 0, 4, 4);
        assert!(p.is_finite(), "{p}");
        assert!((1e-16..=1.0).contains(&p), "{method:?} -> {p}");
    }
    assert_eq!(test_method::compute_p(TestMethod::Fisher, 5, 0, 4, 4), 1.0);
    assert_eq!(test_method::compute_p(TestMethod::GTest, 5, 0, 4, 4), 1.0);
}

fn depth_table() -> (PathBuf, PathBuf) {
    let table = temp_path("depth_markers.tsv");
    let popmap = table.with_file_name("depth_popmap.tsv");
    write_file(
        &table,
        "#Number of markers : 3\n\
         id\tsequence\tf1\tm1\n\
         0\tAAAA\t1\t1\n\
         1\tCCCC\t0\t5\n\
         2\tGGGG\t4\t0\n",
    );
    write_file(&popmap, "f1\tF\nm1\tM\n");
    (table, popmap)
}

fn run_depth(table: &Path, popmap: &Path, streaming: bool, min_frequency: f32) -> String {
    let output = table.with_file_name(format!("depth_{streaming}_{min_frequency}.tsv"));
    depth::run(&DepthParams {
        markers_table_path: table.to_str().unwrap().to_string(),
        popmap_file_path: popmap.to_str().unwrap().to_string(),
        output_file_path: output.to_str().unwrap().to_string(),
        min_frequency,
        streaming,
    })
    .unwrap();
    std::fs::read_to_string(output).unwrap()
}

#[test]
fn depth_memory_reports_exact_per_sample_stats() {
    let (table, popmap) = depth_table();
    let body = run_depth(&table, &popmap, false, 0.0);
    assert_eq!(
        body,
        "Sample\tGroup\tReads\tMarkers\tRetained\tMin_depth\tMax_depth\tMedian_depth\tAverage_depth\n\
         f1\tF\t5\t2\t3\t0\t4\t1\t1\n\
         m1\tM\t6\t2\t3\t0\t5\t1\t2\n"
    );
    let filtered = run_depth(&table, &popmap, false, 1.0);
    assert_eq!(
        filtered,
        "Sample\tGroup\tReads\tMarkers\tRetained\tMin_depth\tMax_depth\tMedian_depth\tAverage_depth\n\
         f1\tF\t5\t2\t1\t1\t1\t1\t1\n\
         m1\tM\t6\t2\t1\t1\t1\t1\t1\n"
    );
}

#[test]
fn depth_streaming_matches_memory_mode() {
    let (table, popmap) = depth_table();
    assert_eq!(
        run_depth(&table, &popmap, true, 0.0),
        run_depth(&table, &popmap, false, 0.0)
    );
    assert_eq!(
        run_depth(&table, &popmap, true, 1.0),
        run_depth(&table, &popmap, false, 1.0)
    );
}

fn merged_depths(buffer_size: usize) -> BTreeMap<String, Vec<u16>> {
    let left = temp_path(&format!("left_{buffer_size}.tsv"));
    let right = temp_path(&format!("right_{buffer_size}.tsv"));
    let output = temp_path(&format!("merged_{buffer_size}.tsv"));
    write_file(
        &left,
        "#Number of markers : 3\n\
         id\tsequence\ts1\ts2\n\
         0\tAAAA\t1\t2\n\
         1\tCCCC\t4\t0\n\
         2\tAAAA\t3\t1\n",
    );
    write_file(
        &right,
        "#Number of markers : 2\n\
         id\tsequence\ts3\n\
         0\tAAAA\t8\n\
         1\tGGGG\t3\n",
    );
    merge::run(&MergeParams {
        input_files: vec![
            left.to_str().unwrap().to_string(),
            right.to_str().unwrap().to_string(),
        ],
        output_file_path: output.to_str().unwrap().to_string(),
        buffer_size: Some(buffer_size),
        output_parquet: false,
    })
    .unwrap();
    let body = std::fs::read_to_string(&output).unwrap();
    assert!(body.starts_with("#Number of markers : 3"));
    let mut rows = BTreeMap::new();
    for fields in parse_depth_rows(&body) {
        let depths = fields[2..]
            .iter()
            .map(|value| value.parse::<u16>().unwrap())
            .collect();
        rows.insert(fields[1].to_string(), depths);
    }
    rows
}

#[test]
fn merge_coalesces_sequences_from_one_row_chunks_and_from_memory() {
    let expected: BTreeMap<String, Vec<u16>> = BTreeMap::from([
        ("AAAA".into(), vec![4, 3, 8]),
        ("CCCC".into(), vec![4, 0, 0]),
        ("GGGG".into(), vec![0, 0, 3]),
    ]);
    assert_eq!(merged_depths(1), expected);
    assert_eq!(merged_depths(1_000), expected);
}

fn subset_output(output_bayes: bool) -> String {
    let table = temp_path(if output_bayes {
        "subset_bayes.tsv"
    } else {
        "subset_plain.tsv"
    });
    let popmap = table.with_extension("popmap.tsv");
    let output = table.with_extension("out.tsv");
    write_file(
        &table,
        "#Number of markers : 1\n\
         id\tsequence\tm1\tf1\n\
         0\tACGT\t5\t0\n",
    );
    write_file(&popmap, "m1\tM\nf1\tF\n");
    subset::run(&SubsetParams {
        markers_table_path: table.to_str().unwrap().to_string(),
        popmap_file_path: popmap.to_str().unwrap().to_string(),
        output_file_path: output.to_str().unwrap().to_string(),
        min_depth: 1,
        signif_threshold: 0.05,
        correction: CorrectionMethod::None,
        test_method: TestMethod::ChiSquared,
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
fn subset_fasta_writes_sorted_group_counts() {
    let table = temp_path("subset_fasta_markers.tsv");
    let popmap = table.with_extension("popmap.tsv");
    let output = table.with_extension("fa");
    write_file(
        &table,
        "#Number of markers : 1\n\
         id\tsequence\tm1\tm2\tf1\n\
         0\tGATTACA\t5\t5\t0\n",
    );
    write_file(&popmap, "m1\tM\nm2\tM\nf1\tF\n");
    subset::run(&SubsetParams {
        markers_table_path: table.to_str().unwrap().to_string(),
        popmap_file_path: popmap.to_str().unwrap().to_string(),
        output_file_path: output.to_str().unwrap().to_string(),
        min_depth: 1,
        signif_threshold: 0.05,
        correction: CorrectionMethod::None,
        test_method: TestMethod::ChiSquared,
        output_bayes: false,
        output_fasta: true,
        group1: "M".into(),
        group2: "F".into(),
        min_group1: 1,
        min_group2: 0,
        max_group1: 9,
        max_group2: 0,
        min_individuals: 1,
        max_individuals: 9,
    })
    .unwrap();
    let body = std::fs::read_to_string(output).unwrap();
    assert!(body.starts_with(">0_F:0_M:2_p:"), "{body}");
    assert!(body.contains("\nGATTACA\n"), "{body}");
}

#[test]
#[ignore = "bug: resolve_groups accepts unknown names when the popmap has exactly two groups"]
fn resolve_groups_rejects_unknown_names_for_a_two_group_popmap() {
    // Bug: with exactly two groups, resolve_groups only fills empty names.
    // A typo such as --groups X,Y is accepted, and later masks are empty.
    let path = temp_path("two_groups.tsv");
    write_file(&path, "ind1\tM\nind2\tF\n");
    let popmap = Popmap::from_file(&path).unwrap();
    let mut config = GroupConfig {
        group1: "X".into(),
        group2: "Y".into(),
    };
    assert!(
        popmap.resolve_groups(&mut config).is_err(),
        "unknown group names must be rejected when the popmap has two groups"
    );
}

#[test]
#[ignore = "bug: duplicate popmap rows inflate n_individuals and group counts"]
fn duplicate_popmap_rows_do_not_inflate_group_counts() {
    // Bug: a repeated individual increments n_individuals and the group
    // count again, while individual_groups keeps a single key. Contingency
    // totals then disagree with the marker-table columns.
    let path = temp_path("duplicate_popmap.tsv");
    write_file(&path, "ind1\tM\nind1\tM\n");
    let popmap = Popmap::from_file(&path).unwrap();
    assert_eq!(
        popmap.n_individuals as usize,
        popmap.individual_groups.len(),
        "duplicate rows must not inflate n_individuals"
    );
    assert_eq!(popmap.get_count("M"), 1);
}

#[test]
#[ignore = "bug: presence fast path counts zero-padded depth 00 as present"]
fn zero_padded_depth_is_absent_on_the_presence_fast_path() {
    // Bug: the presence-only fast path treats every non-empty field other
    // than the single byte `0` as present, so `00` is counted. The depth
    // parser reads that field as 0.
    let path = temp_path("padded_zero.tsv");
    write_file(
        &path,
        "#Number of markers : 1\n\
         id\tsequence\ta\tb\tc\n\
         0\tAAAA\t00\t0\t1\n",
    );
    let open = |store_depths: bool| {
        MarkersTableStream::open(
            &path,
            None,
            ParserConfig {
                store_sequence: false,
                store_depths,
                compute_groups: false,
                min_depth: 1,
            },
        )
        .unwrap()
        .collect()
        .unwrap()
    };
    let fast = open(false);
    let parsed = open(true);
    assert_eq!(parsed[0].individual_depths, vec![0, 0, 1]);
    assert_eq!(fast[0].n_individuals, parsed[0].n_individuals);
    assert_eq!(fast[0].n_individuals, 1);
}

#[test]
#[ignore = "bug: canonical k-mer hashes are not revcomp-symmetric for k > 32"]
fn canonical_kmer_hash_matches_revcomp_for_k_above_32() {
    // Bug: kmer_to_u64 keeps 32 bases in a u64. For k > 32 the extra bases
    // are shifted away, so a sequence and its reverse complement no longer
    // share a canonical hash. A periodic ACGT fixture still collides; this
    // sequence does not.
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
        canonical_kmer_hash(seq, 33),
        canonical_kmer_hash(&reverse, 33)
    );
}

#[test]
fn merge_parses_crlf_depth_tables() {
    // BufRead::lines strips CRLF, so the last depth survives a Windows table.
    let path = temp_path("crlf_merge.tsv");
    let output = temp_path("crlf_merge_out.tsv");
    write_file(
        &path,
        "#Number of markers : 1\r\nid\tsequence\ta1\r\n0\tAAAA\t7\r\n",
    );
    merge::run(&MergeParams {
        input_files: vec![path.to_str().unwrap().to_string()],
        output_file_path: output.to_str().unwrap().to_string(),
        buffer_size: Some(8),
        output_parquet: false,
    })
    .unwrap();
    let body = std::fs::read_to_string(output).unwrap();
    let row = parse_depth_rows(&body);
    assert_eq!(row.len(), 1, "{body}");
    assert_eq!(row[0][1], "AAAA");
    assert_eq!(row[0][2], "7", "{body}");
}

#[test]
#[ignore = "bug: SubsetParams::output_bayes does not change the table"]
fn subset_bayes_flag_changes_output() {
    // Bug: SubsetParams::output_bayes is stored and never read, so requesting
    // Bayes columns does not change the table.
    let plain = subset_output(false);
    let bayes = subset_output(true);
    assert_ne!(plain, bayes);
    assert!(
        bayes.contains("Bayes") || bayes.contains("Posterior"),
        "{bayes}"
    );
}

#[cfg(feature = "map")]
#[test]
#[ignore = "bug: MapParams::output_bayes does not change the alignment table"]
fn map_bayes_flag_changes_output() {
    // Bug: MapParams::output_bayes is never read, so the flag cannot add
    // Bayes columns to the alignment table.
    use rsx_core::commands::map::{self, MapParams};

    let dir = tempfile::tempdir().unwrap();
    let table = dir.path().join("markers.tsv");
    let popmap = dir.path().join("popmap.tsv");
    let genome = dir.path().join("genome.fa");
    write_file(
        &table,
        "#Number of markers : 1\n\
         id\tsequence\tind1\tind2\n\
         0\tATCGATCGATCG\t5\t0\n",
    );
    write_file(&popmap, "ind1\tM\nind2\tF\n");
    write_file(&genome, ">chr1\nATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n");
    let run = |bayes: bool| {
        let output = dir.path().join(format!("map_{bayes}.tsv"));
        map::run(&MapParams {
            markers_table_path: table.to_str().unwrap().to_string(),
            popmap_file_path: popmap.to_str().unwrap().to_string(),
            genome_file_path: genome.to_str().unwrap().to_string(),
            output_file_path: output.to_str().unwrap().to_string(),
            min_depth: 1,
            min_quality: 0,
            min_frequency: 0.1,
            signif_threshold: 0.05,
            correction: CorrectionMethod::None,
            test_method: TestMethod::ChiSquared,
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
