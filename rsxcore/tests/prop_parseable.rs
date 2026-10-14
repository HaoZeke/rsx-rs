// GPL-3.0-or-later
// Property tests over inputs the parsers accept: digit depths, DNA, popmaps,
// FASTA, and marker tables. Strategies stay inside those grammars.

use std::collections::BTreeMap;
use std::io::Write;

use proptest::prelude::*;
use proptest::string::string_regex;
use proptest::test_runner::Config;
use rsx_core::io::seq_reader::{count_sequences, pack_2bit, unpack_2bit};
use rsx_core::io::table_io::fast_parse_u16;
use rsx_core::kmer::canonical_kmer_hash;
use rsx_core::markers_table::{MarkersTableStream, ParserConfig};
use rsx_core::popmap::Popmap;
use rsx_core::stats::{benjamini_hochberg, bonferroni_correct, find_median, group_bias};
use rsx_core::test_method::{CorrectionMethod, TestMethod};

fn cases(n: u32) -> Config {
    Config::with_cases(n)
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|base| match base {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        })
        .collect()
}

fn reference_median(data: &[u16]) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    let mut sorted = data.to_vec();
    sorted.sort_unstable();
    let n = sorted.len();
    if n % 2 == 1 {
        f64::from(sorted[n / 2])
    } else {
        (f64::from(sorted[n / 2 - 1]) + f64::from(sorted[n / 2])) / 2.0
    }
}

proptest! {
    #![proptest_config(cases(64))]

    #[test]
    fn fast_parse_u16_matches_saturating_decimal(digits in "[0-9]{1,12}") {
        let parsed = digits.parse::<u64>().unwrap();
        let expected = u16::try_from(parsed).unwrap_or(u16::MAX);
        prop_assert_eq!(fast_parse_u16(digits.as_bytes()), expected);
    }

    #[test]
    fn pack_unpack_roundtrip_is_case_insensitive(seq in "[ACGT]{0,96}") {
        let packed = pack_2bit(seq.as_bytes());
        prop_assert_eq!(unpack_2bit(&packed), seq.as_bytes());
        prop_assert_eq!(packed, pack_2bit(seq.to_ascii_lowercase().as_bytes()));
    }

    #[test]
    fn canonical_kmer_hash_matches_reverse_complement(seq in "[ACGT]{1,40}") {
        let bytes = seq.as_bytes();
        for k in 1..=bytes.len().min(32) {
            let forward = canonical_kmer_hash(bytes, k);
            let reverse = canonical_kmer_hash(&revcomp(bytes), k);
            prop_assert_eq!(forward, reverse);
        }
    }

    #[test]
    fn canonical_kmer_hash_ignores_case(seq in "[ACGT]{1,24}") {
        let k = (seq.len() / 2).max(1);
        prop_assert_eq!(
            canonical_kmer_hash(seq.as_bytes(), k),
            canonical_kmer_hash(seq.to_ascii_lowercase().as_bytes(), k)
        );
    }

    #[test]
    fn counted_fasta_records_match_packed_keys(
        seqs in prop::collection::vec("[ACGT]{1,30}", 1..=4)
    ) {
        let mut body = String::new();
        let mut expected: BTreeMap<Vec<u8>, u16> = BTreeMap::new();
        for (index, seq) in seqs.iter().enumerate() {
            body.push_str(&format!(">r{index}\n{seq}\n"));
            let key = pack_2bit(seq.as_bytes());
            *expected.entry(key).or_insert(0) += 1;
        }
        let mut file = tempfile::NamedTempFile::new().unwrap();
        write!(file, "{body}").unwrap();
        let counts = count_sequences(file.path()).unwrap();
        prop_assert_eq!(counts.len(), expected.len());
        for (key, count) in expected {
            prop_assert_eq!(counts.get(&key).copied(), Some(count));
        }
    }

    #[test]
    fn popmap_round_trips_unique_assignments(
        assignments in prop::collection::btree_map(
            string_regex("[a-z]{1,8}").unwrap(),
            string_regex("[A-Z]{1,4}").unwrap(),
            1..=8,
        )
    ) {
        let mut body = String::new();
        for (individual, group) in &assignments {
            body.push_str(individual);
            body.push('\t');
            body.push_str(group);
            body.push_str("\r\n");
        }
        body.push_str("\r\n");
        let mut file = tempfile::NamedTempFile::new().unwrap();
        write!(file, "{body}").unwrap();
        let popmap = Popmap::from_file(file.path()).unwrap();
        prop_assert_eq!(popmap.n_individuals as usize, assignments.len());
        prop_assert_eq!(popmap.individual_groups.len(), assignments.len());
        let mut counts: BTreeMap<&str, u32> = BTreeMap::new();
        for (individual, group) in &assignments {
            prop_assert_eq!(popmap.get_group(individual), Some(group.as_str()));
            *counts.entry(group.as_str()).or_insert(0) += 1;
        }
        for (group, count) in counts {
            prop_assert_eq!(popmap.get_count(group), count);
        }
    }

    #[test]
    fn marker_table_parser_paths_agree(
        rows in (1usize..=4).prop_flat_map(|n_individuals| {
            prop::collection::vec(
                ("[ACGT]{4,12}", prop::collection::vec(0u16..=30, n_individuals)),
                1..=5,
            )
        })
    ) {
        let n_individuals = rows[0].1.len();
        let mut body = format!("#Number of markers : {}\nid\tsequence", rows.len());
        for index in 0..n_individuals {
            body.push_str(&format!("\ti{index}"));
        }
        body.push('\n');
        for (index, (sequence, depths)) in rows.iter().enumerate() {
            body.push_str(&format!("{index}\t{sequence}"));
            for depth in depths {
                body.push_str(&format!("\t{depth}"));
            }
            body.push('\n');
        }
        let mut file = tempfile::NamedTempFile::new().unwrap();
        write!(file, "{body}").unwrap();

        let full = MarkersTableStream::open(
            file.path(),
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
        let presence = MarkersTableStream::open(
            file.path(),
            None,
            ParserConfig {
                store_sequence: false,
                store_depths: false,
                compute_groups: false,
                min_depth: 1,
            },
        )
        .unwrap()
        .collect()
        .unwrap();
        let deep = MarkersTableStream::open(
            file.path(),
            None,
            ParserConfig {
                store_sequence: false,
                store_depths: true,
                compute_groups: false,
                min_depth: 5,
            },
        )
        .unwrap()
        .collect()
        .unwrap();

        prop_assert_eq!(full.len(), rows.len());
        for (index, ((sequence, depths), marker)) in rows.iter().zip(full.iter()).enumerate() {
            prop_assert_eq!(&marker.sequence, sequence);
            prop_assert_eq!(&marker.individual_depths, depths);
            let present = depths.iter().filter(|depth| **depth >= 1).count() as u32;
            prop_assert_eq!(marker.n_individuals, present);
            prop_assert_eq!(presence[index].n_individuals, present);
            let present_deep = depths.iter().filter(|depth| **depth >= 5).count() as u32;
            prop_assert_eq!(deep[index].n_individuals, present_deep);
            prop_assert_eq!(&deep[index].individual_depths, depths);
        }
    }

    #[test]
    fn median_matches_the_sorted_definition(mut data in prop::collection::vec(any::<u16>(), 0..=48)) {
        let expected = reference_median(&data);
        prop_assert_eq!(find_median(&mut data), expected);
    }

    #[test]
    fn group_bias_is_the_frequency_difference(
        n1 in 0u32..=30,
        total1 in 1u32..=30,
        n2 in 0u32..=30,
        total2 in 1u32..=30,
    ) {
        prop_assume!(n1 <= total1 && n2 <= total2);
        let bias = group_bias(n1, total1, n2, total2);
        let expected = f64::from(n1) / f64::from(total1) - f64::from(n2) / f64::from(total2);
        prop_assert!((bias - expected).abs() < 1e-12);
        prop_assert!((-1.0..=1.0).contains(&bias));
    }

    #[test]
    fn benjamini_hochberg_dominates_raw_p(p_values in prop::collection::vec(0.0f64..=1.0, 0..=12)) {
        let adjusted = benjamini_hochberg(&p_values);
        prop_assert_eq!(adjusted.len(), p_values.len());
        for (raw, q_value) in p_values.iter().zip(&adjusted) {
            prop_assert!(*q_value >= 0.0 && *q_value <= 1.0);
            prop_assert!(*q_value + 1e-9 >= *raw);
        }
        let mut order: Vec<usize> = (0..p_values.len()).collect();
        order.sort_by(|&left, &right| p_values[left].total_cmp(&p_values[right]));
        for window in order.windows(2) {
            prop_assert!(adjusted[window[0]] <= adjusted[window[1]] + 1e-9);
        }
    }

    #[test]
    fn bonferroni_matches_the_capped_product(p_value in 0.0f64..=1.0, n_markers in 0u64..=400) {
        let corrected = bonferroni_correct(p_value, n_markers);
        let expected = (p_value * n_markers as f64).min(1.0);
        prop_assert!((corrected - expected).abs() < 1e-9);
    }

    #[test]
    fn unknown_test_names_are_rejected(name in "[A-Za-z0-9_-]{1,12}") {
        let lower = name.to_ascii_lowercase();
        prop_assume!(!matches!(
            lower.as_str(),
            "chisq" | "chi-squared" | "chi2" | "fisher" | "fisher-exact" | "gtest" | "g-test" | "g"
        ));
        prop_assert!(TestMethod::parse_str(&name).is_err());
    }

    #[test]
    fn unknown_corrections_are_rejected(name in "[A-Za-z0-9_-]{1,16}") {
        let lower = name.to_ascii_lowercase();
        prop_assume!(!matches!(
            lower.as_str(),
            "bonferroni" | "bonf" | "fdr" | "bh" | "benjamini-hochberg" | "none" | "disabled"
        ));
        prop_assert!(CorrectionMethod::parse_str(&name).is_err());
    }
}
