// GPL-3.0-or-later
// Strict numerical precision tests: verify rsx-rs matches C++ radsex exactly.

use std::io::Write;
use std::path::PathBuf;

fn test_dir() -> PathBuf {
    let dir = std::env::temp_dir().join("rsx_precision_test");
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn create_precision_markers(dir: &std::path::Path, n_ind: u16) -> PathBuf {
    let path = dir.join("precision_markers.tsv");
    let mut f = std::fs::File::create(&path).unwrap();
    let n_markers = 100;
    writeln!(f, "#Number of markers : {n_markers}").unwrap();
    write!(f, "id\tsequence").unwrap();
    for i in 0..n_ind {
        write!(f, "\tind{i}").unwrap();
    }
    writeln!(f).unwrap();

    let half = n_ind / 2;
    for m in 0..n_markers {
        let seq: String = (0..16)
            .map(|i| ["A", "T", "C", "G"][(m * 7 + i) % 4])
            .collect();
        write!(f, "{m}\t{seq}").unwrap();
        for j in 0..n_ind {
            let d = match m {
                // All present (should give p=1)
                0..=9 => 10,
                // Male-only (j < half)
                10..=19 => {
                    if j < half {
                        15
                    } else {
                        0
                    }
                }
                // Female-only (j >= half)
                20..=29 => {
                    if j >= half {
                        20
                    } else {
                        0
                    }
                }
                // Gradient: increasing male presence
                30..=49 => {
                    if j < half && (j as u32) < (m as u32).saturating_sub(30) {
                        8
                    } else {
                        0
                    }
                }
                // Gradient: decreasing female presence
                50..=69 => {
                    if j >= half
                        && (j as u32).saturating_sub(half as u32) < (m as u32).saturating_sub(50)
                    {
                        12
                    } else {
                        0
                    }
                }
                // Mixed presence at various thresholds
                _ => {
                    if (m as u32 + j as u32) % 3 == 0 {
                        5
                    } else {
                        0
                    }
                }
            };
            write!(f, "\t{d}").unwrap();
        }
        writeln!(f).unwrap();
    }
    path
}

fn create_precision_popmap(dir: &std::path::Path, n_ind: u16) -> PathBuf {
    let path = dir.join("precision_popmap.tsv");
    let mut f = std::fs::File::create(&path).unwrap();
    let half = n_ind / 2;
    for i in 0..n_ind {
        let group = if i < half { "M" } else { "F" };
        writeln!(f, "ind{i}\t{group}").unwrap();
    }
    path
}

// === Chi-squared precision tests ===

#[test]
fn test_chi_squared_exact_values() {
    use rsx_core::stats;

    // Known values computed independently (Wolfram Alpha / scipy)
    // chi2 with Yates for 2x2 table: N*(|ad-bc| - N/2)^2 / (a+b)(c+d)(a+c)(b+d)

    // Case 1: 10 males present, 0 females, 10M 10F total
    let chi = stats::chi_squared_yates(10, 0, 10, 10);
    // N=20, |10*10 - 0*10| = 100, Yates = (100-10)^2 = 8100
    // chi = 20 * 8100 / (10*10*10*10) = 16.2
    assert!((chi - 16.2).abs() < 1e-10, "chi={chi}, expected 16.2");

    // Case 2: 5 males, 5 females, 10M 10F -- no association
    let chi = stats::chi_squared_yates(5, 5, 10, 10);
    // N=20, |5*10 - 5*10| = 0, Yates = max(0, 0-10) = 0
    // chi = 0
    assert!(chi == 0.0, "chi={chi}, expected 0.0");

    // Case 3: 3 males, 7 females, 10M 10F
    let chi = stats::chi_squared_yates(3, 7, 10, 10);
    // N=20, |3*10 - 7*10| = 40, Yates = (40-10)^2 = 900
    // chi = 20 * 900 / (10*10*10*10) = 1.8
    assert!((chi - 1.8).abs() < 1e-10, "chi={chi}, expected 1.8");
}

#[test]
fn test_p_value_monotonicity() {
    use rsx_core::stats;

    // P-values must decrease as the imbalance |g1 - expected| increases.
    // With total_g1=total_g2=10, the expected count under H0 is n/2.
    // As g1 increases from 5 to 10 (holding g1+g2=10 fixed), p should decrease.
    let mut prev_p = 1.0;
    for g1 in 5..=10 {
        let g2 = 10 - g1;
        let p = stats::p_association(g1, g2, 10, 10);
        assert!(
            p <= prev_p + 1e-12,
            "p={p} > prev_p={prev_p} at g1={g1}, g2={g2}"
        );
        prev_p = p;
    }
}

#[test]
fn test_p_value_floor() {
    use rsx_core::stats;

    // Extreme case: all males present, no females
    let p = stats::p_association(100, 0, 100, 100);
    assert!(p >= 1e-16, "p={p} below floor 1e-16");
    assert!(p < 0.05, "p={p} should be significant");
}

#[test]
fn test_p_value_symmetry() {
    use rsx_core::stats;

    // p(g1=10, g2=0) should equal p(g1=0, g2=10) with same totals
    let p1 = stats::p_association(10, 0, 20, 20);
    let p2 = stats::p_association(0, 10, 20, 20);
    assert!((p1 - p2).abs() < 1e-15, "asymmetry: p1={p1}, p2={p2}");
}

#[test]
fn test_bonferroni_identity() {
    use rsx_core::stats;

    // With 1 marker, corrected p = original p
    let p = 0.03;
    assert_eq!(stats::bonferroni_correct(p, 1), p);

    // With many markers, corrected p is capped at 1.0
    assert_eq!(stats::bonferroni_correct(0.5, 10), 1.0);
}

// === Bitset correctness tests ===

#[test]
fn test_bitset_group_counts_consistent() {
    use rsx_core::bitset::GroupMask;
    use rsx_core::markers_table::{MarkersTableStream, ParserConfig};
    use rsx_core::popmap::Popmap;

    let dir = test_dir().join("bitset_consistency");
    std::fs::create_dir_all(&dir).unwrap();

    for n_ind in [10u16, 20, 40, 100, 200] {
        let table = create_precision_markers(&dir, n_ind);
        let popmap_path = create_precision_popmap(&dir, n_ind);
        let popmap = Popmap::from_file(&popmap_path).unwrap();

        let config = ParserConfig {
            store_sequence: false,
            compute_groups: true,
            min_depth: 5,
        };

        let stream = MarkersTableStream::open(&table, Some(&popmap), config).unwrap();
        let mask_m = GroupMask::from_columns(&stream.groups, "M", stream.header.n_individuals);
        let mask_f = GroupMask::from_columns(&stream.groups, "F", stream.header.n_individuals);

        let expected_m = mask_m.count();
        let expected_f = mask_f.count();
        let half = n_ind / 2;

        // Verify masks have the right number of bits
        assert_eq!(
            expected_m, half as u32,
            "n_ind={n_ind}: mask_m count {expected_m} != expected {half}"
        );
        assert_eq!(
            expected_f,
            (n_ind - half) as u32,
            "n_ind={n_ind}: mask_f count {expected_f} != expected {}",
            n_ind - half
        );

        let mut marker_idx = 0u32;
        stream
            .for_each(|marker| {
                let g1 = marker.presence.count_masked(&mask_m);
                let g2 = marker.presence.count_masked(&mask_f);
                let total = marker.presence.count_total();

                // Group counts must sum to total (no individual in both groups)
                assert_eq!(
                    g1 + g2,
                    total,
                    "n_ind={n_ind} marker={marker_idx}: g1={g1} + g2={g2} != total={total}"
                );

                // n_individuals must equal bitset total
                assert_eq!(
                    total, marker.n_individuals,
                    "n_ind={n_ind} marker={marker_idx}: bitset total={total} != n_individuals={}",
                    marker.n_individuals
                );

                // Group counts must not exceed group size
                assert!(
                    g1 <= expected_m,
                    "n_ind={n_ind} marker={marker_idx}: g1={g1} > group_size={expected_m}"
                );
                assert!(
                    g2 <= expected_f,
                    "n_ind={n_ind} marker={marker_idx}: g2={g2} > group_size={expected_f}"
                );

                marker_idx += 1;
            })
            .unwrap();

        assert!(marker_idx > 0, "no markers processed for n_ind={n_ind}");
    }
}

// === End-to-end golden tests with known outputs ===

#[test]
fn test_distrib_exact_output() {
    // Known input -> known output, verified by hand calculation
    let dir = test_dir().join("distrib_exact");
    std::fs::create_dir_all(&dir).unwrap();

    // 4 individuals: 2M, 2F. 3 markers.
    let table = dir.join("markers.tsv");
    let mut f = std::fs::File::create(&table).unwrap();
    writeln!(f, "#Number of markers : 3").unwrap();
    writeln!(f, "id\tsequence\tind1\tind2\tind3\tind4").unwrap();
    writeln!(f, "0\tAAAA\t10\t10\t10\t10").unwrap(); // all present
    writeln!(f, "1\tBBBB\t10\t10\t0\t0").unwrap(); // M only
    writeln!(f, "2\tCCCC\t0\t0\t10\t10").unwrap(); // F only

    let popmap = dir.join("popmap.tsv");
    let mut f = std::fs::File::create(&popmap).unwrap();
    writeln!(f, "ind1\tM\nind2\tM\nind3\tF\nind4\tF").unwrap();

    let output = dir.join("distrib.tsv");
    rsx_core::commands::distrib::run(&rsx_core::commands::distrib::DistribParams {
        markers_table_path: table.to_str().unwrap().to_string(),
        popmap_file_path: popmap.to_str().unwrap().to_string(),
        output_file_path: output.to_str().unwrap().to_string(),
        min_depth: 1,
        signif_threshold: 0.05,
        correction: rsx_core::test_method::CorrectionMethod::Bonferroni,
        test_method: rsx_core::test_method::TestMethod::ChiSquared,
        output_bayes: false,
        group1: "M".to_string(),
        group2: "F".to_string(),
    })
    .unwrap();

    let content = std::fs::read_to_string(&output).unwrap();
    let lines: Vec<&str> = content.lines().collect();

    // Line for (M=2, F=2) should have 1 marker, p=1 (chi2=0)
    let line_2_2 = lines
        .iter()
        .find(|l| l.starts_with("2\t2\t"))
        .expect("missing line for M=2 F=2");
    let fields: Vec<&str> = line_2_2.split('\t').collect();
    assert_eq!(fields[2], "1", "marker count for (2,2)");
    assert_eq!(fields[3], "1", "p-value for equal distribution");

    // Line for (M=2, F=0) should have 1 marker, significant
    let line_2_0 = lines
        .iter()
        .find(|l| l.starts_with("2\t0\t"))
        .expect("missing line for M=2 F=0");
    let fields: Vec<&str> = line_2_0.split('\t').collect();
    assert_eq!(fields[2], "1", "marker count for (2,0)");

    // Line for (M=0, F=2) should have 1 marker
    let line_0_2 = lines
        .iter()
        .find(|l| l.starts_with("0\t2\t"))
        .expect("missing line for M=0 F=2");
    let fields: Vec<&str> = line_0_2.split('\t').collect();
    assert_eq!(fields[2], "1", "marker count for (0,2)");
}

#[test]
fn test_freq_exact_output() {
    let dir = test_dir().join("freq_exact");
    std::fs::create_dir_all(&dir).unwrap();

    let table = dir.join("markers.tsv");
    let mut f = std::fs::File::create(&table).unwrap();
    writeln!(f, "#Number of markers : 5").unwrap();
    writeln!(f, "id\tsequence\tind1\tind2\tind3\tind4").unwrap();
    writeln!(f, "0\tAAAA\t10\t10\t10\t10").unwrap(); // freq 4
    writeln!(f, "1\tBBBB\t10\t10\t0\t0").unwrap(); // freq 2
    writeln!(f, "2\tCCCC\t0\t0\t10\t10").unwrap(); // freq 2
    writeln!(f, "3\tDDDD\t10\t0\t0\t0").unwrap(); // freq 1
    writeln!(f, "4\tEEEE\t0\t0\t0\t0").unwrap(); // freq 0 (not counted)

    let output = dir.join("freq.tsv");
    rsx_core::commands::freq::run(&rsx_core::commands::freq::FreqParams {
        markers_table_path: table.to_str().unwrap().to_string(),
        output_file_path: output.to_str().unwrap().to_string(),
        min_depth: 1,
    })
    .unwrap();

    let content = std::fs::read_to_string(&output).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    // Skip comment + header
    assert_eq!(lines[2], "1\t1"); // 1 marker at freq 1
    assert_eq!(lines[3], "2\t2"); // 2 markers at freq 2
    assert_eq!(lines[4], "3\t0"); // 0 markers at freq 3
    assert_eq!(lines[5], "4\t1"); // 1 marker at freq 4
}

// === Cg float formatter precision tests ===

#[test]
fn test_cg_format_matches_cpp_g() {
    use rsx_core::stats::Cg;

    // These are exact values from the C++ radsex output
    assert_eq!(format!("{}", Cg(0.456056540250256)), "0.456057");
    assert_eq!(format!("{}", Cg(0.0000569941162332776)), "5.69941e-05");
    assert_eq!(format!("{}", Cg(1.0)), "1");
    assert_eq!(format!("{}", Cg(0.0)), "0");
    assert_eq!(format!("{}", Cg(-0.2)), "-0.2");
    assert_eq!(format!("{}", Cg(-1.0)), "-1");
    assert_eq!(format!("{}", Cg(0.0113988)), "0.0113988");
}
