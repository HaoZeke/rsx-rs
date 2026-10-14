// GPL-3.0-or-later
// End-to-end CLI coverage on the small panel under tests/fixtures.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn fixture(relative: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(relative)
}

fn output_of(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_rsx"))
        .args(args)
        .output()
        .expect("spawn rsx")
}

fn succeed(args: &[&str]) -> Output {
    let output = output_of(args);
    assert!(
        output.status.success(),
        "rsx {args:?} failed with {}\nstdout:\n{}\nstderr:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn fail(args: &[&str]) {
    let output = output_of(args);
    assert!(
        !output.status.success(),
        "rsx {args:?} should fail\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

fn sequence_rows(body: &str) -> BTreeMap<String, Vec<u16>> {
    let mut rows = BTreeMap::new();
    for line in body.lines() {
        if line.is_empty() || line.starts_with('#') || line.starts_with("id\t") {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let depths = fields[2..]
            .iter()
            .map(|value| value.parse::<u16>().expect("depth"))
            .collect();
        rows.insert(fields[1].to_string(), depths);
    }
    rows
}

fn assert_marker_count(body: &str, rows: &BTreeMap<String, Vec<u16>>) {
    let count_line = body.lines().next().expect("marker count");
    assert!(count_line.starts_with("#Number of markers"), "{count_line}");
    let claimed: usize = count_line
        .split(':')
        .nth(1)
        .expect("count field")
        .trim()
        .parse()
        .expect("count integer");
    assert_eq!(rows.len(), claimed, "{body}");
}

fn path_str(path: &Path) -> &str {
    path.to_str().expect("utf-8 path")
}

#[test]
fn help_and_version_succeed() {
    let help = succeed(&["--help"]);
    let help_text = String::from_utf8(help.stdout).unwrap();
    assert!(help_text.contains("process"));
    assert!(help_text.contains("signif"));

    let freq_help = succeed(&["freq", "--help"]);
    let freq_text = String::from_utf8(freq_help.stdout).unwrap();
    assert!(freq_text.contains("markers-table"));

    let version = succeed(&["--version"]);
    assert_eq!(
        String::from_utf8(version.stdout).unwrap(),
        format!("rsx {}\n", env!("CARGO_PKG_VERSION"))
    );
}

#[test]
fn missing_inputs_and_bad_groups_exit_nonzero() {
    let dir = tempfile::tempdir().unwrap();
    let output = dir.path().join("out.tsv");
    fail(&[
        "freq",
        "-t",
        path_str(&dir.path().join("missing.tsv")),
        "-o",
        path_str(&output),
    ]);
    fail(&[
        "process",
        "-i",
        path_str(dir.path()),
        "-o",
        path_str(&output),
    ]);
    fail(&[
        "distrib",
        "-t",
        path_str(&fixture("panel/popmap.tsv")),
        "-p",
        path_str(&fixture("panel/popmap.tsv")),
        "-o",
        path_str(&output),
        "-G",
        "M",
    ]);
}

#[test]
fn panel_pipeline_matches_the_small_fixture() {
    let dir = tempfile::tempdir().unwrap();
    let markers = dir.path().join("markers.tsv");
    let freq = dir.path().join("freq.tsv");
    let depth = dir.path().join("depth.tsv");
    let distrib = dir.path().join("distrib.tsv");
    let signif = dir.path().join("signif.tsv");
    let subset = dir.path().join("subset.tsv");
    let triage = dir.path().join("triage.tsv");
    let pca = dir.path().join("pca");
    let merged = dir.path().join("merged.tsv");
    let reads = fixture("panel/reads");
    let popmap = fixture("panel/popmap.tsv");

    succeed(&[
        "process",
        "-i",
        path_str(&reads),
        "-o",
        path_str(&markers),
        "-d",
        "1",
        "-T",
        "1",
    ]);
    let processed = std::fs::read_to_string(&markers).unwrap();
    let header = processed.lines().nth(1).unwrap();
    assert_eq!(header, "id\tsequence\tf1\tf2\tf3\tf4\tm1\tm2\tm3\tm4");
    let rows = sequence_rows(&processed);
    assert_marker_count(&processed, &rows);
    assert_eq!(
        rows.get("ATCGATCGATCG").map(Vec::as_slice),
        Some([1, 1, 1, 1, 1, 1, 1, 1].as_slice())
    );
    assert_eq!(
        rows.get("GGGGCCCC").map(Vec::as_slice),
        Some([0, 0, 0, 0, 1, 1, 1, 1].as_slice())
    );
    assert_eq!(
        rows.get("TTTTAAAA").map(Vec::as_slice),
        Some([1, 1, 1, 1, 0, 0, 0, 0].as_slice())
    );

    succeed(&[
        "freq",
        "-t",
        path_str(&markers),
        "-o",
        path_str(&freq),
        "-d",
        "1",
    ]);
    let freq_body = std::fs::read_to_string(&freq).unwrap();
    assert!(freq_body.starts_with("#source:rsx-freq;min_depth:1\n"));
    let frequencies: BTreeMap<usize, u32> = freq_body
        .lines()
        .skip(2)
        .map(|line| {
            let (frequency, count) = line.split_once('\t').unwrap();
            (frequency.parse().unwrap(), count.parse().unwrap())
        })
        .collect();
    assert_eq!(frequencies.get(&4), Some(&2));
    assert_eq!(frequencies.get(&8), Some(&1));
    assert!(frequencies.iter().all(|(frequency, count)| {
        (*frequency == 4 && *count == 2) || (*frequency == 8 && *count == 1) || *count == 0
    }));

    succeed(&[
        "depth",
        "-t",
        path_str(&markers),
        "-p",
        path_str(&popmap),
        "-o",
        path_str(&depth),
        "--streaming-mode",
        "memory",
    ]);
    let depth_body = std::fs::read_to_string(&depth).unwrap();
    for name in ["f1", "f2", "f3", "f4"] {
        assert!(
            depth_body.contains(&format!("{name}\tF\t2\t2\t1\t1\t1\t1\t1\n")),
            "{depth_body}"
        );
    }
    for name in ["m1", "m2", "m3", "m4"] {
        assert!(
            depth_body.contains(&format!("{name}\tM\t2\t2\t1\t1\t1\t1\t1\n")),
            "{depth_body}"
        );
    }

    succeed(&[
        "distrib",
        "-t",
        path_str(&markers),
        "-p",
        path_str(&popmap),
        "-o",
        path_str(&distrib),
        "-G",
        "M,F",
        "--correction",
        "none",
    ]);
    let distrib_body = std::fs::read_to_string(&distrib).unwrap();
    assert!(distrib_body.contains("#source:rsx-distrib"));
    assert!(distrib_body.contains("\n4\t0\t1\t"), "{distrib_body}");
    assert!(distrib_body.contains("\n0\t4\t1\t"), "{distrib_body}");
    assert!(distrib_body.contains("\n4\t4\t1\t"), "{distrib_body}");

    succeed(&[
        "signif",
        "-t",
        path_str(&markers),
        "-p",
        path_str(&popmap),
        "-o",
        path_str(&signif),
        "-G",
        "M,F",
        "--correction",
        "none",
        "--test",
        "fisher",
    ]);
    let signif_body = std::fs::read_to_string(&signif).unwrap();
    assert!(signif_body.contains("#source:rsx-signif"));
    assert!(signif_body.contains("GGGGCCCC"), "{signif_body}");
    assert!(signif_body.contains("TTTTAAAA"), "{signif_body}");
    assert!(
        !signif_body
            .lines()
            .filter(|line| !line.starts_with('#'))
            .any(|line| line.contains("ATCGATCGATCG")),
        "{signif_body}"
    );

    succeed(&[
        "subset",
        "-t",
        path_str(&markers),
        "-p",
        path_str(&popmap),
        "-o",
        path_str(&subset),
        "-G",
        "M,F",
        "-m",
        "4",
        "-N",
        "0",
    ]);
    let subset_body = std::fs::read_to_string(&subset).unwrap();
    assert!(subset_body.contains("#source:rsx-subset"), "{subset_body}");
    let subset_rows = sequence_rows(&subset_body);
    assert_eq!(subset_rows.len(), 1);
    assert!(subset_rows.contains_key("GGGGCCCC"));

    succeed(&[
        "triage",
        "-t",
        path_str(&markers),
        "-p",
        path_str(&popmap),
        "-o",
        path_str(&triage),
        "-G",
        "M,F",
    ]);
    let triage_body = std::fs::read_to_string(&triage).unwrap();
    assert!(triage_body.contains("#source:rsx-triage"));
    assert!(triage_body.contains("Candidate_Class"));

    succeed(&[
        "pca",
        "-t",
        path_str(&markers),
        "-o",
        path_str(&pca),
        "-r",
        "2",
    ]);
    let summary = std::fs::read_to_string(pca.join("summary.txt")).unwrap();
    assert!(summary.contains("Markers: 3"), "{summary}");
    assert!(summary.contains("Individuals: 8"), "{summary}");
    assert!(pca.join("eigenvalues.tsv").is_file());
    assert!(pca.join("loadings.tsv").is_file());
    assert!(pca.join("sample_scores.tsv").is_file());

    succeed(&[
        "merge",
        "-o",
        path_str(&merged),
        path_str(&fixture("merge/left.tsv")),
        path_str(&fixture("merge/right.tsv")),
    ]);
    let merged_body = std::fs::read_to_string(&merged).unwrap();
    let merged_rows = sequence_rows(&merged_body);
    assert_marker_count(&merged_body, &merged_rows);
    assert_eq!(
        merged_rows.get("AAAA").map(Vec::as_slice),
        Some([1, 2, 8].as_slice())
    );
    assert_eq!(
        merged_rows.get("CCCC").map(Vec::as_slice),
        Some([4, 0, 0].as_slice())
    );
    assert_eq!(
        merged_rows.get("GGGG").map(Vec::as_slice),
        Some([0, 0, 3].as_slice())
    );
}
