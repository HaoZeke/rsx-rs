// GPL-3.0-or-later

use std::fs;
use std::io::Read;
use std::process::Command;
use std::thread;
use std::time::Duration;

use rsx_core::run_profile::{CommandProfile, RunProfile};
use rsx_core::stats::{bayes_factor_2x2, bayes_factor_2x2_with_model};
use sha2::{Digest, Sha256};

fn write_markers_table(path: &std::path::Path) {
    fs::write(
        path,
        "#Number of markers : 2\nid\tsequence\tind1\tind2\n0\tACGT\t4\t0\n1\tTGCA\t2\t3\n",
    )
    .unwrap();
}

fn read_archive_member(path: &std::path::Path, name: &str) -> Vec<u8> {
    let mut archive = zip::ZipArchive::new(fs::File::open(path).unwrap()).unwrap();
    let mut contents = Vec::new();
    archive
        .by_name(name)
        .unwrap()
        .read_to_end(&mut contents)
        .unwrap();
    contents
}

fn assert_archive_checksums(path: &std::path::Path) {
    let checksum_manifest = String::from_utf8(read_archive_member(path, "SHA256SUMS")).unwrap();
    for line in checksum_manifest.lines() {
        let (expected, name) = line.split_once("  ").unwrap();
        let actual = format!("{:x}", Sha256::digest(read_archive_member(path, name)));
        assert_eq!(actual, expected, "checksum mismatch for {name}");
    }
}

#[test]
fn triage_profile_priors_reach_the_runtime_calculation() {
    let directory = tempfile::tempdir().unwrap();
    let profile_path = directory.path().join("input.toml");
    let hydrated_path = directory.path().join("hydrated.toml");
    let markers_path = directory.path().join("markers.tsv");
    let popmap_path = directory.path().join("popmap.tsv");
    let output_path = directory.path().join("triage.tsv");
    fs::write(
        &markers_path,
        "#Number of markers : 1\nid\tsequence\tind1\tind2\n0\tACGT\t4\t0\n",
    )
    .unwrap();
    fs::write(&popmap_path, "ind1\tM\nind2\tF\n").unwrap();
    fs::write(
        &profile_path,
        format!(
            r#"
schema_version = 1
profile_name = "custom-priors-v1"
write_hydrated_profile = "{}"

[run]
command = "triage"
markers_table = "{}"
popmap = "{}"
output_file = "{}"
min_depth = 1
groups = ["M", "F"]
signif_threshold = 0.05
posterior_threshold = 0.01
bayes_factor_threshold = 0.01

[run.bayes_model]
linkage_prior = 0.5
linked_prevalence = 0.9
null_prevalence = 0.5
group1_linked_weight = 0.5

[run.bayes_model.bayes_factor.alternative_group1]
alpha = 8.0
beta = 2.0

[run.bayes_model.bayes_factor.alternative_group2]
alpha = 2.0
beta = 8.0

[run.bayes_model.bayes_factor.null]
alpha = 10.0
beta = 10.0
"#,
            hydrated_path.display(),
            markers_path.display(),
            popmap_path.display(),
            output_path.display(),
        ),
    )
    .unwrap();

    let status = Command::new(env!("CARGO_BIN_EXE_rsx"))
        .arg("--profile")
        .arg(&profile_path)
        .status()
        .unwrap();
    assert!(status.success());

    let hydrated = RunProfile::parse_toml(&fs::read_to_string(&hydrated_path).unwrap()).unwrap();
    let model = match hydrated.run {
        CommandProfile::Triage(profile) => profile.bayes_model.to_runtime().unwrap(),
        _ => panic!("triage profile changed command variants"),
    };
    let expected = bayes_factor_2x2_with_model(1, 0, 1, 1, &model.bayes_factor).unwrap();
    let compatibility = bayes_factor_2x2(1, 0, 1, 1);
    assert_ne!(format!("{expected:.4}"), format!("{compatibility:.4}"));

    let output = fs::read_to_string(&output_path).unwrap();
    let row = output.lines().nth(2).unwrap();
    let observed: f64 = row.split('\t').nth(14).unwrap().parse().unwrap();
    assert_eq!(observed, format!("{expected:.4}").parse::<f64>().unwrap());
}

#[test]
fn failed_analysis_keeps_the_hydrated_profile_written_before_execution() {
    let directory = tempfile::tempdir().unwrap();
    let profile_path = directory.path().join("input.toml");
    let hydrated_path = directory.path().join("hydrated.toml");
    let archive_path = directory.path().join("reproducibility.zip");
    let missing_markers = directory.path().join("missing-markers.tsv");
    let output_path = directory.path().join("freq.tsv");
    fs::write(
        &profile_path,
        format!(
            r#"
schema_version = 1
profile_name = "failure-capture-v1"
write_hydrated_profile = "{}"
reproducibility_archive = "{}"

[run]
command = "freq"
markers_table = "{}"
output_file = "{}"
min_depth = 1
"#,
            hydrated_path.display(),
            archive_path.display(),
            missing_markers.display(),
            output_path.display(),
        ),
    )
    .unwrap();

    let status = Command::new(env!("CARGO_BIN_EXE_rsx"))
        .arg("--profile")
        .arg(&profile_path)
        .status()
        .unwrap();
    assert!(!status.success());

    let hydrated = RunProfile::parse_toml(&fs::read_to_string(&hydrated_path).unwrap()).unwrap();
    assert_eq!(hydrated.profile_name, "failure-capture-v1");
    match hydrated.run {
        CommandProfile::Freq(profile) => {
            assert_eq!(profile.markers_table, missing_markers.to_string_lossy());
            assert_eq!(profile.output_file, output_path.to_string_lossy());
            assert_eq!(profile.min_depth, 1);
        }
        _ => panic!("hydrated profile changed command variants"),
    }

    let mut archive = zip::ZipArchive::new(fs::File::open(&archive_path).unwrap()).unwrap();
    for member in [
        "profile.input.toml",
        "profile.hydrated.toml",
        "build-manifest.toml",
        "run-manifest.toml",
        "sbom.cdx.json",
        "Cargo.lock",
        "CITATION.cff",
        "LICENSE",
        "bin/rsx",
        "SHA256SUMS",
    ] {
        assert!(archive.by_name(member).is_ok(), "missing {member}");
    }
    let mut manifest = String::new();
    archive
        .by_name("run-manifest.toml")
        .unwrap()
        .read_to_string(&mut manifest)
        .unwrap();
    assert!(manifest.contains("status = \"failed\""));
    assert!(manifest.contains("created_before_execution = true"));
}

#[test]
fn invalid_profile_syntax_still_creates_a_resolution_failure_archive() {
    let directory = tempfile::tempdir().unwrap();
    let profile_path = directory.path().join("invalid.toml");
    let archive_path = directory.path().join("resolution-failure.zip");
    let invalid_profile = "schema_version = 1\nprofile_name = [\n";
    fs::write(&profile_path, invalid_profile).unwrap();

    let status = Command::new(env!("CARGO_BIN_EXE_rsx"))
        .arg("--profile")
        .arg(&profile_path)
        .arg("--reproducibility-archive")
        .arg(&archive_path)
        .status()
        .unwrap();
    assert!(!status.success());

    let mut archive = zip::ZipArchive::new(fs::File::open(&archive_path).unwrap()).unwrap();
    for member in [
        "profile.input.toml",
        "build-manifest.toml",
        "run-manifest.toml",
        "sbom.cdx.json",
        "bin/rsx",
        "SHA256SUMS",
    ] {
        assert!(archive.by_name(member).is_ok(), "missing {member}");
    }
    assert!(archive.by_name("profile.hydrated.toml").is_err());

    let mut stored_profile = String::new();
    archive
        .by_name("profile.input.toml")
        .unwrap()
        .read_to_string(&mut stored_profile)
        .unwrap();
    assert_eq!(stored_profile, invalid_profile);

    let mut manifest = String::new();
    archive
        .by_name("run-manifest.toml")
        .unwrap()
        .read_to_string(&mut manifest)
        .unwrap();
    assert!(manifest.contains("status = \"resolution-failed\""));
    assert!(manifest.contains("error_category = \"configuration\""));
}

#[test]
fn successful_analysis_creates_a_completed_repeatable_archive() {
    let directory = tempfile::tempdir().unwrap();
    let profile_path = directory.path().join("input.toml");
    let archive_path = directory.path().join("completed.zip");
    let markers_path = directory.path().join("markers.tsv");
    let output_path = directory.path().join("freq.tsv");
    write_markers_table(&markers_path);
    fs::write(
        &profile_path,
        format!(
            r#"
schema_version = 1
profile_name = "completed-v1"
reproducibility_archive = "{}"

[run]
command = "freq"
markers_table = "{}"
output_file = "{}"
min_depth = 1
"#,
            archive_path.display(),
            markers_path.display(),
            output_path.display(),
        ),
    )
    .unwrap();

    let run = || {
        Command::new(env!("CARGO_BIN_EXE_rsx"))
            .arg("--profile")
            .arg(&profile_path)
            .env("SOURCE_DATE_EPOCH", "0")
            .status()
            .unwrap()
    };
    assert!(run().success());
    let first = fs::read(&archive_path).unwrap();
    assert!(
        String::from_utf8(read_archive_member(&archive_path, "run-manifest.toml"))
            .unwrap()
            .contains("status = \"completed\"")
    );
    assert_archive_checksums(&archive_path);

    assert!(run().success());
    assert_eq!(fs::read(&archive_path).unwrap(), first);
}

#[cfg(unix)]
#[test]
fn unhandled_termination_leaves_the_started_archive() {
    let directory = tempfile::tempdir().unwrap();
    let profile_path = directory.path().join("input.toml");
    let archive_path = directory.path().join("interrupted.zip");
    let markers_path = directory.path().join("markers.tsv");
    let output_fifo = directory.path().join("blocked-output.tsv");
    write_markers_table(&markers_path);
    assert!(Command::new("mkfifo")
        .arg(&output_fifo)
        .status()
        .unwrap()
        .success());
    fs::write(
        &profile_path,
        format!(
            r#"
schema_version = 1
profile_name = "interrupted-v1"
reproducibility_archive = "{}"

[run]
command = "freq"
markers_table = "{}"
output_file = "{}"
min_depth = 1
"#,
            archive_path.display(),
            markers_path.display(),
            output_fifo.display(),
        ),
    )
    .unwrap();

    let mut child = Command::new(env!("CARGO_BIN_EXE_rsx"))
        .arg("--profile")
        .arg(&profile_path)
        .spawn()
        .unwrap();
    for _ in 0..400 {
        if archive_path.exists() {
            break;
        }
        thread::sleep(Duration::from_millis(50));
    }
    assert!(archive_path.exists(), "initial archive was not written");
    child.kill().unwrap();
    let status = child.wait().unwrap();
    assert!(!status.success());

    let manifest =
        String::from_utf8(read_archive_member(&archive_path, "run-manifest.toml")).unwrap();
    assert!(manifest.contains("status = \"started\""));
    assert!(manifest.contains("created_before_execution = true"));
    assert_archive_checksums(&archive_path);
}
