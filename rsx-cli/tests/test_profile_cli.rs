// GPL-3.0-or-later

use std::fs;
use std::io::Read;
use std::process::Command;

use rsx_core::run_profile::{CommandProfile, RunProfile};

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
