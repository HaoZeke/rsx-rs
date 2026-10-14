// GPL-3.0-or-later

use std::fs;
use std::process::Command;

fn rsx() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsx"))
}

fn stderr_of(command: &mut Command) -> (bool, String) {
    let output = command.output().expect("rsx binary runs");
    let text = String::from_utf8(output.stderr).expect("stderr is utf-8");
    (output.status.success(), text)
}

#[test]
fn help_names_the_commands_and_the_signif_columns() {
    let (ok, _) = stderr_of(rsx().arg("--help"));
    assert!(ok, "rsx --help exits 0");
    let help = String::from_utf8(rsx().arg("--help").output().unwrap().stdout).unwrap();
    assert!(help.contains("signif"), "{help}");
    assert!(help.contains("Examples:"), "{help}");
    assert!(help.contains("--profile"), "{help}");

    let signif =
        String::from_utf8(rsx().args(["signif", "--help"]).output().unwrap().stdout).unwrap();
    assert!(signif.contains("Bayes_Factor"), "{signif}");
    assert!(signif.contains("fdr"), "{signif}");
    assert!(signif.contains("--groups"), "{signif}");
}

#[test]
fn unknown_test_names_the_flag_and_the_accepted_values() {
    let (ok, err) = stderr_of(
        rsx()
            .args([
                "signif",
                "-t",
                "markers.tsv",
                "-p",
                "popmap.tsv",
                "-o",
                "out.tsv",
            ])
            .args(["--test", "anova"]),
    );
    assert!(!ok);
    assert!(err.contains("error:"), "{err}");
    assert!(err.contains("signif:"), "{err}");
    assert!(err.contains("anova"), "{err}");
    assert!(err.contains("chisq"), "{err}");
    assert!(err.contains("fisher"), "{err}");
    assert!(err.contains("gtest"), "{err}");
    assert!(!err.to_lowercase().contains("panic"), "{err}");
}

#[test]
fn a_missing_markers_table_names_the_path_and_the_expected_header() {
    let directory = tempfile::tempdir().unwrap();
    let missing = directory.path().join("markers.tsv");
    let output = directory.path().join("freq.tsv");
    let (ok, err) = stderr_of(
        rsx()
            .arg("freq")
            .arg("-t")
            .arg(&missing)
            .arg("-o")
            .arg(&output),
    );
    assert!(!ok);
    assert!(err.contains("cannot open markers table"), "{err}");
    assert!(err.contains(&missing.display().to_string()), "{err}");
    assert!(err.contains("id"), "{err}");
    assert!(!output.exists());
}

#[test]
fn one_group_name_explains_the_required_pair() {
    let directory = tempfile::tempdir().unwrap();
    let markers = directory.path().join("markers.tsv");
    let popmap = directory.path().join("popmap.tsv");
    let output = directory.path().join("signif.tsv");
    fs::write(
        &markers,
        "#Number of markers : 1\nid\tsequence\tind1\tind2\n0\tACGT\t4\t0\n",
    )
    .unwrap();
    fs::write(&popmap, "ind1\tM\nind2\tF\n").unwrap();
    let (ok, err) = stderr_of(
        rsx()
            .arg("signif")
            .arg("-t")
            .arg(&markers)
            .arg("-p")
            .arg(&popmap)
            .arg("-o")
            .arg(&output)
            .args(["--groups", "M"]),
    );
    assert!(!ok);
    assert!(err.contains("exactly two non-empty group names"), "{err}");
    assert!(err.contains("--groups M,F"), "{err}");
}

#[test]
fn a_zero_kmer_length_is_rejected_before_reading_inputs() {
    let (ok, err) = stderr_of(
        rsx()
            .args(["process", "-i", "does-not-need-to-exist", "-o", "out.tsv"])
            .args(["-k", "0"]),
    );
    assert!(!ok);
    assert!(err.contains("kmer_dedup must be at least 1"), "{err}");
    assert!(!err.to_lowercase().contains("panic"), "{err}");
}
