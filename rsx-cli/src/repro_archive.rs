// GPL-3.0-or-later

use std::collections::BTreeMap;
use std::error::Error;
use std::fs::{self, File, OpenOptions};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use rsx_core::run_profile::RunProfile;
use serde::Serialize;
use serde_json::json;
use sha2::{Digest, Sha256};
use zip::write::SimpleFileOptions;
use zip::{CompressionMethod, DateTime, ZipWriter};

const CARGO_LOCK: &[u8] = include_bytes!("../../Cargo.lock");
const PIXI_LOCK: &[u8] = include_bytes!("../../pixi.lock");
const CITATION: &[u8] = include_bytes!("../../CITATION.cff");
const LICENSE: &[u8] = include_bytes!("../../LICENSE");
const RUN_PROFILE_SCHEMA: &[u8] =
    include_bytes!("../../rsx-python/schema/run-profile-v1.schema.json");

#[derive(Clone, Copy, Debug)]
pub enum ArchiveStatus<'a> {
    Started,
    Completed,
    Failed(&'a str),
}

#[derive(Serialize)]
struct RunManifest {
    manifest_version: u32,
    status: String,
    created_before_execution: bool,
    profile_resolved: bool,
    command: Vec<String>,
    profile_name: String,
    profile_sha256: String,
    working_directory: String,
    environment: BTreeMap<String, String>,
    error_category: Option<String>,
    error_message: Option<String>,
}

#[derive(Serialize)]
struct BuildManifest {
    manifest_version: u32,
    rsx_version: &'static str,
    git_commit: &'static str,
    git_dirty: &'static str,
    rustc: &'static str,
    target: &'static str,
    enabled_features: Vec<&'static str>,
    executable_sha256: String,
    sbom_format: &'static str,
    run_profile_schema: &'static str,
}

pub fn write_archive(
    destination: &Path,
    input_profile: &str,
    hydrated: &RunProfile,
    status: ArchiveStatus<'_>,
) -> Result<(), Box<dyn Error>> {
    let hydrated_toml = hydrated.to_toml()?;

    let (status_name, error_category, error_message) = match status {
        ArchiveStatus::Started => ("started", None, None),
        ArchiveStatus::Completed => ("completed", None, None),
        ArchiveStatus::Failed(message) => (
            "failed",
            Some("analysis".to_owned()),
            Some(message.to_owned()),
        ),
    };
    let run_manifest = RunManifest {
        manifest_version: 1,
        status: status_name.to_owned(),
        created_before_execution: true,
        profile_resolved: true,
        command: std::env::args().collect(),
        profile_name: hydrated.profile_name.clone(),
        profile_sha256: sha256(hydrated_toml.as_bytes()),
        working_directory: std::env::current_dir()?.to_string_lossy().into_owned(),
        environment: allowed_environment(),
        error_category,
        error_message,
    };
    let mut members = common_members()?;
    members.insert("profile.hydrated.toml".into(), hydrated_toml.into_bytes());
    members.insert(
        "profile.input.toml".into(),
        input_profile.as_bytes().to_vec(),
    );
    members.insert(
        "run-manifest.toml".into(),
        toml::to_string_pretty(&run_manifest)?.into_bytes(),
    );
    finish_archive(destination, members)
}

pub fn write_resolution_failure(
    destination: &Path,
    input_profile: &str,
    error: &str,
) -> Result<(), Box<dyn Error>> {
    let run_manifest = RunManifest {
        manifest_version: 1,
        status: "resolution-failed".to_owned(),
        created_before_execution: true,
        profile_resolved: false,
        command: std::env::args().collect(),
        profile_name: "unresolved".to_owned(),
        profile_sha256: sha256(input_profile.as_bytes()),
        working_directory: std::env::current_dir()?.to_string_lossy().into_owned(),
        environment: allowed_environment(),
        error_category: Some("configuration".to_owned()),
        error_message: Some(error.to_owned()),
    };
    let mut members = common_members()?;
    members.insert(
        "profile.input.toml".into(),
        input_profile.as_bytes().to_vec(),
    );
    members.insert(
        "run-manifest.toml".into(),
        toml::to_string_pretty(&run_manifest)?.into_bytes(),
    );
    finish_archive(destination, members)
}

fn common_members() -> Result<BTreeMap<String, Vec<u8>>, Box<dyn Error>> {
    let executable = fs::read(std::env::current_exe()?)?;
    let build_manifest = BuildManifest {
        manifest_version: 1,
        rsx_version: env!("CARGO_PKG_VERSION"),
        git_commit: env!("RSX_GIT_COMMIT"),
        git_dirty: env!("RSX_GIT_DIRTY"),
        rustc: env!("RSX_RUSTC_VERSION"),
        target: env!("RSX_BUILD_TARGET"),
        enabled_features: enabled_features(),
        executable_sha256: sha256(&executable),
        sbom_format: "CycloneDX 1.5",
        run_profile_schema: "run-profile-v1.schema.json",
    };
    let mut members = BTreeMap::<String, Vec<u8>>::new();
    members.insert("Cargo.lock".into(), CARGO_LOCK.to_vec());
    members.insert("CITATION.cff".into(), CITATION.to_vec());
    members.insert("LICENSE".into(), LICENSE.to_vec());
    members.insert("README.txt".into(), readme().into_bytes());
    members.insert("bin/rsx".into(), executable);
    members.insert(
        "build-manifest.toml".into(),
        toml::to_string_pretty(&build_manifest)?.into_bytes(),
    );
    members.insert("pixi.lock".into(), PIXI_LOCK.to_vec());
    members.insert(
        "run-profile-v1.schema.json".into(),
        RUN_PROFILE_SCHEMA.to_vec(),
    );
    members.insert("sbom.cdx.json".into(), cyclonedx_sbom()?);
    Ok(members)
}

fn finish_archive(
    destination: &Path,
    mut members: BTreeMap<String, Vec<u8>>,
) -> Result<(), Box<dyn Error>> {
    let mut checksums = String::new();
    for (name, contents) in &members {
        checksums.push_str(&format!("{}  {name}\n", sha256(contents)));
    }
    members.insert("SHA256SUMS".into(), checksums.into_bytes());
    write_zip_atomic(destination, &members)
}

fn cyclonedx_sbom() -> Result<Vec<u8>, Box<dyn Error>> {
    let lock: toml::Value = toml::from_str(std::str::from_utf8(CARGO_LOCK)?)?;
    let components: Vec<_> = lock
        .get("package")
        .and_then(toml::Value::as_array)
        .into_iter()
        .flatten()
        .filter_map(|package| {
            let name = package.get("name")?.as_str()?;
            let version = package.get("version")?.as_str()?;
            Some(json!({
                "type": "library",
                "bom-ref": format!("pkg:cargo/{name}@{version}"),
                "name": name,
                "version": version,
                "purl": format!("pkg:cargo/{name}@{version}")
            }))
        })
        .collect();
    let document = json!({
        "bomFormat": "CycloneDX",
        "specVersion": "1.5",
        "serialNumber": cyclonedx_serial(CARGO_LOCK),
        "version": 1,
        "metadata": {
            "component": {
                "type": "application",
                "name": "rsx",
                "version": env!("CARGO_PKG_VERSION")
            }
        },
        "components": components
    });
    Ok(serde_json::to_vec_pretty(&document)?)
}

fn cyclonedx_serial(content: &[u8]) -> String {
    let digest = sha256(content);
    format!(
        "urn:uuid:{}-{}-{}-{}-{}",
        &digest[0..8],
        &digest[8..12],
        &digest[12..16],
        &digest[16..20],
        &digest[20..32]
    )
}

fn write_zip_atomic(
    destination: &Path,
    members: &BTreeMap<String, Vec<u8>>,
) -> Result<(), Box<dyn Error>> {
    let temporary = sibling_temporary(destination);
    let file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&temporary)?;
    let mut archive = ZipWriter::new(file);
    let options = SimpleFileOptions::default()
        .compression_method(CompressionMethod::Stored)
        .last_modified_time(DateTime::default())
        .unix_permissions(0o644);
    for (name, contents) in members {
        archive.start_file(name, options)?;
        archive.write_all(contents)?;
    }
    let file = archive.finish()?;
    file.sync_all()?;
    verify_zip(&temporary, members.len())?;
    fs::rename(temporary, destination)?;
    Ok(())
}

fn verify_zip(path: &Path, expected_members: usize) -> Result<(), Box<dyn Error>> {
    let mut archive = zip::ZipArchive::new(File::open(path)?)?;
    if archive.len() != expected_members {
        return Err(format!(
            "archive has {} members; expected {expected_members}",
            archive.len()
        )
        .into());
    }
    let mut checksums = String::new();
    archive
        .by_name("SHA256SUMS")?
        .read_to_string(&mut checksums)?;
    if checksums.is_empty() {
        return Err("archive checksum manifest is empty".into());
    }
    Ok(())
}

fn sibling_temporary(destination: &Path) -> PathBuf {
    static SEQUENCE: AtomicU64 = AtomicU64::new(0);
    let sequence = SEQUENCE.fetch_add(1, Ordering::Relaxed);
    let parent = destination.parent().unwrap_or_else(|| Path::new("."));
    let name = destination
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("reproducibility.zip");
    parent.join(format!(".{name}.rsx-{}-{sequence}.tmp", std::process::id()))
}

fn sha256(contents: &[u8]) -> String {
    format!("{:x}", Sha256::digest(contents))
}

fn allowed_environment() -> BTreeMap<String, String> {
    [
        "CUDA_VISIBLE_DEVICES",
        "RAYON_NUM_THREADS",
        "RUST_LOG",
        "SOURCE_DATE_EPOCH",
    ]
    .into_iter()
    .filter_map(|name| {
        std::env::var(name)
            .ok()
            .map(|value| (name.to_owned(), value))
    })
    .collect()
}

fn enabled_features() -> Vec<&'static str> {
    let mut features = Vec::new();
    if cfg!(feature = "cuda") {
        features.push("cuda");
    }
    if cfg!(feature = "map") {
        features.push("map");
    }
    if cfg!(feature = "mpi") {
        features.push("mpi");
    }
    if cfg!(feature = "parquet-io") {
        features.push("parquet-io");
    }
    features
}

fn readme() -> String {
    "rsx software reproducibility archive\n\nVerify SHA256SUMS, inspect run-manifest.toml, and replay profile.hydrated.toml with the bundled executable. Input datasets and result tables are identified by the profile but are not stored in this archive.\n".to_owned()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn atomic_writer_rejects_a_member_checksum_mismatch() {
        let directory = tempfile::tempdir().unwrap();
        let destination = directory.path().join("invalid.zip");
        let mut members = BTreeMap::new();
        members.insert("payload.txt".to_owned(), b"payload".to_vec());
        members.insert(
            "SHA256SUMS".to_owned(),
            format!("{}  payload.txt\n", "0".repeat(64)).into_bytes(),
        );

        let error = write_zip_atomic(&destination, &members).unwrap_err();
        assert!(error
            .to_string()
            .contains("checksum mismatch for payload.txt"));
        assert!(!destination.exists());
    }
}
