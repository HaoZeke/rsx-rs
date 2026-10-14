#[cfg(feature = "gen-header")]
use std::{env, path::PathBuf};

#[cfg(feature = "gen-header")]
fn track_rust_sources(path: &std::path::Path) {
    println!("cargo:rerun-if-changed={}", path.display());
    if !path.is_dir() {
        return;
    }
    for entry in std::fs::read_dir(path).expect("Unable to inspect Rust sources") {
        let entry = entry.expect("Unable to inspect Rust source entry");
        let entry_path = entry.path();
        if entry_path.is_dir() || entry_path.extension().is_some_and(|value| value == "rs") {
            track_rust_sources(&entry_path);
        }
    }
}

#[cfg(feature = "gen-header")]
fn generate_c_header(crate_dir: &str) {
    let output_dir = PathBuf::from(crate_dir).join("include");
    std::fs::create_dir_all(&output_dir).unwrap();

    let mut config =
        cbindgen::Config::from_file("cbindgen.toml").expect("Unable to find cbindgen.toml");

    let version = env::var("CARGO_PKG_VERSION").unwrap();
    let parts: Vec<&str> = version.split('.').collect();
    let version_block = format!(
        "#define RSX_VERSION \"{version}\"\n\
         #define RSX_VERSION_MAJOR {major}\n\
         #define RSX_VERSION_MINOR {minor}\n\
         #define RSX_VERSION_PATCH {patch}",
        version = version,
        major = parts[0],
        minor = parts[1],
        patch = parts[2],
    );
    config.after_includes = Some(match config.after_includes {
        Some(existing) => format!("{existing}\n\n{version_block}"),
        None => version_block,
    });

    cbindgen::Builder::new()
        .with_crate(crate_dir)
        .with_config(config)
        .generate()
        .expect("Unable to generate C bindings")
        .write_to_file(output_dir.join("rsx.h"));

    let path = output_dir.join("rsx.h");
    let content = std::fs::read_to_string(&path).unwrap();
    std::fs::write(&path, format!("{}\n", content.trim_end())).unwrap();
}

fn main() {
    println!("cargo::rustc-check-cfg=cfg(rsx_avx512_intrinsics)");
    // AVX-512 packed-double helpers in std::arch stabilized in 1.89.
    // The crate MSRV stays 1.85; the Gram prototype uses them only when the
    // compiler is new enough.
    if rustc_minor() >= 89 {
        println!("cargo:rustc-cfg=rsx_avx512_intrinsics");
    }

    #[cfg(feature = "gen-header")]
    {
        let crate_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
        let crate_path = PathBuf::from(&crate_dir);
        track_rust_sources(&crate_path.join("src"));
        println!(
            "cargo:rerun-if-changed={}",
            crate_path.join("cbindgen.toml").display()
        );
        println!(
            "cargo:rerun-if-changed={}",
            crate_path.join("Cargo.toml").display()
        );
        generate_c_header(&crate_dir);
    }
}

fn rustc_minor() -> u32 {
    let rustc = std::env::var("RUSTC").unwrap_or_else(|_| "rustc".into());
    let Ok(output) = std::process::Command::new(rustc).arg("--version").output() else {
        return 0;
    };
    let text = String::from_utf8_lossy(&output.stdout);
    let Some(version) = text.split_whitespace().nth(1) else {
        return 0;
    };
    version
        .split('.')
        .nth(1)
        .and_then(|minor| minor.parse().ok())
        .unwrap_or(0)
}
