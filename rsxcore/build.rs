#[cfg(feature = "gen-header")]
use std::{env, path::PathBuf};

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
    #[cfg(feature = "gen-header")]
    {
        let crate_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
        generate_c_header(&crate_dir);
    }
}
