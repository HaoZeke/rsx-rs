use std::process::Command;

fn output(command: &str, arguments: &[&str]) -> String {
    Command::new(command)
        .args(arguments)
        .output()
        .ok()
        .filter(|result| result.status.success())
        .map(|result| String::from_utf8_lossy(&result.stdout).trim().to_owned())
        .filter(|value| !value.is_empty())
        .unwrap_or_else(|| "unknown".to_owned())
}

fn output_allow_empty(command: &str, arguments: &[&str]) -> Option<String> {
    Command::new(command)
        .args(arguments)
        .output()
        .ok()
        .filter(|result| result.status.success())
        .map(|result| String::from_utf8_lossy(&result.stdout).trim().to_owned())
}

fn main() {
    let commit = output("git", &["-C", "..", "rev-parse", "HEAD"]);
    let dirty = match output_allow_empty("git", &["-C", "..", "status", "--porcelain"]) {
        None => "unknown",
        Some(status) if status.is_empty() => "false",
        Some(_) => "true",
    };
    println!("cargo:rustc-env=RSX_GIT_COMMIT={commit}");
    println!("cargo:rustc-env=RSX_GIT_DIRTY={dirty}");
    println!(
        "cargo:rustc-env=RSX_RUSTC_VERSION={}",
        output("rustc", &["--version"])
    );
    println!(
        "cargo:rustc-env=RSX_BUILD_TARGET={}",
        std::env::var("TARGET").unwrap_or_else(|_| "unknown".to_owned())
    );
    println!("cargo:rerun-if-changed=../.git/HEAD");
    println!("cargo:rerun-if-changed=../.git/index");
}
