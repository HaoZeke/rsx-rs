// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Sequence file reader wrapping needletail for FASTQ/FASTA (optionally gzipped).

use std::path::{Path, PathBuf};

/// Supported file extensions for demultiplexed read files.
const EXTENSIONS: &[&str] = &[
    ".fq",
    ".fq.gz",
    ".fastq",
    ".fastq.gz",
    ".fasta",
    ".fasta.gz",
    ".fa",
    ".fa.gz",
    ".fna",
    ".fna.gz",
];

/// An input file with its individual name derived from the filename.
#[derive(Debug, Clone)]
pub struct InputFile {
    pub path: PathBuf,
    pub individual_name: String,
}

/// Scan a directory for supported sequence files and extract individual names.
pub fn get_input_files(dir: &Path) -> std::io::Result<Vec<InputFile>> {
    let mut files = Vec::new();

    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        if !path.is_file() {
            continue;
        }

        let filename = match path.file_name().and_then(|n| n.to_str()) {
            Some(n) => n.to_string(),
            None => continue,
        };

        // Check if the file has a supported extension
        let mut matched_ext = None;
        for ext in EXTENSIONS {
            if filename.ends_with(ext) {
                matched_ext = Some(*ext);
                break;
            }
        }

        if let Some(ext) = matched_ext {
            // Individual name = filename without the matched extension
            let individual_name = filename[..filename.len() - ext.len()].to_string();
            files.push(InputFile {
                path,
                individual_name,
            });
        }
    }

    files.sort_by(|a, b| a.individual_name.cmp(&b.individual_name));

    if files.is_empty() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            format!("No valid input files found in {}", dir.display()),
        ));
    }

    log::info!("Found {} reads files", files.len());
    Ok(files)
}

/// Count occurrences of each unique sequence in a single file.
/// Returns a map of sequence -> count.
pub fn count_sequences(
    path: &Path,
) -> Result<ahash::AHashMap<Vec<u8>, u16>, Box<dyn std::error::Error + Send + Sync>> {
    use needletail::parse_fastx_file;

    let mut counts: ahash::AHashMap<Vec<u8>, u16> = ahash::AHashMap::new();
    let mut reader = parse_fastx_file(path)?;

    while let Some(record) = reader.next() {
        let record = record?;
        let seq = record.seq();
        let entry = counts.entry(seq.to_vec()).or_insert(0);
        *entry = entry.saturating_add(1);
    }

    Ok(counts)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extension_matching() {
        // Verify extension list is correct
        assert!(EXTENSIONS.contains(&".fq.gz"));
        assert!(EXTENSIONS.contains(&".fastq"));
        assert!(!EXTENSIONS.contains(&".bam"));
    }
}
