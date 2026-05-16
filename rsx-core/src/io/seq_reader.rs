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

/// Pack a DNA sequence into 2-bit encoding: A=00, C=01, G=10, T=11.
/// 4 bases per byte, big-endian within each byte.
/// Returns the packed bytes.
#[inline(always)]
pub fn pack_2bit(seq: &[u8]) -> Vec<u8> {
    let n = seq.len();
    let packed_len = n.div_ceil(4);
    let mut packed = vec![0u8; packed_len + 1]; // +1 byte for length % 4
    packed[0] = (n % 4) as u8; // store remainder for unpack

    for (i, &base) in seq.iter().enumerate() {
        let bits = match base {
            b'A' | b'a' => 0b00,
            b'C' | b'c' => 0b01,
            b'G' | b'g' => 0b10,
            b'T' | b't' => 0b11,
            _ => 0b00, // N and others map to A
        };
        let byte_idx = 1 + i / 4;
        let bit_pos = 6 - 2 * (i % 4); // 6, 4, 2, 0
        packed[byte_idx] |= bits << bit_pos;
    }
    packed
}

/// Unpack a 2-bit encoded DNA sequence back to ASCII.
#[inline(always)]
pub fn unpack_2bit(packed: &[u8]) -> Vec<u8> {
    if packed.is_empty() {
        return Vec::new();
    }
    let remainder = packed[0] as usize;
    let data = &packed[1..];
    let n = if remainder == 0 {
        data.len() * 4
    } else {
        (data.len() - 1) * 4 + remainder
    };

    let mut seq = Vec::with_capacity(n);
    for i in 0..n {
        let byte_idx = i / 4;
        let bit_pos = 6 - 2 * (i % 4);
        let bits = (data[byte_idx] >> bit_pos) & 0b11;
        seq.push(match bits {
            0b00 => b'A',
            0b01 => b'C',
            0b10 => b'G',
            0b11 => b'T',
            _ => unreachable!(),
        });
    }
    seq
}

/// Count occurrences of each unique sequence in a single file.
/// Uses 2-bit packed DNA keys for 4x memory reduction.
/// Returns a map of packed_sequence -> count.
pub fn count_sequences(
    path: &Path,
) -> Result<ahash::AHashMap<Vec<u8>, u16>, Box<dyn std::error::Error + Send + Sync>> {
    use needletail::parse_fastx_file;

    let mut counts: ahash::AHashMap<Vec<u8>, u16> = ahash::AHashMap::new();
    let mut reader = parse_fastx_file(path)?;

    while let Some(record) = reader.next() {
        let record = record?;
        let seq = record.seq();
        let packed = pack_2bit(&seq);
        let entry = counts.entry(packed).or_insert(0);
        *entry = entry.saturating_add(1); // already saturating; high-depth tags cap at 65535 per individual
    }

    Ok(counts)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extension_matching() {
        assert!(EXTENSIONS.contains(&".fq.gz"));
        assert!(EXTENSIONS.contains(&".fastq"));
        assert!(!EXTENSIONS.contains(&".bam"));
    }

    #[test]
    fn test_pack_unpack_roundtrip() {
        let seqs = [
            b"ATCG".as_slice(),
            b"AAAA",
            b"TTTTCCCCGGGGAAAA",
            b"A",
            b"AT",
            b"ATC",
            b"ATCGATCGATCGATCGATCGATCGATCGATCG", // 32bp
        ];
        for seq in &seqs {
            let packed = pack_2bit(seq);
            let unpacked = unpack_2bit(&packed);
            assert_eq!(
                &unpacked,
                seq,
                "roundtrip failed for {}",
                std::str::from_utf8(seq).unwrap()
            );
        }
    }

    #[test]
    fn test_pack_2bit_compression() {
        // 100bp -> should be ~26 bytes (1 header + 25 data)
        let seq = vec![b'A'; 100];
        let packed = pack_2bit(&seq);
        assert_eq!(packed.len(), 26); // 1 + ceil(100/4)
        assert!(packed.len() < seq.len()); // 4x compression
    }

    #[test]
    fn test_pack_2bit_case_insensitive() {
        let upper = pack_2bit(b"ATCG");
        let lower = pack_2bit(b"atcg");
        assert_eq!(upper, lower);
    }
}
