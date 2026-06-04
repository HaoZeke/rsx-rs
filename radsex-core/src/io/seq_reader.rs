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
    pack_2bit_vec(seq)
}

#[inline(always)]
fn pack_2bit_vec(seq: &[u8]) -> Vec<u8> {
    let n = seq.len();
    let packed_len = n.div_ceil(4);
    let mut packed = vec![0u8; packed_len + 1];
    packed[0] = (n % 4) as u8;

    for (i, &base) in seq.iter().enumerate() {
        let byte_idx = 1 + i / 4;
        let bit_pos = 6 - 2 * (i % 4);
        packed[byte_idx] |= base_to_2bit(base) << bit_pos;
    }
    packed
}

#[inline(always)]
fn base_to_2bit(base: u8) -> u8 {
    match base {
        b'A' | b'a' => 0b00,
        b'C' | b'c' => 0b01,
        b'G' | b'g' => 0b10,
        b'T' | b't' => 0b11,
        _ => {
            // Ambiguous bases use the same two-bit code path as the byte-key API.
            0b00
        }
    }
}

#[inline(always)]
fn pack_2bit_key(seq: &[u8]) -> PackedDnaKey {
    let n = seq.len();
    let packed_len = n.div_ceil(4);
    if packed_len + 1 > PACKED_DNA_INLINE_CAPACITY {
        return PackedDnaKey::Heap(pack_2bit_vec(seq));
    }

    let mut packed = [0u8; 48];
    packed[0] = (n % 4) as u8;

    for (i, &base) in seq.iter().enumerate() {
        let byte_idx = 1 + i / 4;
        let bit_pos = 6 - 2 * (i % 4);
        packed[byte_idx] |= base_to_2bit(base) << bit_pos;
    }

    let used = packed_len + 1;
    PackedDnaKey::Inline(PackedDna {
        data: packed,
        len: used as u8,
    })
}

const PACKED_DNA_INLINE_CAPACITY: usize = 48;

/// Compact stack-allocated key for 2-bit packed DNA (up to 188 bp / 47 bytes packed + len).
/// Used in the hot per-file counting map to avoid per-read heap allocation.
/// Only the used prefix participates in Hash/Eq. Converts to Vec<u8> for
/// public APIs and serialization boundaries.
#[derive(Clone, Copy)]
pub(crate) struct PackedDna {
    data: [u8; 48],
    len: u8,
}

impl PackedDna {
    #[inline(always)]
    fn as_slice(&self) -> &[u8] {
        &self.data[..self.len as usize]
    }
}

impl PartialEq for PackedDna {
    #[inline(always)]
    fn eq(&self, other: &Self) -> bool {
        self.as_slice() == other.as_slice()
    }
}

impl Eq for PackedDna {}

impl std::hash::Hash for PackedDna {
    #[inline(always)]
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.as_slice().hash(state);
    }
}

impl From<PackedDna> for Vec<u8> {
    fn from(p: PackedDna) -> Self {
        p.as_slice().to_vec()
    }
}

#[derive(Clone)]
pub(crate) enum PackedDnaKey {
    Inline(PackedDna),
    Heap(Vec<u8>),
}

impl PackedDnaKey {
    #[inline(always)]
    pub(crate) fn as_slice(&self) -> &[u8] {
        match self {
            PackedDnaKey::Inline(packed) => packed.as_slice(),
            PackedDnaKey::Heap(packed) => packed,
        }
    }
}

impl PartialEq for PackedDnaKey {
    #[inline(always)]
    fn eq(&self, other: &Self) -> bool {
        self.as_slice() == other.as_slice()
    }
}

impl Eq for PackedDnaKey {}

impl std::hash::Hash for PackedDnaKey {
    #[inline(always)]
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.as_slice().hash(state);
    }
}

impl From<PackedDnaKey> for Vec<u8> {
    fn from(p: PackedDnaKey) -> Self {
        match p {
            PackedDnaKey::Inline(packed) => packed.into(),
            PackedDnaKey::Heap(packed) => packed,
        }
    }
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

/// Count occurrences of each unique sequence in a single file using compact
/// internal keys.
pub(crate) fn count_sequences_packed(
    path: &Path,
) -> Result<ahash::AHashMap<PackedDnaKey, u16>, Box<dyn std::error::Error + Send + Sync>> {
    use needletail::parse_fastx_file;

    // Use stack-allocated keys for common read lengths and heap-backed keys for
    // longer reads. This preserves exact sequence identity for every input.
    let mut counts: ahash::AHashMap<PackedDnaKey, u16> = ahash::AHashMap::new();
    let mut reader = parse_fastx_file(path)?;

    while let Some(record) = reader.next() {
        let record = record?;
        let seq = record.seq();
        let packed = pack_2bit_key(&seq);
        let entry = counts.entry(packed).or_insert(0);
        *entry = entry.saturating_add(1);
    }

    Ok(counts)
}

/// Count occurrences of each unique sequence in a single file.
/// Uses 2-bit packed DNA keys for 4x memory reduction.
/// Returns a map of packed_sequence -> count.
pub fn count_sequences(
    path: &Path,
) -> Result<ahash::AHashMap<Vec<u8>, u16>, Box<dyn std::error::Error + Send + Sync>> {
    Ok(count_sequences_packed(path)?
        .into_iter()
        .map(|(k, v)| (k.into(), v))
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_extension_matching() {
        assert!(EXTENSIONS.contains(&".fq.gz"));
        assert!(EXTENSIONS.contains(&".fastq"));
        assert!(!EXTENSIONS.contains(&".bam"));
    }

    #[test]
    fn test_pack_unpack_roundtrip() {
        let long = vec![b'T'; 260];
        let seqs: [&[u8]; 8] = [
            b"ATCG".as_slice(),
            b"AAAA",
            b"TTTTCCCCGGGGAAAA",
            b"A",
            b"AT",
            b"ATC",
            b"ATCGATCGATCGATCGATCGATCGATCGATCG", // 32bp
            long.as_slice(),
        ];
        for seq in seqs {
            let packed = pack_2bit(seq);
            let unpacked = unpack_2bit(&packed);
            assert_eq!(
                unpacked.as_slice(),
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

    #[test]
    fn test_count_sequences_keeps_distinct_long_reads() {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        let long_a = vec![b'A'; 260];
        let mut long_b = vec![b'A'; 260];
        long_b[220] = b'T';

        writeln!(f, ">read1").unwrap();
        writeln!(f, "{}", std::str::from_utf8(&long_a).unwrap()).unwrap();
        writeln!(f, ">read2").unwrap();
        writeln!(f, "{}", std::str::from_utf8(&long_b).unwrap()).unwrap();

        let counts = count_sequences(f.path()).unwrap();
        assert_eq!(counts.len(), 2);
        assert_eq!(counts.get(&pack_2bit(&long_a)), Some(&1));
        assert_eq!(counts.get(&pack_2bit(&long_b)), Some(&1));
    }

    #[test]
    fn test_count_sequences_packed_uses_inline_short_read_key() {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, ">read1").unwrap();
        writeln!(f, "ATCGATCG").unwrap();
        writeln!(f, ">read2").unwrap();
        writeln!(f, "ATCGATCG").unwrap();

        let counts = count_sequences_packed(f.path()).unwrap();
        assert_eq!(counts.len(), 1);
        let (key, count) = counts.iter().next().unwrap();
        assert!(matches!(key, PackedDnaKey::Inline(_)));
        assert_eq!(*count, 2);
    }
}
