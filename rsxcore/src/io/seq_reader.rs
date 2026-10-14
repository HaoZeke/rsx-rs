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

/// Branchless base encoding. Ambiguous bases stay 0, matching the previous match.
static BASE_TO_2BIT: [u8; 256] = {
    let mut table = [0_u8; 256];
    table[b'C' as usize] = 0b01;
    table[b'c' as usize] = 0b01;
    table[b'G' as usize] = 0b10;
    table[b'g' as usize] = 0b10;
    table[b'T' as usize] = 0b11;
    table[b't' as usize] = 0b11;
    table
};

#[inline(always)]
fn pack_2bit_vec(seq: &[u8]) -> Vec<u8> {
    let packed_len = seq.len().div_ceil(4);
    let mut packed = vec![0_u8; packed_len + 1];
    pack_into(seq, &mut packed);
    packed
}

/// Pack `seq` into `packed`, which must hold `1 + ceil(len/4)` zeroed bytes.
#[inline(always)]
fn pack_into(seq: &[u8], packed: &mut [u8]) {
    let n = seq.len();
    packed[0] = (n % 4) as u8;
    let mut index = 0;
    let mut byte_idx = 1;
    while index + 4 <= n {
        packed[byte_idx] = (BASE_TO_2BIT[seq[index] as usize] << 6)
            | (BASE_TO_2BIT[seq[index + 1] as usize] << 4)
            | (BASE_TO_2BIT[seq[index + 2] as usize] << 2)
            | BASE_TO_2BIT[seq[index + 3] as usize];
        byte_idx += 1;
        index += 4;
    }
    let mut acc = 0_u8;
    let mut shift = 6_u32;
    while index < n {
        acc |= BASE_TO_2BIT[seq[index] as usize] << shift;
        shift -= 2;
        index += 1;
    }
    if n % 4 != 0 {
        packed[byte_idx] = acc;
    }
}

#[inline(always)]
fn pack_2bit_key(seq: &[u8]) -> PackedDnaKey {
    let n = seq.len();
    let packed_len = n.div_ceil(4);
    if packed_len + 1 > PACKED_DNA_INLINE_CAPACITY {
        return PackedDnaKey::Heap(pack_2bit_vec(seq));
    }

    let mut packed = [0_u8; 48];
    pack_into(seq, &mut packed);
    PackedDnaKey::Inline(PackedDna {
        data: packed,
        len: (packed_len + 1) as u8,
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
    let mut seq = Vec::new();
    unpack_2bit_into(packed, &mut seq);
    seq
}

/// Unpack into a reused buffer. `seq` is cleared first.
#[inline]
pub fn unpack_2bit_into(packed: &[u8], seq: &mut Vec<u8>) {
    seq.clear();
    if packed.is_empty() {
        return;
    }
    let remainder = packed[0] as usize;
    let data = &packed[1..];
    let n = if remainder == 0 {
        data.len() * 4
    } else {
        (data.len() - 1) * 4 + remainder
    };

    seq.reserve(n);
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
}

/// Count occurrences of each unique sequence in a single file using compact
/// internal keys.
pub(crate) fn count_sequences_packed(
    path: &Path,
) -> Result<ahash::AHashMap<PackedDnaKey, u16>, Box<dyn std::error::Error + Send + Sync>> {
    // Four-line FASTQ (plain or gzip) is the demultiplexed RAD layout. A
    // dedicated scanner avoids per-read record objects. FASTA, wrapped FASTQ,
    // and anything the scanner rejects fall back to needletail so the counts
    // stay identical.
    if let Some(counts) = try_count_fastq_4line(path)? {
        return Ok(counts);
    }
    count_sequences_needletail(path)
}

fn count_sequences_needletail(
    path: &Path,
) -> Result<ahash::AHashMap<PackedDnaKey, u16>, Box<dyn std::error::Error + Send + Sync>> {
    use needletail::parse_fastx_file;

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

/// `Some(true)` for gzip FASTQ, `Some(false)` for plain FASTQ, `None` otherwise.
fn fastq_gzip_kind(path: &Path) -> Option<bool> {
    let name = path.file_name()?.to_str()?.to_ascii_lowercase();
    if name.ends_with(".fastq.gz") || name.ends_with(".fq.gz") {
        Some(true)
    } else if name.ends_with(".fastq") || name.ends_with(".fq") {
        Some(false)
    } else {
        None
    }
}

fn try_count_fastq_4line(
    path: &Path,
) -> Result<Option<ahash::AHashMap<PackedDnaKey, u16>>, Box<dyn std::error::Error + Send + Sync>> {
    let Some(gzip) = fastq_gzip_kind(path) else {
        return Ok(None);
    };
    let file = std::fs::File::open(path)?;
    if gzip {
        let decoder = flate2::read::GzDecoder::new(file);
        let reader = std::io::BufReader::with_capacity(1 << 20, decoder);
        scan_fastq_4line(reader)
    } else {
        let reader = std::io::BufReader::with_capacity(1 << 20, file);
        scan_fastq_4line(reader)
    }
}

fn scan_fastq_4line<R: std::io::BufRead>(
    mut reader: R,
) -> Result<Option<ahash::AHashMap<PackedDnaKey, u16>>, Box<dyn std::error::Error + Send + Sync>> {
    let mut counts: ahash::AHashMap<PackedDnaKey, u16> = ahash::AHashMap::new();
    let mut pending: Vec<u8> = Vec::new();
    let mut stage = 0_u8;
    let mut seq_len = 0_usize;
    loop {
        let step = {
            let buf = reader.fill_buf()?;
            if buf.is_empty() {
                None
            } else {
                let mut start = 0;
                let mut invalid = false;
                while let Some(rel) = memchr::memchr(b'\n', &buf[start..]) {
                    let end = start + rel;
                    let mut line_end = end;
                    if line_end > start && buf[line_end - 1] == b'\r' {
                        line_end -= 1;
                    }
                    let valid = if pending.is_empty() {
                        accept_fastq_line(
                            &buf[start..line_end],
                            &mut stage,
                            &mut seq_len,
                            &mut counts,
                        )
                    } else {
                        pending.extend_from_slice(&buf[start..line_end]);
                        let valid =
                            accept_fastq_line(&pending, &mut stage, &mut seq_len, &mut counts);
                        pending.clear();
                        valid
                    };
                    if !valid {
                        invalid = true;
                        break;
                    }
                    start = end + 1;
                }
                if !invalid && start < buf.len() {
                    pending.extend_from_slice(&buf[start..]);
                }
                Some((buf.len(), invalid))
            }
        };
        let Some((consumed, invalid)) = step else {
            break;
        };
        reader.consume(consumed);
        if invalid {
            return Ok(None);
        }
    }
    if !pending.is_empty() && !accept_fastq_line(&pending, &mut stage, &mut seq_len, &mut counts) {
        return Ok(None);
    }
    if stage != 0 {
        return Ok(None);
    }
    Ok(Some(counts))
}

fn accept_fastq_line(
    line: &[u8],
    stage: &mut u8,
    seq_len: &mut usize,
    counts: &mut ahash::AHashMap<PackedDnaKey, u16>,
) -> bool {
    match *stage {
        0 => {
            if line.first() != Some(&b'@') {
                return false;
            }
            *stage = 1;
        }
        1 => {
            *seq_len = line.len();
            let packed = pack_2bit_key(line);
            let entry = counts.entry(packed).or_insert(0);
            *entry = entry.saturating_add(1);
            *stage = 2;
        }
        2 => {
            if line.first() != Some(&b'+') {
                return false;
            }
            *stage = 3;
        }
        _ => {
            if line.len() != *seq_len {
                return false;
            }
            *stage = 0;
        }
    }
    true
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

    #[test]
    fn fastq_4line_matches_needletail_including_gzip() {
        let dir = tempfile::tempdir().unwrap();
        let plain = dir.path().join("sample.fq");
        let mut body = Vec::new();
        body.extend_from_slice(b"@r1\nACGT\n+\nIIII\n@r2\nacgt\n+\nIIII\n@r3\nACGT\n+\nIIII\n");
        std::fs::write(&plain, &body).unwrap();
        let gz_path = dir.path().join("sample.fq.gz");
        let mut encoder = flate2::write::GzEncoder::new(
            std::fs::File::create(&gz_path).unwrap(),
            flate2::Compression::fast(),
        );
        encoder.write_all(&body).unwrap();
        encoder.finish().unwrap();

        for path in [&plain, &gz_path] {
            let fast = try_count_fastq_4line(path).unwrap().expect("4-line fastq");
            let slow = count_sequences_needletail(path).unwrap();
            assert_eq!(fast.len(), slow.len());
            for (key, count) in slow {
                assert_eq!(fast.get(&key), Some(&count));
            }
        }
    }

    #[test]
    fn non_fastq_extension_and_fasta_named_fq_fall_back() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("sample.fa");
        std::fs::write(&fasta, b">r1\nACGTACGT\n").unwrap();
        assert!(try_count_fastq_4line(&fasta).unwrap().is_none());
        let counts = count_sequences(&fasta).unwrap();
        assert_eq!(counts.get(&pack_2bit(b"ACGTACGT")), Some(&1));

        // A .fq path whose first record is FASTA is rejected by the 4-line
        // scanner and counted by needletail instead.
        let disguised = dir.path().join("sample.fq");
        std::fs::write(&disguised, b">r1\nTTTT\n").unwrap();
        assert!(try_count_fastq_4line(&disguised).unwrap().is_none());
        let counts = count_sequences(&disguised).unwrap();
        assert_eq!(counts.get(&pack_2bit(b"TTTT")), Some(&1));
    }
}
