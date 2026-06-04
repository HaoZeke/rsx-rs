// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! TSV table I/O for markers depth tables.

use std::io::{self, BufRead};
use std::path::Path;

/// Header information parsed from a markers depth table.
#[derive(Debug, Clone)]
pub struct TableHeader {
    /// Total number of markers (from the comment line).
    pub n_markers: u64,
    /// Number of individuals (columns - 2: id and sequence).
    pub n_individuals: u16,
    /// All column names from the header line.
    pub columns: Vec<String>,
}

impl TableHeader {
    /// Parse the header from a markers depth table file.
    /// The file starts with an optional comment line `#Number of markers : N`
    /// followed by a tab-separated header line `id\tsequence\tind1\t...\tindN`.
    pub fn from_file(path: &Path) -> io::Result<Self> {
        let file = std::fs::File::open(path)?;
        let mut reader = io::BufReader::new(file);
        let mut line = String::new();

        // Read first line
        reader.read_line(&mut line)?;
        let mut n_markers = 0u64;

        if line.starts_with('#') {
            // Parse marker count from comment: "#Number of markers : 12345"
            if let Some(count_str) = line.split(" : ").nth(1) {
                n_markers = count_str.trim().parse().unwrap_or(0);
            }
            line.clear();
            reader.read_line(&mut line)?;
        }

        let line = line.trim_end_matches('\n').trim_end_matches('\r');
        let columns: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
        let n_individuals = if columns.len() >= 2 {
            (columns.len() - 2) as u16
        } else {
            0
        };

        Ok(TableHeader {
            n_markers,
            n_individuals,
            columns,
        })
    }
}

/// Fast integer parsing for non-negative integers (matching C++ `fast_stoi`).
/// Saturates at u16::MAX instead of wrapping, to avoid silent corruption on
/// high-depth tags (>65535 reads of one RAD marker in one individual).
#[inline(always)]
pub fn fast_parse_u16(bytes: &[u8]) -> u16 {
    let mut val: u16 = 0;
    for &b in bytes {
        if val > u16::MAX / 10 {
            return u16::MAX;
        }
        val = val.saturating_mul(10).saturating_add((b - b'0') as u16);
    }
    val
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_fast_parse_u16() {
        assert_eq!(fast_parse_u16(b"0"), 0);
        assert_eq!(fast_parse_u16(b"42"), 42);
        assert_eq!(fast_parse_u16(b"65535"), 65535);
        // Saturates instead of wrapping (prevents silent depth corruption)
        assert_eq!(fast_parse_u16(b"65536"), u16::MAX);
        assert_eq!(fast_parse_u16(b"999999"), u16::MAX);
    }

    #[test]
    fn test_table_header() {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, "#Number of markers : 100").unwrap();
        writeln!(f, "id\tsequence\tind1\tind2\tind3").unwrap();
        let header = TableHeader::from_file(f.path()).unwrap();
        assert_eq!(header.n_markers, 100);
        assert_eq!(header.n_individuals, 3);
        assert_eq!(header.columns.len(), 5);
    }
}
