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
        let n_individuals_usize = columns.len().saturating_sub(2);
        let n_individuals = u16::try_from(n_individuals_usize).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "markers table has {n_individuals_usize} individuals; the limit is {}",
                    u16::MAX
                ),
            )
        })?;

        Ok(TableHeader {
            n_markers,
            n_individuals,
            columns,
        })
    }
}

/// Fast integer parsing for non-negative integers (matching C++ `fast_stoi`
/// on digit strings). Saturates at u16::MAX instead of wrapping, so a depth
/// above 65535 cannot silently alias to a small count.
///
/// A field that contains any non-digit is not a depth: it parses as 0 rather
/// than as a wrapped code-point difference. That keeps `"null"`, `"10a"`, and
/// `"-1"` from panicking in the 3–4 byte unrolled path and from being counted
/// as present.
///
/// Depth fields in RAD tables are almost always 1–3 ASCII digits; a short
/// unrolled path avoids the saturating mul loop for that common case.
#[inline(always)]
pub fn fast_parse_u16(bytes: &[u8]) -> u16 {
    match bytes.len() {
        0 => 0,
        1 => digit_u16(bytes[0]).unwrap_or(0),
        2 => match (digit_u16(bytes[0]), digit_u16(bytes[1])) {
            (Some(tens), Some(ones)) => tens * 10 + ones,
            _ => 0,
        },
        3 => match (
            digit_u16(bytes[0]),
            digit_u16(bytes[1]),
            digit_u16(bytes[2]),
        ) {
            (Some(hundreds), Some(tens), Some(ones)) => hundreds * 100 + tens * 10 + ones,
            _ => 0,
        },
        _ => parse_u16_saturating(bytes),
    }
}

#[inline(always)]
fn digit_u16(byte: u8) -> Option<u16> {
    if byte.is_ascii_digit() {
        Some(u16::from(byte - b'0'))
    } else {
        None
    }
}

/// Digit strings of 4 or more characters. Non-digits reject the whole field.
#[inline]
fn parse_u16_saturating(bytes: &[u8]) -> u16 {
    let mut value = 0u32;
    for &byte in bytes {
        let Some(digit) = digit_u16(byte) else {
            return 0;
        };
        value = value.saturating_mul(10).saturating_add(u32::from(digit));
        if value >= u32::from(u16::MAX) {
            return u16::MAX;
        }
    }
    u16::try_from(value).unwrap_or(u16::MAX)
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
