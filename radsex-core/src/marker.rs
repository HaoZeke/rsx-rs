// GPL-3.0-or-later
// Copyright 2024--present radsex-rs developers

//! Marker: a DNA sequence with per-individual depth counts.

use std::io::{self, Write};

/// A single RAD-seq marker with its depth across all individuals.
#[derive(Debug, Clone)]
pub struct Marker {
    /// Marker ID (row number in the depth table).
    pub id: String,
    /// DNA sequence.
    pub sequence: String,
    /// Depth of this marker in each individual (ordered by table columns).
    pub individual_depths: Vec<u16>,
    /// Count of individuals per group where marker depth >= min_depth.
    /// Uses AHashMap for fast hashing (the hot path in table parsing).
    pub group_counts: ahash::AHashMap<String, u32>,
    /// Total number of individuals where marker is present (depth >= min_depth).
    pub n_individuals: u32,
    /// P-value of association with group.
    pub p: f64,
    /// Bonferroni-corrected p-value.
    pub p_corrected: f64,
}

impl Marker {
    /// Create a new marker with space for `n_individuals` depth slots.
    pub fn new(n_individuals: u16) -> Self {
        Marker {
            id: String::new(),
            sequence: String::new(),
            individual_depths: vec![0; n_individuals as usize],
            group_counts: ahash::AHashMap::new(),
            n_individuals: 0,
            p: 0.0,
            p_corrected: 0.0,
        }
    }

    /// Reset marker fields for reuse (avoids reallocation).
    pub fn reset(&mut self, keep_sequence: bool) {
        if !keep_sequence {
            self.id.clear();
            self.sequence.clear();
        }
        for d in &mut self.individual_depths {
            *d = 0;
        }
        for v in self.group_counts.values_mut() {
            *v = 0;
        }
        self.n_individuals = 0;
        self.p = 0.0;
        self.p_corrected = 0.0;
    }

    /// Write this marker in TSV table format:
    /// id\tsequence\tdepth1\t...\tdepthN\n
    pub fn write_as_table<W: Write>(&self, w: &mut W) -> io::Result<()> {
        write!(w, "{}\t{}", self.id, self.sequence)?;
        for d in &self.individual_depths {
            write!(w, "\t{d}")?;
        }
        writeln!(w)
    }

    /// Write this marker in FASTA format:
    /// `>id_group1:count_group2:count_p:pval_pcorr:pcorr_mindepth:md`
    /// followed by the sequence on the next line.
    pub fn write_as_fasta<W: Write>(&self, w: &mut W, min_depth: u32) -> io::Result<()> {
        write!(w, ">{}", self.id)?;
        // Sort group names for deterministic output
        let mut groups: Vec<(&String, &u32)> = self.group_counts.iter().collect();
        groups.sort_by_key(|(k, _)| k.as_str());
        for (group, count) in groups {
            write!(w, "_{group}:{count}")?;
        }
        writeln!(w, "_p:{}_pcorr:{}_mindepth:{min_depth}", self.p, self.p_corrected)?;
        writeln!(w, "{}", self.sequence)
    }
}

/// A marker that has been aligned to a reference genome.
#[derive(Debug, Clone)]
pub struct AlignedMarker {
    pub id: String,
    pub contig: String,
    pub position: i64,
    pub bias: f64,
    pub p: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_marker_table_output() {
        let mut m = Marker::new(3);
        m.id = "42".to_string();
        m.sequence = "ATCG".to_string();
        m.individual_depths = vec![10, 0, 5];
        let mut buf = Vec::new();
        m.write_as_table(&mut buf).unwrap();
        assert_eq!(String::from_utf8(buf).unwrap(), "42\tATCG\t10\t0\t5\n");
    }

    #[test]
    fn test_marker_fasta_output() {
        let mut m = Marker::new(2);
        m.id = "1".to_string();
        m.sequence = "GATTACA".to_string();
        m.group_counts.insert("M".to_string(), 5);
        m.group_counts.insert("F".to_string(), 1);
        m.p = 0.001;
        m.p_corrected = 0.01;
        let mut buf = Vec::new();
        m.write_as_fasta(&mut buf, 5).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.starts_with(">1_F:1_M:5_p:0.001_pcorr:0.01_mindepth:5\n"));
        assert!(output.contains("GATTACA\n"));
    }

    #[test]
    fn test_marker_reset() {
        let mut m = Marker::new(2);
        m.id = "42".to_string();
        m.individual_depths = vec![10, 5];
        m.n_individuals = 2;
        m.reset(false);
        assert!(m.id.is_empty());
        assert_eq!(m.individual_depths, vec![0, 0]);
        assert_eq!(m.n_individuals, 0);
    }
}
