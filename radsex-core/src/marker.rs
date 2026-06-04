// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Marker: a DNA sequence with per-individual depth counts.

use crate::bitset::{BitsetRow, GroupMask};
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
    /// Bitset: bit `i` set iff individual `i` has depth >= min_depth.
    /// Group counting via `presence.count_masked(&group_mask)`.
    pub presence: BitsetRow,
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
            presence: BitsetRow::new(n_individuals),
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
        self.presence.clear();
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

    /// Write this marker in FASTA format with group counts computed from bitset.
    /// `>id_group1:count_group2:count_p:pval_pcorr:pcorr_mindepth:md`
    /// followed by the sequence on the next line.
    pub fn write_as_fasta_bitset<W: Write>(
        &self,
        w: &mut W,
        min_depth: u32,
        group_names: &[(String, &GroupMask)],
    ) -> io::Result<()> {
        write!(w, ">{}", self.id)?;
        // Sort by group name for deterministic output
        let mut groups: Vec<(&str, u32)> = group_names
            .iter()
            .map(|(name, mask)| (name.as_str(), self.presence.count_masked(mask)))
            .collect();
        groups.sort_by_key(|(k, _)| *k);
        for (group, count) in groups {
            write!(w, "_{group}:{count}")?;
        }
        writeln!(
            w,
            "_p:{}_pcorr:{}_mindepth:{min_depth}",
            self.p, self.p_corrected
        )?;
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
    fn test_marker_fasta_bitset_output() {
        let mut m = Marker::new(4);
        m.id = "1".to_string();
        m.sequence = "GATTACA".to_string();
        // Set presence: individuals 0,1,2 present (3 males), individual 3 present (1 female)
        m.presence.set(0);
        m.presence.set(1);
        m.presence.set(2);
        m.presence.set(3);
        m.p = 0.001;
        m.p_corrected = 0.01;

        let mask_m = GroupMask::from_columns(
            &[
                "".into(),
                "".into(),
                "M".into(),
                "M".into(),
                "M".into(),
                "F".into(),
            ],
            "M",
            4,
        );
        let mask_f = GroupMask::from_columns(
            &[
                "".into(),
                "".into(),
                "M".into(),
                "M".into(),
                "M".into(),
                "F".into(),
            ],
            "F",
            4,
        );

        let group_names = vec![("M".to_string(), &mask_m), ("F".to_string(), &mask_f)];

        let mut buf = Vec::new();
        m.write_as_fasta_bitset(&mut buf, 5, &group_names).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains(">1_F:1_M:3_p:0.001_pcorr:0.01_mindepth:5\n"));
        assert!(output.contains("GATTACA\n"));
    }

    #[test]
    fn test_marker_reset() {
        let mut m = Marker::new(2);
        m.id = "42".to_string();
        m.individual_depths = vec![10, 5];
        m.presence.set(0);
        m.presence.set(1);
        m.n_individuals = 2;
        m.reset(false);
        assert!(m.id.is_empty());
        assert_eq!(m.individual_depths, vec![0, 0]);
        assert_eq!(m.n_individuals, 0);
        assert_eq!(m.presence.count_total(), 0);
    }
}
