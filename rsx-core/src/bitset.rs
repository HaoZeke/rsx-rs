// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Bitset representation for marker presence/absence.
//!
//! Stores `depth >= min_depth` as a single bit per (marker, individual),
//! enabling group counts via `popcount(row & group_mask)` instead of
//! HashMap lookups. This gives 10-16x memory reduction and eliminates
//! all hashing overhead in the hot path.

/// A fixed-width bitset row representing one marker across N individuals.
/// Internally stored as a `Vec<u64>` where bit `i` = individual `i` present.
#[derive(Debug, Clone)]
pub struct BitsetRow {
    words: Vec<u64>,
}

impl BitsetRow {
    /// Create a new zeroed bitset for `n_individuals`.
    #[inline]
    pub fn new(n_individuals: u16) -> Self {
        let n_words = (n_individuals as usize).div_ceil(64);
        BitsetRow {
            words: vec![0u64; n_words],
        }
    }

    /// Set bit `i` (marking individual `i` as present).
    #[inline(always)]
    pub fn set(&mut self, i: usize) {
        let word = i / 64;
        let bit = i % 64;
        if word < self.words.len() {
            self.words[word] |= 1u64 << bit;
        }
    }

    /// Count bits set in `self & mask` (group count via popcount).
    #[inline]
    pub fn count_masked(&self, mask: &GroupMask) -> u32 {
        debug_assert_eq!(self.words.len(), mask.words.len());
        let mut count = 0u32;
        for (w, m) in self.words.iter().zip(mask.words.iter()) {
            count += (w & m).count_ones();
        }
        count
    }

    /// Total number of set bits (n_individuals present).
    #[inline]
    pub fn count_total(&self) -> u32 {
        let mut count = 0u32;
        for w in &self.words {
            count += w.count_ones();
        }
        count
    }

    /// Clear all bits to zero (for reuse).
    #[inline]
    pub fn clear(&mut self) {
        for w in &mut self.words {
            *w = 0;
        }
    }
}

/// Pre-computed bitmask for a group (e.g. all males or all females).
/// Bit `i` is set if individual `i` belongs to this group.
#[derive(Debug, Clone)]
pub struct GroupMask {
    words: Vec<u64>,
}

impl GroupMask {
    /// Build a group mask from per-column group labels.
    /// `column_groups[i]` is the group name for individual at column index `i`.
    /// Only columns matching `group_name` get their bit set.
    pub fn from_columns(column_groups: &[String], group_name: &str, n_individuals: u16) -> Self {
        let n_words = (n_individuals as usize).div_ceil(64);
        let mut words = vec![0u64; n_words];

        // column_groups has 2 extra entries (id, sequence) at indices 0,1
        for (col_idx, group) in column_groups.iter().enumerate().skip(2) {
            if group == group_name {
                let ind_idx = col_idx - 2;
                let word = ind_idx / 64;
                let bit = ind_idx % 64;
                if word < n_words {
                    words[word] |= 1u64 << bit;
                }
            }
        }

        GroupMask { words }
    }

    /// Number of set bits (total individuals in this group).
    pub fn count(&self) -> u32 {
        self.words.iter().map(|w| w.count_ones()).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bitset_set_and_count() {
        let mut row = BitsetRow::new(10);
        row.set(0);
        row.set(3);
        row.set(9);
        assert_eq!(row.count_total(), 3);
    }

    #[test]
    fn test_bitset_masked_count() {
        let mut row = BitsetRow::new(10);
        row.set(0);
        row.set(1);
        row.set(2);
        row.set(5);
        row.set(6);

        // Mask: individuals 0-4 are group1
        let mut mask = GroupMask { words: vec![0u64; 1] };
        mask.words[0] = 0b11111; // bits 0-4
        assert_eq!(row.count_masked(&mask), 3); // 0,1,2 match

        // Mask: individuals 5-9 are group2
        let mut mask2 = GroupMask { words: vec![0u64; 1] };
        mask2.words[0] = 0b1111100000; // bits 5-9
        assert_eq!(row.count_masked(&mask2), 2); // 5,6 match
    }

    #[test]
    fn test_bitset_clear() {
        let mut row = BitsetRow::new(10);
        row.set(0);
        row.set(5);
        assert_eq!(row.count_total(), 2);
        row.clear();
        assert_eq!(row.count_total(), 0);
    }

    #[test]
    fn test_bitset_large() {
        // 200 individuals = 4 u64 words
        let mut row = BitsetRow::new(200);
        for i in 0..200 {
            row.set(i);
        }
        assert_eq!(row.count_total(), 200);

        // Mask first 100
        let n_words = (200 + 63) / 64;
        let mut mask = GroupMask { words: vec![0u64; n_words] };
        for i in 0..100 {
            mask.words[i / 64] |= 1u64 << (i % 64);
        }
        assert_eq!(row.count_masked(&mask), 100);
    }

    #[test]
    fn test_group_mask_from_columns() {
        let _columns: Vec<String> = vec![
            "id".into(), "sequence".into(),
            "ind1".into(), "ind2".into(), "ind3".into(), "ind4".into(),
        ];
        let groups: Vec<String> = vec![
            "".into(), "".into(),
            "M".into(), "M".into(), "F".into(), "F".into(),
        ];

        let mask_m = GroupMask::from_columns(&groups, "M", 4);
        assert_eq!(mask_m.count(), 2); // ind1, ind2

        let mask_f = GroupMask::from_columns(&groups, "F", 4);
        assert_eq!(mask_f.count(), 2); // ind3, ind4
    }
}
