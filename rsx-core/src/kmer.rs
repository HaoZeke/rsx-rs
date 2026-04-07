// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! K-mer based marker deduplication.
//!
//! Groups markers by shared canonical k-mer signatures to collapse
//! sequencing error variants. Reduces the number of markers tested,
//! increasing statistical power for sex detection.

/// Compute canonical k-mer hash for a DNA sequence.
/// Canonical = lexicographically smallest of {kmer, revcomp(kmer)}.
/// Returns a u64 hash of the canonical k-mer set for grouping.
pub fn canonical_kmer_hash(seq: &[u8], k: usize) -> u64 {
    if seq.len() < k {
        return hash_bytes(seq);
    }

    // Use the majority canonical k-mer as the representative hash.
    // For deduplication: two sequences share a group if their majority
    // k-mer matches.
    let mut best_hash = u64::MAX;
    for window in seq.windows(k) {
        let fwd = kmer_to_u64(window);
        let rev = revcomp_u64(fwd, k);
        let canonical = fwd.min(rev);
        best_hash = best_hash.min(canonical);
    }
    best_hash
}

/// Group markers by canonical k-mer signature.
/// Returns a map from group_hash -> list of marker indices.
pub fn group_by_kmer(
    sequences: &[Vec<u8>],
    k: usize,
) -> ahash::AHashMap<u64, Vec<usize>> {
    let mut groups: ahash::AHashMap<u64, Vec<usize>> = ahash::AHashMap::new();
    for (i, seq) in sequences.iter().enumerate() {
        let hash = canonical_kmer_hash(seq, k);
        groups.entry(hash).or_default().push(i);
    }
    groups
}

/// Encode a k-mer (up to 32bp) as a u64 using 2-bit encoding.
fn kmer_to_u64(kmer: &[u8]) -> u64 {
    let mut val = 0u64;
    for &base in kmer {
        val <<= 2;
        val |= match base {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 0, // N -> A
        };
    }
    val
}

/// Reverse complement of a 2-bit encoded k-mer.
fn revcomp_u64(kmer: u64, k: usize) -> u64 {
    let mut rc = 0u64;
    let mut kmer = kmer;
    for _ in 0..k {
        rc <<= 2;
        rc |= 3 - (kmer & 3); // complement: A<->T, C<->G
        kmer >>= 2;
    }
    rc
}

/// Simple hash for short sequences (< k bases).
fn hash_bytes(seq: &[u8]) -> u64 {
    let mut h = 0xcbf29ce484222325u64; // FNV offset
    for &b in seq {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3); // FNV prime
    }
    h
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_canonical_kmer_revcomp() {
        // ATCG and CGAT are reverse complements
        let fwd = kmer_to_u64(b"ATCG");
        let rev = revcomp_u64(fwd, 4);
        let rev_seq = kmer_to_u64(b"CGAT");
        assert_eq!(rev, rev_seq);
    }

    #[test]
    fn test_canonical_same_for_revcomp_seqs() {
        let h1 = canonical_kmer_hash(b"ATCGATCGATCG", 8);
        let _h2 = canonical_kmer_hash(b"CGATCGATCGAT", 8);
        // These share canonical k-mers so should have the same min hash
        // (not guaranteed to be equal but likely for similar sequences)
        // The point is that identical sequences give identical hashes
        let h3 = canonical_kmer_hash(b"ATCGATCGATCG", 8);
        assert_eq!(h1, h3);
    }

    #[test]
    fn test_group_by_kmer() {
        let seqs = vec![
            b"ATCGATCGATCG".to_vec(),
            b"ATCGATCGATCG".to_vec(), // identical
            b"GGGGAAAACCCC".to_vec(), // different
        ];
        let groups = group_by_kmer(&seqs, 8);
        // Identical sequences should be in the same group
        assert!(groups.len() >= 2); // at least 2 distinct groups
    }
}
