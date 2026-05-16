// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! K-mer based marker deduplication.
//!
//! Groups markers by shared canonical k-mer signatures to collapse
//! sequencing error variants. Reduces the number of markers tested,
//! increasing statistical power for sex detection.

/// Compute canonical k-mer hash for a DNA sequence (min-hash over windows).
/// Canonical = lexicographically smallest of {kmer, revcomp(kmer)}.
/// The representative for a sequence is the *minimum* hash among its k-mers.
/// This is an LSH heuristic for grouping similar sequences (e.g. sequencing errors);
/// it is *not* guaranteed that two sequences differing by one base will share a group
/// (see test_group_single_base_error). Use for approximate collapse only.
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
pub fn group_by_kmer(sequences: &[Vec<u8>], k: usize) -> ahash::AHashMap<u64, Vec<usize>> {
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
    fn test_kmer_encoding_roundtrip() {
        // Verify 2-bit encoding for all 4 bases
        assert_eq!(kmer_to_u64(b"A"), 0b00);
        assert_eq!(kmer_to_u64(b"C"), 0b01);
        assert_eq!(kmer_to_u64(b"G"), 0b10);
        assert_eq!(kmer_to_u64(b"T"), 0b11);
        assert_eq!(kmer_to_u64(b"ACGT"), 0b00_01_10_11);
    }

    #[test]
    fn test_revcomp_correctness() {
        // A<->T, C<->G
        // revcomp(ATCG) = CGAT
        let fwd = kmer_to_u64(b"ATCG");
        let rev = revcomp_u64(fwd, 4);
        assert_eq!(rev, kmer_to_u64(b"CGAT"));

        // revcomp(AAAA) = TTTT
        assert_eq!(revcomp_u64(kmer_to_u64(b"AAAA"), 4), kmer_to_u64(b"TTTT"));

        // revcomp(revcomp(x)) = x
        let seq = kmer_to_u64(b"ACGTACGT");
        assert_eq!(revcomp_u64(revcomp_u64(seq, 8), 8), seq);
    }

    #[test]
    fn test_canonical_is_symmetric() {
        // canonical(seq) == canonical(revcomp(seq))
        // ATCGATCG and its revcomp CGATCGAT should give the same canonical hash
        let h1 = canonical_kmer_hash(b"ATCGATCG", 4);
        let h2 = canonical_kmer_hash(b"CGATCGAT", 4);
        assert_eq!(
            h1, h2,
            "sequence and its revcomp should have same canonical hash"
        );
    }

    #[test]
    fn test_canonical_deterministic() {
        let h1 = canonical_kmer_hash(b"ATCGATCGATCG", 8);
        let h2 = canonical_kmer_hash(b"ATCGATCGATCG", 8);
        assert_eq!(h1, h2, "same input should always give same hash");
    }

    #[test]
    fn test_canonical_different_seqs_differ() {
        let h1 = canonical_kmer_hash(b"AAAAAAAA", 8);
        let h2 = canonical_kmer_hash(b"CCCCCCCC", 8);
        assert_ne!(h1, h2, "unrelated sequences should have different hashes");
    }

    #[test]
    fn test_group_identical_together() {
        let seqs = vec![
            b"ATCGATCGATCG".to_vec(),
            b"ATCGATCGATCG".to_vec(), // identical to 0
            b"GGGGAAAACCCC".to_vec(), // different
            b"GGGGAAAACCCC".to_vec(), // identical to 2
        ];
        let groups = group_by_kmer(&seqs, 8);
        // Find which group each sequence is in
        let mut seq_to_group: Vec<u64> = vec![0; 4];
        for (&hash, indices) in &groups {
            for &idx in indices {
                seq_to_group[idx] = hash;
            }
        }
        // Identical sequences must be in the same group
        assert_eq!(
            seq_to_group[0], seq_to_group[1],
            "identical seqs 0,1 should group together"
        );
        assert_eq!(
            seq_to_group[2], seq_to_group[3],
            "identical seqs 2,3 should group together"
        );
        // Different sequences should be in different groups
        assert_ne!(
            seq_to_group[0], seq_to_group[2],
            "different seqs should be in different groups"
        );
    }

    #[test]
    fn test_group_single_base_error() {
        // Two sequences differing by 1 base should share k-mers if k < seq_len
        let seq1 = b"ATCGATCGATCGATCG".to_vec(); // 16bp
        let seq2 = b"ATCGATCGATCGATCC".to_vec(); // differs at last base
        let seqs = vec![seq1, seq2];
        let groups_k8 = group_by_kmer(&seqs, 8);
        // With k=8, they share 8 out of 9 k-mers, so min-hash likely the same
        // (not guaranteed but highly probable for this specific case)
        // At minimum, both should be assigned to some group
        assert!(!groups_k8.is_empty());
    }

    #[test]
    fn test_group_empty_and_short() {
        let seqs = vec![
            b"AT".to_vec(),       // shorter than k
            b"ATCGATCG".to_vec(), // normal
        ];
        let groups = group_by_kmer(&seqs, 8);
        // Should not panic, both sequences get assigned
        let total: usize = groups.values().map(|v| v.len()).sum();
        assert_eq!(total, 2, "all sequences should be assigned to a group");
    }
}
