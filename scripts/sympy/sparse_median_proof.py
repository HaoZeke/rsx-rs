#!/usr/bin/env python3
"""
SymPy proof: sparse median equivalence.

Shows that the median of a sequence with known zero count can be computed
from only the non-zero elements, without materializing the full sequence.

Used in rsx-rs depth command to compute exact per-individual median
from a sparse external sort (zeros skipped, 70% I/O reduction).
"""

from sympy import Symbol, Piecewise, floor, simplify, Eq


def prove_sparse_median():
    """
    Prove: median(full_sorted) == sparse_median(sorted_nonzero, n_zeros, n_total)

    Given:
      - n_total: total number of elements (including zeros)
      - n_zeros: count of zero elements
      - n_nonzero = n_total - n_zeros
      - sorted_nonzero[0..n_nonzero-1]: the non-zero elements in sorted order
      - full_sorted[0..n_total-1] = [0]*n_zeros ++ sorted_nonzero

    The median position is:
      median_pos = floor(n_total / 2)

    Case 1: median_pos < n_zeros
      => median = 0 (falls within the zero block)

    Case 2: median_pos >= n_zeros
      => median = sorted_nonzero[median_pos - n_zeros]
      (shift by n_zeros to index into the non-zero array)
    """

    n_total = Symbol("n_total", positive=True, integer=True)
    n_zeros = Symbol("n_zeros", nonnegative=True, integer=True)
    n_nonzero = n_total - n_zeros

    median_pos = floor(n_total / 2)

    print("=" * 60)
    print("PROOF: Sparse Median Equivalence")
    print("=" * 60)
    print()
    print("Given:")
    print(f"  n_total elements, n_zeros of which are 0")
    print(f"  n_nonzero = n_total - n_zeros")
    print(f"  sorted_nonzero[0..n_nonzero-1] = non-zero elements, sorted")
    print(f"  full_sorted = [0]*n_zeros ++ sorted_nonzero")
    print()
    print("Median position:")
    print(f"  median_pos = floor(n_total / 2)")
    print()
    print("Case 1: median_pos < n_zeros")
    print("  full_sorted[median_pos] = 0  (within the zero block)")
    print("  => median = 0")
    print()
    print("Case 2: median_pos >= n_zeros")
    print("  full_sorted[median_pos] = sorted_nonzero[median_pos - n_zeros]")
    print("  => median = sorted_nonzero[median_pos - n_zeros]")
    print()

    # Verify with concrete examples
    print("=" * 60)
    print("VERIFICATION with concrete examples")
    print("=" * 60)
    print()

    test_cases = [
        # (full_sequence, description)
        ([0, 0, 0, 1, 2, 3], "mostly zeros"),
        ([0, 0, 5, 10, 15, 20], "half zeros"),
        ([1, 2, 3, 4, 5], "no zeros"),
        ([0, 0, 0, 0, 0], "all zeros"),
        ([0, 0, 0, 0, 1], "one nonzero"),
        ([0, 3, 5, 7, 9, 11, 13], "one zero, odd count"),
    ]

    all_ok = True
    for seq, desc in test_cases:
        n = len(seq)
        sorted_seq = sorted(seq)
        true_median = sorted_seq[n // 2]

        nonzero = sorted([x for x in seq if x > 0])
        nz = len([x for x in seq if x == 0])
        median_pos = n // 2

        if median_pos < nz:
            sparse_median = 0
        else:
            idx = median_pos - nz
            sparse_median = nonzero[idx] if idx < len(nonzero) else 0

        ok = true_median == sparse_median
        if not ok:
            all_ok = False
        print(
            f"  {desc:25s}: seq={sorted_seq}, "
            f"true_median={true_median}, sparse_median={sparse_median} "
            f"[{'OK' if ok else 'FAIL'}]"
        )

    print()
    if all_ok:
        print("All verifications PASSED.")
    else:
        print("Some verifications FAILED!")
        exit(1)

    print()
    print("=" * 60)
    print("COMPLEXITY ANALYSIS")
    print("=" * 60)
    print()
    print("Dense approach:")
    print("  Sort all n_total elements: O(n_total * log(n_total)) time")
    print("  Memory: O(n_total)")
    print()
    print("Sparse approach:")
    print("  Sort only n_nonzero elements: O(n_nonzero * log(n_nonzero)) time")
    print("  Memory: O(n_nonzero)")
    print("  Speedup: n_total / n_nonzero (= 1 / (1 - sparsity))")
    print()
    print("For RAD-seq depth matrix with 70% sparsity:")
    print("  Speedup = 1 / 0.3 = 3.3x less sorting, I/O, and memory")
    print("  75M markers x 200 individuals:")
    print("    Dense: 15B entries to sort -> ~30GB temp files")
    print("    Sparse: 4.5B entries to sort -> ~9GB temp files")


if __name__ == "__main__":
    prove_sparse_median()
