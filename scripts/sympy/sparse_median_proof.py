#!/usr/bin/env python3
"""
SymPy proof: sparse median equivalence.

Shows that the median of a sequence with known zero count can be computed
from only the non-zero elements, without materializing the full sequence.

Used in rsx-rs depth command to compute exact per-individual median
from a sparse external sort (zeros skipped, 70% I/O reduction).
"""

from fractions import Fraction
from itertools import product

from sympy import floor, simplify, symbols


def sparse_order_statistic(sorted_nonzero: list[int], n_zeros: int, rank: int) -> int:
    """Read one rank from `[0] * n_zeros + sorted_nonzero`."""
    if rank < n_zeros:
        return 0
    return sorted_nonzero[rank - n_zeros]


def sparse_median(sorted_nonzero: list[int], n_zeros: int) -> Fraction:
    """Compute the mathematical median without materializing the zero prefix."""
    n_total = n_zeros + len(sorted_nonzero)
    assert n_total > 0
    lower_rank = (n_total - 1) // 2
    upper_rank = n_total // 2
    lower = sparse_order_statistic(sorted_nonzero, n_zeros, lower_rank)
    upper = sparse_order_statistic(sorted_nonzero, n_zeros, upper_rank)
    return Fraction(lower + upper, 2)


def prove_sparse_median():
    """
    Prove: median(full_sorted) == sparse_median(sorted_nonzero, n_zeros, n_total)

    Given:
      - n_total: total number of elements (including zeros)
      - n_zeros: count of zero elements
      - n_nonzero = n_total - n_zeros
      - sorted_nonzero[0..n_nonzero-1]: the non-zero elements in sorted order
      - full_sorted[0..n_total-1] = [0]*n_zeros ++ sorted_nonzero

    The middle positions are:
      lower_rank = floor((n_total - 1) / 2)
      upper_rank = floor(n_total / 2)

    Each rank maps independently through the implicit zero prefix:
      value(rank) = 0, if rank < n_zeros
                  = sorted_nonzero[rank - n_zeros], otherwise

    The mathematical median is (value(lower_rank) + value(upper_rank)) / 2.
    """

    k = symbols("k", integer=True, nonnegative=True)
    assert simplify(floor(((2 * k + 1) - 1) / 2) - k) == 0
    assert simplify(floor((2 * k + 1) / 2) - k) == 0
    positive_k = symbols("positive_k", integer=True, positive=True)
    assert simplify(floor((2 * positive_k - 1) / 2) - (positive_k - 1)) == 0
    assert simplify(floor((2 * positive_k) / 2) - positive_k) == 0

    print("=" * 60)
    print("PROOF: Sparse Median Equivalence")
    print("=" * 60)
    print()
    print("Given:")
    print("  n_total elements, n_zeros of which are 0")
    print("  n_nonzero = n_total - n_zeros")
    print("  sorted_nonzero[0..n_nonzero-1] = non-zero elements, sorted")
    print("  full_sorted = [0]*n_zeros ++ sorted_nonzero")
    print()
    print("Middle positions:")
    print("  lower_rank = floor((n_total - 1) / 2)")
    print("  upper_rank = floor(n_total / 2)")
    print()
    print("Map each rank through the zero prefix and average both values.")
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
        ([0, 0, 1, 2], "half-integer median"),
        ([1, 2, 3, 4, 5], "no zeros"),
        ([0, 0, 0, 0, 0], "all zeros"),
        ([0, 0, 0, 0, 1], "one nonzero"),
        ([0, 3, 5, 7, 9, 11, 13], "one zero, odd count"),
    ]

    all_ok = True
    for seq, desc in test_cases:
        n = len(seq)
        sorted_seq = sorted(seq)
        lower_rank = (n - 1) // 2
        upper_rank = n // 2
        true_median = Fraction(
            sorted_seq[lower_rank] + sorted_seq[upper_rank], 2
        )

        nonzero = sorted([x for x in seq if x > 0])
        nz = len([x for x in seq if x == 0])
        sparse_result = sparse_median(nonzero, nz)

        ok = true_median == sparse_result
        if not ok:
            all_ok = False
        print(
            f"  {desc:25s}: seq={sorted_seq}, "
            f"true_median={true_median}, sparse_median={sparse_result} "
            f"[{'OK' if ok else 'FAIL'}]"
        )

    print()
    if all_ok:
        print("All verifications PASSED.")
    else:
        print("Some verifications FAILED!")
        raise SystemExit(1)

    exhaustive = 0
    for n in range(1, 8):
        for seq in product(range(4), repeat=n):
            sorted_seq = sorted(seq)
            lower_rank = (n - 1) // 2
            upper_rank = n // 2
            expected = Fraction(
                sorted_seq[lower_rank] + sorted_seq[upper_rank], 2
            )
            nonzero = sorted(value for value in seq if value > 0)
            n_zeros = sum(value == 0 for value in seq)
            assert sparse_median(nonzero, n_zeros) == expected
            exhaustive += 1
    print(f"Exhaustive n=1..7 |A|=4: {exhaustive} PASSED.")

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
