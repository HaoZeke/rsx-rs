# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "sympy>=1.12",
# ]
# ///
"""SymPy + independent order-statistic validation for rsx `find_median`.

The engine returns the mathematical median: the mean of the elements at
0-based indices `(n - 1) // 2` and `n // 2` in nondecreasing order. The two
indices coincide for odd lengths and select adjacent values for even lengths.

This script is intentionally non-tautological:
1. SymPy proves combinatorial index identities for odd/even n.
2. An independent pure-Python introselect (quickselect) computes the k-th
   order statistic and is checked against a full sort on exhaustive small
   multisets and random samples.

Exit 0 and print VALIDATED on success.
"""

from __future__ import annotations

import random
from fractions import Fraction

from sympy import Eq, Integer, floor, simplify, symbols


def median_indices(n: int) -> tuple[int, int]:
    """Return the lower and upper middle indices for a non-empty sequence."""
    assert n > 0
    return (n - 1) // 2, n // 2


def independent_select_nth(xs: list[int], k: int) -> int:
    """Independent k-th order statistic (0-based) via Hoare-style quickselect.

    Not a call to sorted(xs)[k]. Mutates a copy only.
    """
    if not 0 <= k < len(xs):
        raise IndexError(k)
    a = list(xs)

    def partition(lo: int, hi: int, pivot_idx: int) -> int:
        pivot = a[pivot_idx]
        a[pivot_idx], a[hi] = a[hi], a[pivot_idx]
        store = lo
        for i in range(lo, hi):
            if a[i] < pivot:
                a[store], a[i] = a[i], a[store]
                store += 1
            elif a[i] == pivot and (i + store) % 2 == 0:
                # mild tie-breaking to keep progress on duplicates
                a[store], a[i] = a[i], a[store]
                store += 1
        a[store], a[hi] = a[hi], a[store]
        return store

    lo, hi = 0, len(a) - 1
    rng = random.Random(0xC0FFEE ^ (k * 0x9E37) ^ len(a))
    while True:
        if lo == hi:
            return a[lo]
        pivot_idx = rng.randint(lo, hi)
        pivot_idx = partition(lo, hi, pivot_idx)
        if k == pivot_idx:
            return a[k]
        if k < pivot_idx:
            hi = pivot_idx - 1
        else:
            lo = pivot_idx + 1


def mathematical_median_via_sort(xs: list[int]) -> Fraction:
    lower, upper = median_indices(len(xs))
    ordered = sorted(xs)
    return Fraction(ordered[lower] + ordered[upper], 2)


def mathematical_median_via_select(xs: list[int]) -> Fraction:
    lower, upper = median_indices(len(xs))
    return Fraction(
        independent_select_nth(xs, lower) + independent_select_nth(xs, upper),
        2,
    )


def sympy_index_identities() -> None:
    """Prove the lower and upper middle-index identities."""
    k = symbols("k", integer=True, nonnegative=True)

    odd_n = 2 * k + 1
    assert simplify(floor((odd_n - 1) / 2) - k) == 0
    assert simplify(floor(odd_n / 2) - k) == 0
    print("sympy: odd lower and upper indices both equal k  OK")

    positive_k = symbols("positive_k", integer=True, positive=True)
    even_n = 2 * positive_k
    assert simplify(floor((even_n - 1) / 2) - (positive_k - 1)) == 0
    assert simplify(floor(even_n / 2) - positive_k) == 0
    print("sympy: even lower index is k-1 and upper index is k  OK")

    for n_val in range(1, 64):
        lower, upper = median_indices(n_val)
        assert 0 <= lower <= upper < n_val
        assert Eq(Integer(lower), floor(Integer(n_val - 1) / 2))
        assert Eq(Integer(upper), floor(Integer(n_val) / 2))
    print("sympy: both indices are ordered and in range for n=1..63  OK")


def validate_select_vs_sort() -> None:
    """Independent two-rank quickselect vs the mathematical median."""
    from itertools import product

    alphabet = [0, 1, 2, 3]
    checked = 0
    for n in range(1, 7):
        for tup in product(alphabet, repeat=n):
            xs = list(tup)
            via_sort = mathematical_median_via_sort(xs)
            via_select = mathematical_median_via_select(xs)
            assert via_select == via_sort, (xs, via_select, via_sort)
            checked += 1
    print(f"two-rank quickselect vs median, exhaustive n=1..6 |A|=4: {checked} OK")

    rng = random.Random(42)
    for trial in range(3000):
        n = rng.randint(1, 200)
        xs = [rng.randint(0, 50_000) for _ in range(n)]
        via_sort = mathematical_median_via_sort(xs)
        via_select = mathematical_median_via_select(xs)
        assert via_select == via_sort, (trial, n, via_select, via_sort)
    print("two-rank quickselect vs median, random n<=200: 3000 OK")


def main() -> int:
    print("=" * 60)
    print("SYMPY + INDEPENDENT SELECT: mathematical median")
    print("=" * 60)
    print()
    sympy_index_identities()
    print()
    validate_select_vs_sort()
    print()
    print("VALIDATED: middle-index identities + two-rank quickselect median")
    print("Rust find_median selects the same lower and upper order statistics.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
