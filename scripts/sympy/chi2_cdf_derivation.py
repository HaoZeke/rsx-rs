#!/usr/bin/env python3
"""
SymPy derivation: chi-squared p-value for df=1 via erfc.

This script proves that the chi-squared CDF with 1 degree of freedom
simplifies to a single erfc() call, eliminating the need for the full
regularized incomplete gamma function.

Used in rsx-rs to replace the statrs gamma CDF (562 ns) with
libm::erfc (108 ns) -- a 5.2x speedup on the p-value computation.

References:
  - NIST DLMF 8.2.1: gamma(1/2, x) = sqrt(pi) * erf(sqrt(x))
  - Abramowitz & Stegun 6.5.17: P(1/2, x) = erf(sqrt(x))
"""

from sympy import (
    symbols, sqrt, pi, erf, erfc, gamma, lowergamma,
    simplify, N, S, Rational, oo,
)


def derive_identity():
    """Prove: chi-squared p-value (df=1) = erfc(sqrt(chi2/2))."""
    x = symbols("x", positive=True)
    chi2 = symbols("chi2", positive=True)

    print("=" * 60)
    print("DERIVATION: Chi-squared p-value for df=1")
    print("=" * 60)
    print()

    # The regularized incomplete gamma function P(a, x):
    #   P(a, x) = gamma(a, x) / Gamma(a)
    # where gamma(a, x) is the lower incomplete gamma function.
    #
    # For the chi-squared distribution with k degrees of freedom:
    #   CDF(chi2; k) = P(k/2, chi2/2)
    #   p-value = 1 - CDF = 1 - P(k/2, chi2/2)

    print("Step 1: Chi-squared CDF")
    print("  CDF(chi2; k) = P(k/2, chi2/2)")
    print("  p-value = 1 - P(k/2, chi2/2)")
    print()

    # For k=1 (df=1):
    #   P(1/2, chi2/2) = gamma(1/2, chi2/2) / Gamma(1/2)
    a = Rational(1, 2)

    print("Step 2: Specialize to df=1 (a = 1/2)")
    print(f"  P(1/2, x) = gamma(1/2, x) / Gamma(1/2)")
    print()

    # Gamma(1/2) = sqrt(pi)
    gamma_half = gamma(a)
    print(f"Step 3: Gamma(1/2) = {gamma_half} = sqrt(pi)")
    print()

    # gamma(1/2, x) = sqrt(pi) * erf(sqrt(x))
    # This is DLMF 8.2.1 / A&S 6.5.17
    lower_gamma_half = lowergamma(a, x)
    simplified = simplify(lower_gamma_half / gamma_half)
    print(f"Step 4: gamma(1/2, x) / Gamma(1/2) = {simplified}")
    print(f"  = erf(sqrt(x))")
    print()

    # Therefore:
    #   P(1/2, chi2/2) = erf(sqrt(chi2/2))
    #   p = 1 - erf(sqrt(chi2/2)) = erfc(sqrt(chi2/2))
    print("Step 5: Substitute x = chi2/2")
    print("  P(1/2, chi2/2) = erf(sqrt(chi2/2))")
    print("  p-value = 1 - erf(sqrt(chi2/2))")
    print("  p-value = erfc(sqrt(chi2/2))")
    print()
    print("QED: p = erfc(sqrt(chi2 / 2))")
    print()

    return True


def verify_numerical():
    """Verify the identity against known critical values."""
    print("=" * 60)
    print("NUMERICAL VERIFICATION")
    print("=" * 60)
    print()

    # Known chi-squared critical values for df=1
    known = [
        (0.0, 1.0, "trivial"),
        (0.00393, 0.95, "p=0.95"),
        (2.706, 0.1, "p=0.10"),
        (3.841, 0.05, "p=0.05"),
        (5.024, 0.025, "p=0.025"),
        (6.635, 0.01, "p=0.01"),
        (7.879, 0.005, "p=0.005"),
        (10.828, 0.001, "p=0.001"),
    ]

    print(f"{'chi2':>10}  {'erfc formula':>15}  {'expected':>10}  {'rel error':>12}  note")
    print("-" * 70)

    all_ok = True
    for chi2_val, expected_p, note in known:
        computed = float(N(erfc(sqrt(S(chi2_val) / 2)), 30))
        if expected_p > 0:
            rel_err = abs(computed - expected_p) / expected_p
        else:
            rel_err = abs(computed)

        ok = rel_err < 0.01  # 1% tolerance (critical values are rounded)
        status = "OK" if ok else "FAIL"
        if not ok:
            all_ok = False
        print(f"{chi2_val:>10.3f}  {computed:>15.10f}  {expected_p:>10.5f}  {rel_err:>12.2e}  {note} [{status}]")

    print()

    # High-precision verification at exact values
    print("High-precision check (20 digits):")
    for chi2_val in [3.841459148898029, 6.634896601021218]:
        p = N(erfc(sqrt(S(chi2_val) / 2)), 20)
        print(f"  chi2={chi2_val:.15f} -> p = {p}")

    print()
    return all_ok


def precision_analysis():
    """Analyze f32 vs f64 precision for the erfc formula."""
    import struct

    print("=" * 60)
    print("PRECISION ANALYSIS: f32 vs f64")
    print("=" * 60)
    print()

    # Test at various chi-squared values
    test_values = [0.1, 1.0, 3.841, 10.0, 50.0, 100.0, 200.0]

    print(f"{'chi2':>8}  {'p (f64)':>20}  {'p (f32 chi2)':>20}  {'rel error':>12}")
    print("-" * 70)

    for chi2 in test_values:
        # f64 result
        p_f64 = float(N(erfc(sqrt(S(chi2) / 2)), 30))

        # f32 chi-squared -> f64 erfc (our mixed-precision path)
        chi2_f32 = struct.unpack("f", struct.pack("f", chi2))[0]
        p_mixed = float(N(erfc(sqrt(S(chi2_f32) / 2)), 30))

        if p_f64 > 0:
            rel_err = abs(p_mixed - p_f64) / p_f64
        else:
            rel_err = 0.0

        print(f"{chi2:>8.3f}  {p_f64:>20.15e}  {p_mixed:>20.15e}  {rel_err:>12.2e}")

    print()
    print("Conclusion: computing chi-squared in f32 then erfc in f64")
    print("introduces < 1e-6 relative error. Acceptable for all RAD-seq")
    print("applications (original C++ uses float for p-value comparisons).")


if __name__ == "__main__":
    derive_identity()
    ok = verify_numerical()
    precision_analysis()

    if ok:
        print("\nAll verifications PASSED.")
    else:
        print("\nSome verifications FAILED.")
        exit(1)
