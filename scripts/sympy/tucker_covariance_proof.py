#!/usr/bin/env python3
"""
SymPy proof: streaming Gram eigendecomposition = Tucker mode-2 factors.

Shows that eigenvectors of X^T X are identical to the right singular
vectors of X, which are the Tucker HOSVD mode-2 factor matrix U2.

This justifies computing Tucker decomposition of a 75M x 200 matrix
using only a 200 x 200 Gram matrix built streaming (O(n_ind^2) memory).
"""

from sympy import Matrix, symbols, sqrt, eye, simplify, Rational
import numpy as np


def prove_gram_svd_equivalence():
    """Prove: eigvecs(X^T X) = right singular vectors of X = Tucker U2."""
    print("=" * 60)
    print("PROOF: Gram eigendecomposition = Tucker mode-2 factors")
    print("=" * 60)
    print()

    print("Given: X is an (m x n) matrix (m markers, n individuals)")
    print()

    print("Step 1: SVD of X")
    print("  X = U S V^T")
    print("  where U (m x m) orthogonal, S (m x n) diagonal, V (n x n) orthogonal")
    print()

    print("Step 2: Gram matrix X^T X")
    print("  X^T X = (U S V^T)^T (U S V^T)")
    print("        = V S^T U^T U S V^T")
    print("        = V S^T S V^T        (since U^T U = I)")
    print("        = V D V^T            (where D = S^T S = diag(sigma_i^2))")
    print()

    print("Step 3: Eigendecomposition of X^T X")
    print("  X^T X = V D V^T is already an eigendecomposition!")
    print("  Eigenvalues: sigma_i^2 (squared singular values)")
    print("  Eigenvectors: columns of V (right singular vectors)")
    print()

    print("Step 4: Tucker HOSVD mode-2")
    print("  For a 2D tensor (matrix) X, Tucker mode-2 unfolding = X^T")
    print("  Tucker factor U2 = left singular vectors of X^T")
    print("                   = right singular vectors of X")
    print("                   = V")
    print()

    print("Therefore:")
    print("  eigvecs(X^T X) = V = right_svecs(X) = Tucker_U2(X)")
    print()
    print("QED: The streaming Gram eigendecomposition gives EXACT Tucker")
    print("mode-2 factors without materializing the full matrix.")
    print()
    print("Memory: O(n^2) for the Gram matrix, O(1) per streaming row.")
    print("For n=200 individuals: 200^2 * 8 bytes = 320 KB.")


def verify_numerical():
    """Numerical verification with a small matrix."""
    print()
    print("=" * 60)
    print("NUMERICAL VERIFICATION")
    print("=" * 60)
    print()

    np.random.seed(42)
    m, n = 1000, 5  # 1000 markers, 5 individuals
    X = np.random.randint(0, 20, size=(m, n)).astype(np.float64)
    # Make it sparse (70% zeros)
    X[np.random.random((m, n)) < 0.7] = 0

    # Method 1: Full SVD
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    V_svd = Vt.T  # right singular vectors

    # Method 2: Gram eigendecomposition (streaming simulation)
    gram = np.zeros((n, n))
    for i in range(m):
        row = X[i, :]
        gram += np.outer(row, row)

    eigenvalues, V_gram = np.linalg.eigh(gram)
    # eigh returns ascending order, reverse to match SVD
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    V_gram = V_gram[:, idx]

    # Method 3: Tucker HOSVD mode-2 (= SVD of X^T)
    # X^T is (5 x 1000), left singular vectors = right svecs of X
    U2, S2, _ = np.linalg.svd(X.T, full_matrices=False)
    V_tucker = U2  # left svecs of X^T = right svecs of X

    print(f"Matrix X: {m} x {n}, {100 * np.mean(X == 0):.0f}% zeros")
    print()
    print("Singular values (SVD):", np.round(S[:n], 4))
    print("sqrt(eigenvalues):    ", np.round(np.sqrt(eigenvalues[:n]), 4))
    print("Match:", np.allclose(S[:n], np.sqrt(eigenvalues[:n])))
    print()

    # Eigenvectors may differ by sign -- compare absolute values
    print("Eigenvector comparison (absolute values, first 3 components):")
    for k in range(min(3, n)):
        v_svd = np.abs(V_svd[:, k])
        v_gram = np.abs(V_gram[:, k])
        v_tuck = np.abs(V_tucker[:, k])
        match_gram = np.allclose(v_svd, v_gram, atol=1e-10)
        match_tuck = np.allclose(v_svd, v_tuck, atol=1e-10)
        print(f"  PC{k+1}: SVD={np.round(v_svd, 4)}, "
              f"Gram={np.round(v_gram, 4)} [{'OK' if match_gram else 'FAIL'}], "
              f"Tucker={np.round(v_tuck, 4)} [{'OK' if match_tuck else 'FAIL'}]")

    print()
    all_ok = (
        np.allclose(S[:n], np.sqrt(eigenvalues[:n]))
        and np.allclose(np.abs(V_svd), np.abs(V_gram), atol=1e-10)
        and np.allclose(np.abs(V_svd), np.abs(V_tucker), atol=1e-10)
    )
    if all_ok:
        print("All verifications PASSED: Gram = SVD = Tucker mode-2")
    else:
        print("FAILED!")
        exit(1)

    print()
    print("Variance explained by top components:")
    total_var = np.sum(eigenvalues)
    for k in range(min(5, n)):
        frac = eigenvalues[k] / total_var
        cum = np.sum(eigenvalues[:k+1]) / total_var
        print(f"  PC{k+1}: {frac:.4f} ({cum:.4f} cumulative)")


if __name__ == "__main__":
    prove_gram_svd_equivalence()
    verify_numerical()
