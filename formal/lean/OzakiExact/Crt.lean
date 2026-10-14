module

import Mathlib.Data.Nat.ChineseRemainder

/-!
# Exact reconstruction from pairwise coprime residues

Ozaki scheme II reduces an integer product modulo several pairwise coprime
moduli and rebuilds it with the Chinese Remainder Theorem. Two residues
determine at most one integer in the half-open range below the modulus product,
so a reconstruction that meets the congruences and the range bound is the
integer product.
-/

namespace OzakiExact

theorem crt_pair_unique (x y m n : Nat) (hcop : Nat.Coprime m n)
    (hx : x < m * n) (hy : y < m * n) (h1 : x ≡ y [MOD m]) (h2 : x ≡ y [MOD n]) : x = y := by
  have hmod : x ≡ y [MOD m * n] := (Nat.modEq_and_modEq_iff_modEq_mul hcop).1 ⟨h1, h2⟩
  exact Nat.ModEq.eq_of_lt_of_lt hmod hx hy

/-- A whole list of pairwise coprime moduli still has a unique representative. -/
theorem crt_list_unique (moduli : List Nat) (x y : Nat) (hcop : moduli.Pairwise Nat.Coprime)
    (hx : x < moduli.prod) (hy : y < moduli.prod)
    (h : ∀ i : Fin moduli.length, x ≡ y [MOD moduli.get i]) : x = y := by
  have hmod : x ≡ y [MOD moduli.prod] := (Nat.modEq_list_prod_iff hcop).2 h
  exact Nat.ModEq.eq_of_lt_of_lt hmod hx hy

end OzakiExact
