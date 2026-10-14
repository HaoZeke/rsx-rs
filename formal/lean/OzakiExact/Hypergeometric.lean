module

import Mathlib.Data.Nat.Choose.Basic
import Mathlib.Tactic.Ring

/-!
# Exact hypergeometric step

The weight of a 2×2 table with fixed margins is a product of two binomial
coefficients. Stepping the observed cell by one multiplies that weight by a
ratio of four positive integers, and the division is exact. That is the
recurrence Fisher's tail walks, with no rounding.
-/

namespace OzakiExact

theorem choose_move_right (n k : Nat) :
    Nat.choose n (k + 1) * (k + 1) = Nat.choose n k * (n - k) :=
  Nat.choose_succ_right_eq n k

theorem choose_move_left (m b : Nat) (hb : 0 < b) :
    Nat.choose m (b - 1) * (m - (b - 1)) = Nat.choose m b * b := by
  have hb' : b - 1 + 1 = b := Nat.sub_add_cancel hb
  simpa [hb'] using (Nat.choose_succ_right_eq m (b - 1)).symm

/-- One step of the Fisher recurrence, written with the four cell factors.

`row1` is the present-margin and `col1` is the group-1 margin. `a` is the
present count in group 1. The identity holds in `Nat`, so every division in
the implementation is exact. -/
theorem fisher_step (row1 col1 N a : Nat) (_ha : a < row1) (hb : a < col1)
    (hroom : col1 - a ≤ N - row1) :
    Nat.choose row1 (a + 1) * Nat.choose (N - row1) (col1 - (a + 1)) *
        ((a + 1) * (N - row1 - (col1 - a) + 1)) =
      Nat.choose row1 a * Nat.choose (N - row1) (col1 - a) *
        ((row1 - a) * (col1 - a)) := by
  set M := N - row1
  set b := col1 - a
  have hbpos : 0 < b := Nat.sub_pos_of_lt hb
  have hindex : col1 - (a + 1) = b - 1 := by
    simp only [b]
    omega
  have hden : M - (b - 1) = M - b + 1 := by
    have hble : b ≤ M := by simpa [M, b] using hroom
    omega
  have h1 : Nat.choose row1 (a + 1) * (a + 1) = Nat.choose row1 a * (row1 - a) :=
    Nat.choose_succ_right_eq row1 a
  have h2 : Nat.choose M (b - 1) * (M - (b - 1)) = Nat.choose M b * b :=
    choose_move_left M b hbpos
  rw [hindex, ← hden]
  calc
    Nat.choose row1 (a + 1) * Nat.choose M (b - 1) * ((a + 1) * (M - (b - 1))) =
        (Nat.choose row1 (a + 1) * (a + 1)) * (Nat.choose M (b - 1) * (M - (b - 1))) := by ring
    _ = (Nat.choose row1 a * (row1 - a)) * (Nat.choose M b * b) := by rw [h1, h2]
    _ = Nat.choose row1 a * Nat.choose M b * ((row1 - a) * b) := by ring

end OzakiExact
