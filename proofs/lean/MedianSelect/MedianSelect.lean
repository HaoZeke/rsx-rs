/-!
# Mathematical median indices (rsx `find_median`)

The mathematical median is the mean of the order statistics at indices
`(n - 1) / 2` and `n / 2` (0-based). They coincide for odd lengths and
identify the adjacent middle values for even lengths.
-/

namespace MedianSelect

def lowerMedianIndex (n : Nat) : Nat := (n - 1) / 2

def upperMedianIndex (n : Nat) : Nat := n / 2

theorem lowerMedianIndex_lt_of_pos (n : Nat) (h : 0 < n) :
    lowerMedianIndex n < n := by
  unfold lowerMedianIndex
  omega

theorem upperMedianIndex_lt_of_pos (n : Nat) (h : 0 < n) :
    upperMedianIndex n < n := by
  unfold upperMedianIndex
  exact Nat.div_lt_self h (by decide : 1 < 2)

theorem upperMedianIndex_odd (k : Nat) :
    upperMedianIndex (2 * k + 1) = k := by
  unfold upperMedianIndex
  omega

theorem lowerMedianIndex_odd (k : Nat) :
    lowerMedianIndex (2 * k + 1) = k := by
  unfold lowerMedianIndex
  omega

theorem upperMedianIndex_even (k : Nat) :
    upperMedianIndex (2 * k) = k := by
  unfold upperMedianIndex
  omega

theorem lowerMedianIndex_even (k : Nat) (h : 0 < k) :
    lowerMedianIndex (2 * k) = k - 1 := by
  unfold lowerMedianIndex
  omega

theorem medianIndices_ordered (n : Nat) :
    lowerMedianIndex n ≤ upperMedianIndex n := by
  unfold lowerMedianIndex upperMedianIndex
  omega

/-- Both selected middle indices are ordered and in range for non-empty input. -/
theorem select_indices_valid (n : Nat) (h : 0 < n) :
    lowerMedianIndex n ≤ upperMedianIndex n ∧
      lowerMedianIndex n < n ∧ upperMedianIndex n < n := by
  exact ⟨medianIndices_ordered n, lowerMedianIndex_lt_of_pos n h,
    upperMedianIndex_lt_of_pos n h⟩

end MedianSelect
