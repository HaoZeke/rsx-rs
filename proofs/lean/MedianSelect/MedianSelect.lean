/-!
# Upper median index (rsx `find_median`)

The engine uses the order statistic at index `k = n / 2` (0-based) — the
**upper median** for even length. Replacing a full sort with
`select_nth_unstable(k)` is justified iff that index is the intended
order statistic; these theorems pin `k`.
-/

namespace MedianSelect

def upperMedianIndex (n : Nat) : Nat := n / 2

theorem upperMedianIndex_lt_of_pos (n : Nat) (h : 0 < n) :
    upperMedianIndex n < n := by
  unfold upperMedianIndex
  exact Nat.div_lt_self h (by decide : 1 < 2)

theorem upperMedianIndex_odd (k : Nat) :
    upperMedianIndex (2 * k + 1) = k := by
  unfold upperMedianIndex
  omega

theorem upperMedianIndex_even (k : Nat) :
    upperMedianIndex (2 * k) = k := by
  unfold upperMedianIndex
  -- 2k / 2 = k
  omega

/-- Selection at `upperMedianIndex n` is in-range for any non-empty length. -/
theorem select_index_valid (n : Nat) (h : 0 < n) :
    upperMedianIndex n < n :=
  upperMedianIndex_lt_of_pos n h

end MedianSelect
