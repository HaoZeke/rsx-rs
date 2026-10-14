module

/-!
# Word-slice exactness

A factor that fits in 64 bits, multiplied by a 64-bit limb and added to a 64-bit
carry, is strictly smaller than `2^128`. The low and high halves of that sum are
therefore an exact split: no rounding and no lost carry. The same counting bound
says a nonempty product of `s`-bit slices fits in `length * s` bits.
-/

namespace OzakiExact

/-- `a < n` is the same statement as `a ≤ n - 1` for positive `n`. -/
theorem le_sub_one_of_lt {a n : Nat} (h : a < n) : a ≤ n - 1 := by
  omega

/-- Strict factors give a strict product, which is the slice-width invariant. -/
theorem mul_lt_of_lt_lt {a b m n : Nat} (ha : a < m) (hb : b < n) (hn : 0 < n) :
    a * b < m * n := by
  have hm : 0 < m := by omega
  calc
    a * b ≤ (m - 1) * b := Nat.mul_le_mul_right _ (le_sub_one_of_lt ha)
    _ ≤ (m - 1) * (n - 1) := Nat.mul_le_mul_left _ (le_sub_one_of_lt hb)
    _ ≤ (m - 1) * n := Nat.mul_le_mul_left _ (Nat.sub_le n 1)
    _ < m * n :=
      Nat.mul_lt_mul_of_pos_right
        (Nat.sub_lt_of_pos_le (by decide : 0 < 1) (by omega : 1 ≤ m)) hn

/-- The schoolbook step `limb * factor + carry` fits in an unsigned 128-bit word. -/
theorem limb_mul_add_lt (a b c : Nat) (ha : a < 2 ^ 64) (hb : b < 2 ^ 64) (hc : c < 2 ^ 64) :
    a * b + c < 2 ^ 128 := by
  have ha' : a ≤ 2 ^ 64 - 1 := le_sub_one_of_lt ha
  have hb' : b ≤ 2 ^ 64 - 1 := le_sub_one_of_lt hb
  have hc' : c ≤ 2 ^ 64 - 1 := le_sub_one_of_lt hc
  have hle : a * b + c ≤ (2 ^ 64 - 1) * (2 ^ 64 - 1) + (2 ^ 64 - 1) :=
    Nat.add_le_add (Nat.mul_le_mul ha' hb') hc'
  have hs : Nat.succ (2 ^ 64 - 1) = 2 ^ 64 := by
    rw [Nat.succ_eq_add_one, Nat.sub_add_cancel (by decide : 1 ≤ 2 ^ 64)]
  have heq : (2 ^ 64 - 1) * (2 ^ 64 - 1) + (2 ^ 64 - 1) = (2 ^ 64 - 1) * 2 ^ 64 := by
    rw [← Nat.mul_succ, hs]
  have hlt : (2 ^ 64 - 1) * 2 ^ 64 < 2 ^ 64 * 2 ^ 64 :=
    Nat.mul_lt_mul_of_pos_right
      (Nat.sub_lt_of_pos_le (by decide : 0 < 1) (by decide : 1 ≤ 2 ^ 64))
      (Nat.two_pow_pos 64)
  have hpow : 2 ^ 64 * 2 ^ 64 = 2 ^ 128 := by
    rw [← Nat.pow_add]
  exact Nat.lt_of_le_of_lt (heq ▸ hle) (hpow ▸ hlt)

/-- Splitting off the low `s` bits is an identity, and the low half is an `s`-bit slice. -/
theorem split_reconstruct (s x : Nat) :
    (x / 2 ^ s) * 2 ^ s + x % 2 ^ s = x ∧ x % 2 ^ s < 2 ^ s := by
  refine ⟨Nat.div_add_mod' x (2 ^ s), Nat.mod_lt x (Nat.two_pow_pos s)⟩

/-- The high half of a limb product is itself a single 64-bit carry. -/
theorem carry_hi_lt (a b c : Nat) (ha : a < 2 ^ 64) (hb : b < 2 ^ 64) (hc : c < 2 ^ 64) :
    (a * b + c) / 2 ^ 64 < 2 ^ 64 := by
  have hwide : a * b + c < 2 ^ 128 := limb_mul_add_lt a b c ha hb hc
  have hpow : 2 ^ 128 = 2 ^ 64 * 2 ^ 64 := by
    rw [← Nat.pow_add]
  exact Nat.div_lt_of_lt_mul (hpow ▸ hwide)

theorem carry_split (a b c : Nat) (ha : a < 2 ^ 64) (hb : b < 2 ^ 64) (hc : c < 2 ^ 64) :
    let wide := a * b + c
    wide / 2 ^ 64 * 2 ^ 64 + wide % 2 ^ 64 = wide ∧ wide / 2 ^ 64 < 2 ^ 64 := by
  refine ⟨Nat.div_add_mod' (a * b + c) (2 ^ 64), carry_hi_lt a b c ha hb hc⟩

/-- Product of a nonempty list of natural numbers. The empty product is unused. -/
def natProd : List Nat → Nat
  | [] => 1
  | x :: xs => x * natProd xs

theorem natProd_cons (x : Nat) (xs : List Nat) : natProd (x :: xs) = x * natProd xs := rfl

theorem natProd_singleton (x : Nat) : natProd [x] = x := by simp [natProd]

/-- A nonempty list of `s`-bit slices multiplies to an integer of at most `length * s` bits. -/
theorem natProd_lt_pow {s : Nat} :
    ∀ xs : List Nat, xs ≠ [] → (∀ x, x ∈ xs → x < 2 ^ s) → natProd xs < 2 ^ (xs.length * s)
  | [], hne, _ => absurd rfl hne
  | [x], _, h => by
      have hx : x < 2 ^ s := h x (by simp)
      simpa [natProd] using hx
  | x :: y :: ys, _, h => by
      have hx : x < 2 ^ s := h x (by simp)
      have htail : natProd (y :: ys) < 2 ^ ((y :: ys).length * s) :=
        natProd_lt_pow (y :: ys) (by simp) (fun z hz => h z (List.mem_cons_of_mem x hz))
      have hmul : x * natProd (y :: ys) < 2 ^ s * 2 ^ ((y :: ys).length * s) :=
        mul_lt_of_lt_lt hx htail (Nat.two_pow_pos _)
      have hlen : (x :: y :: ys).length = (y :: ys).length + 1 := by simp
      have hpow : 2 ^ s * 2 ^ ((y :: ys).length * s) = 2 ^ ((x :: y :: ys).length * s) := by
        rw [← Nat.pow_add, hlen, Nat.add_mul, Nat.one_mul, Nat.add_comm]
      rw [hpow] at hmul
      simpa [natProd] using hmul

/-- If the bit budget `w` covers every slice, the product is an exact `w`-bit integer. -/
theorem natProd_fits {s w : Nat} (xs : List Nat) (hne : xs ≠ [])
    (h : ∀ x, x ∈ xs → x < 2 ^ s) (hw : xs.length * s ≤ w) : natProd xs < 2 ^ w := by
  have hlt := natProd_lt_pow xs hne h
  have hmono : 2 ^ (xs.length * s) ≤ 2 ^ w := Nat.pow_le_pow_right (by decide : 0 < 2) hw
  exact Nat.lt_of_lt_of_le hlt hmono

end OzakiExact
