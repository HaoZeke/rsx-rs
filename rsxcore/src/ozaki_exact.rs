// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Exact slice arithmetic for the association tests.
//!
//! The production kernels stay on their current paths. This module is the
//! measured prototype: hypergeometric weights for Fisher's exact test, rising
//! factorials for integer Beta priors, and a correctly rounded rational for
//! the Yates chi-squared statistic. None of them replaces the device kernels,
//! which still evaluate `lgamma`.

use crate::ln_reduced::{self, LN2_BITS};

pub use crate::ln_reduced::LN_POLY_ABS_BOUND;
use crate::stats::BayesFactorModel;

/// Largest positive integer prior accepted by the exact Bayes factor.
/// Every integer in `1..=4096` is an exact `f64`, and RAD-seq Beta priors sit
/// far below this cutoff. Non-integers, including a Jeffreys prior, return
/// [`None`] so the caller keeps the `lgamma` kernel.
pub const MAX_INTEGER_PRIOR: u64 = 4096;

/// A count vector could not be merged.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShardMergeError {
    /// Two shards do not have the same number of loci.
    Shape,
    /// A coordinate sum does not fit in `u64`.
    Overflow,
}

impl std::fmt::Display for ShardMergeError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Shape => write!(formatter, "count shards have different lengths"),
            Self::Overflow => write!(formatter, "merged count does not fit in u64"),
        }
    }
}

impl std::error::Error for ShardMergeError {}

/// Little-endian base-2^64 integer. The zero value has no limbs.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Limbs {
    d: Vec<u64>,
}

impl Limbs {
    fn zero() -> Self {
        Self { d: Vec::new() }
    }

    fn one() -> Self {
        Self { d: vec![1] }
    }

    fn from_u64(value: u64) -> Self {
        if value == 0 {
            Self::zero()
        } else {
            Self { d: vec![value] }
        }
    }

    fn from_u128(value: u128) -> Self {
        let lo = value as u64;
        let hi = (value >> 64) as u64;
        if hi == 0 {
            Self::from_u64(lo)
        } else {
            Self { d: vec![lo, hi] }
        }
    }

    fn is_zero(&self) -> bool {
        self.d.is_empty()
    }

    fn normalize(mut d: Vec<u64>) -> Self {
        while d.last() == Some(&0) {
            d.pop();
        }
        Self { d }
    }

    fn bit_len(&self) -> u32 {
        match self.d.last() {
            None => 0,
            Some(&hi) => {
                let top = 64 - hi.leading_zeros();
                (self.d.len() as u32 - 1) * 64 + top
            }
        }
    }

    fn to_u64(&self) -> Option<u64> {
        match self.d.as_slice() {
            [] => Some(0),
            [limb] => Some(*limb),
            _ => None,
        }
    }

    fn bit(&self, index: u32) -> bool {
        let word = (index / 64) as usize;
        let shift = index % 64;
        self.d
            .get(word)
            .is_some_and(|limb| ((limb >> shift) & 1) == 1)
    }

    fn cmp_limbs(&self, other: &Self) -> std::cmp::Ordering {
        match self.d.len().cmp(&other.d.len()) {
            std::cmp::Ordering::Equal => {
                for (left, right) in self.d.iter().rev().zip(other.d.iter().rev()) {
                    match left.cmp(right) {
                        std::cmp::Ordering::Equal => {}
                        order => return order,
                    }
                }
                std::cmp::Ordering::Equal
            }
            order => order,
        }
    }

    fn add(&self, other: &Self) -> Self {
        let n = self.d.len().max(other.d.len());
        let mut d = Vec::with_capacity(n + 1);
        let mut carry = 0u64;
        for i in 0..n {
            let left = u128::from(*self.d.get(i).unwrap_or(&0));
            let right = u128::from(*other.d.get(i).unwrap_or(&0));
            let sum = left + right + u128::from(carry);
            d.push(sum as u64);
            carry = (sum >> 64) as u64;
        }
        if carry != 0 {
            d.push(carry);
        }
        Self { d }
    }

    fn add_u64(&self, value: u64) -> Self {
        if value == 0 {
            return self.clone();
        }
        self.add(&Self::from_u64(value))
    }

    fn sub(&self, other: &Self) -> Self {
        debug_assert!(self.cmp_limbs(other) != std::cmp::Ordering::Less);
        let mut d = Vec::with_capacity(self.d.len());
        let mut borrow = 0u64;
        for i in 0..self.d.len() {
            let left = u128::from(self.d[i]);
            let right = u128::from(*other.d.get(i).unwrap_or(&0)) + u128::from(borrow);
            if left >= right {
                d.push((left - right) as u64);
                borrow = 0;
            } else {
                d.push(((1u128 << 64) + left - right) as u64);
                borrow = 1;
            }
        }
        debug_assert_eq!(borrow, 0);
        Self::normalize(d)
    }

    /// `limb * factor + carry` fits in 128 bits, so the high half is a `u64`.
    fn mul_u64(&self, factor: u64) -> Self {
        if factor == 0 || self.is_zero() {
            return Self::zero();
        }
        let mut d = Vec::with_capacity(self.d.len() + 1);
        let mut carry = 0u64;
        for &limb in &self.d {
            let prod = u128::from(limb) * u128::from(factor) + u128::from(carry);
            d.push(prod as u64);
            carry = (prod >> 64) as u64;
        }
        if carry != 0 {
            d.push(carry);
        }
        Self { d }
    }

    fn mul(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return Self::zero();
        }
        let mut acc = vec![0u64; self.d.len() + other.d.len()];
        for (i, &left) in self.d.iter().enumerate() {
            let mut carry = 0u128;
            for (j, &right) in other.d.iter().enumerate() {
                let cur = u128::from(acc[i + j]) + u128::from(left) * u128::from(right) + carry;
                acc[i + j] = cur as u64;
                carry = cur >> 64;
            }
            let mut k = i + other.d.len();
            while carry > 0 {
                if k == acc.len() {
                    acc.push(0);
                }
                let cur = u128::from(acc[k]) + carry;
                acc[k] = cur as u64;
                carry = cur >> 64;
                k += 1;
            }
        }
        Self::normalize(acc)
    }

    fn shl(&self, bits: u32) -> Self {
        if self.is_zero() || bits == 0 {
            return self.clone();
        }
        let words = (bits / 64) as usize;
        let shift = bits % 64;
        let mut d = vec![0u64; words];
        if shift == 0 {
            d.extend_from_slice(&self.d);
            return Self { d };
        }
        let mut carry = 0u64;
        for &limb in &self.d {
            let wide = (u128::from(limb) << shift) | u128::from(carry);
            d.push(wide as u64);
            carry = (wide >> 64) as u64;
        }
        if carry != 0 {
            d.push(carry);
        }
        Self { d }
    }

    fn div_u64(&self, divisor: u64) -> (Self, u64) {
        assert!(divisor != 0, "division by zero");
        if self.is_zero() {
            return (Self::zero(), 0);
        }
        let mut quotient = vec![0u64; self.d.len()];
        let mut remainder = 0u128;
        for (index, &limb) in self.d.iter().enumerate().rev() {
            let cur = (remainder << 64) | u128::from(limb);
            quotient[index] = (cur / u128::from(divisor)) as u64;
            remainder = cur % u128::from(divisor);
        }
        (Self::normalize(quotient), remainder as u64)
    }

    fn div_exact_u64(&self, divisor: u64) -> Self {
        let (quotient, remainder) = self.div_u64(divisor);
        debug_assert_eq!(remainder, 0, "expected an exact division");
        quotient
    }

    fn rem_u64(&self, modulus: u64) -> u64 {
        if modulus == 0 || self.is_zero() {
            return 0;
        }
        let mut remainder = 0u128;
        for &limb in self.d.iter().rev() {
            remainder = ((remainder << 64) | u128::from(limb)) % u128::from(modulus);
        }
        remainder as u64
    }

    /// Binary long division. The quotient used for rounding has fewer than 56 bits;
    /// general products here are a few hundred bits, so a bit loop is enough.
    fn div_mod(numerator: &Self, denominator: &Self) -> (Self, Self) {
        assert!(!denominator.is_zero(), "division by zero");
        if numerator.cmp_limbs(denominator) == std::cmp::Ordering::Less {
            return (Self::zero(), numerator.clone());
        }
        let mut quotient = Self::zero();
        let mut remainder = Self::zero();
        let bits = numerator.bit_len();
        for index in (0..bits).rev() {
            remainder = remainder.shl(1);
            if numerator.bit(index) {
                remainder = remainder.add_u64(1);
            }
            if remainder.cmp_limbs(denominator) != std::cmp::Ordering::Less {
                remainder = remainder.sub(denominator);
                let word = (index / 64) as usize;
                let shift = index % 64;
                if quotient.d.len() <= word {
                    quotient.d.resize(word + 1, 0);
                }
                quotient.d[word] |= 1u64 << shift;
            }
        }
        (quotient, remainder)
    }
}

fn gcd_u64(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let next = a % b;
        a = b;
        b = next;
    }
    a
}

fn gcd_u128(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let next = a % b;
        a = b;
        b = next;
    }
    a
}

fn mod_inverse(value: u64, modulus: u64) -> Option<u64> {
    if modulus <= 1 {
        return None;
    }
    let mut t_curr: i128 = 0;
    let mut t_next: i128 = 1;
    let mut r_curr: i128 = i128::from(modulus);
    let mut r_next: i128 = i128::from(value % modulus);
    while r_next != 0 {
        let q = r_curr / r_next;
        let next_t = t_curr - q * t_next;
        t_curr = t_next;
        t_next = next_t;
        let next_r = r_curr - q * r_next;
        r_curr = r_next;
        r_next = next_r;
    }
    if r_curr != 1 {
        return None;
    }
    if t_curr < 0 {
        t_curr += i128::from(modulus);
    }
    u64::try_from(t_curr).ok()
}

/// Correctly rounded `numerator / denominator` in round-to-nearest, ties to even.
///
/// The scaling puts the integer quotient in `[2^53, 2^55)`. One extra bit is
/// the round bit; anything below it, including a nonzero remainder, is sticky.
fn ratio_to_f64(numerator: &Limbs, denominator: &Limbs) -> f64 {
    if denominator.is_zero() {
        return f64::NAN;
    }
    if numerator.is_zero() {
        return 0.0;
    }
    let num_bits = i64::from(numerator.bit_len());
    let den_bits = i64::from(denominator.bit_len());
    let diff = num_bits - den_bits;
    // Strictly below the tie between 0 and the smallest subnormal.
    if diff <= -1076 {
        return 0.0;
    }
    // Strictly above 2^1024.
    if diff > 1024 {
        return f64::INFINITY;
    }
    let k = 54 + den_bits - num_bits;
    let (quotient, remainder) = if k >= 0 {
        Limbs::div_mod(&numerator.shl(k as u32), denominator)
    } else {
        Limbs::div_mod(numerator, &denominator.shl((-k) as u32))
    };
    let Some(q) = quotient.to_u64() else {
        return f64::NAN;
    };
    debug_assert!(
        ((1u64 << 53)..(1u64 << 55)).contains(&q),
        "scaled quotient {q} left the proved binade"
    );
    round_scaled(q, k, !remainder.is_zero())
}

fn round_scaled(q: u64, k: i64, sticky: bool) -> f64 {
    let adjust = i64::from(u8::from(q >= (1u64 << 54)));
    let exponent = 53 - k + adjust;
    if exponent < -1022 {
        round_subnormal(q, k, sticky)
    } else {
        round_normal(q, k, sticky)
    }
}

fn round_normal(mut q: u64, mut k: i64, mut sticky: bool) -> f64 {
    if q >= (1u64 << 54) {
        sticky |= (q & 1) != 0;
        q >>= 1;
        k -= 1;
    }
    let mut significand = q >> 1;
    let round_bit = (q & 1) != 0;
    let mut exponent = 53 - k;
    if round_bit && (sticky || (significand & 1) == 1) {
        significand += 1;
    }
    if significand == (1u64 << 53) {
        significand = 1u64 << 52;
        exponent += 1;
    }
    if exponent > 1023 {
        return f64::INFINITY;
    }
    if exponent < -1022 {
        return 0.0;
    }
    let biased = (exponent + 1023) as u64;
    let fraction = significand - (1u64 << 52);
    f64::from_bits((biased << 52) | fraction)
}

fn round_subnormal(q: u64, k: i64, sticky: bool) -> f64 {
    let shift = k - 1074;
    if !(1..64).contains(&shift) {
        return 0.0;
    }
    let shift = shift as u32;
    let mut significand = q >> shift;
    let lost = q & ((1u64 << shift) - 1);
    let half = 1u64 << (shift - 1);
    let round_bit = (lost & half) != 0;
    let sticky_bits = (lost & (half - 1)) != 0 || sticky;
    if round_bit && (sticky_bits || (significand & 1) == 1) {
        significand += 1;
    }
    if significand == 0 {
        return 0.0;
    }
    if significand == (1u64 << 52) {
        // Rounded up to the smallest normal, 2^-1022.
        return f64::from_bits(1u64 << 52);
    }
    debug_assert!(significand < (1u64 << 52));
    f64::from_bits(significand)
}

fn rational_u128_to_f64(numerator: u128, denominator: u128) -> f64 {
    if numerator <= (1u128 << 53) && denominator <= (1u128 << 53) {
        return (numerator as f64) / (denominator as f64);
    }
    ratio_to_f64(&Limbs::from_u128(numerator), &Limbs::from_u128(denominator))
}

/// Natural log of a positive rational, via a correctly rounded significand
/// and the certified panel polynomial on `[1, 2)`.
///
/// The polynomial error is [`ln_reduced::LN_POLY_ABS_BOUND`]. Rounding the
/// significand adds at most about `2^-53` before the polynomial, and `e * ln 2`
/// is one `f64` product. This is not a correctly rounded `ln`.
fn ln_ratio(numerator: &Limbs, denominator: &Limbs) -> f64 {
    let rounded = ratio_to_f64(numerator, denominator);
    if rounded == 0.0 {
        return f64::NEG_INFINITY;
    }
    if rounded.is_infinite() {
        return f64::INFINITY;
    }
    if rounded.is_nan() || rounded < 0.0 {
        return f64::NAN;
    }
    let bits = rounded.to_bits();
    let biased = ((bits >> 52) & 0x7ff) as i32;
    let fraction = bits & ((1u64 << 52) - 1);
    let ln2 = f64::from_bits(LN2_BITS);
    if biased == 0 {
        if fraction == 0 {
            return f64::NEG_INFINITY;
        }
        let exp = 63 - fraction.leading_zeros() as i32;
        let scale = (fraction as f64) / ((1u64 << (exp as u32)) as f64);
        return ln_reduced::ln_reduced(scale) + f64::from(exp - 1074) * ln2;
    }
    let significand = f64::from_bits((1023u64 << 52) | fraction);
    let exponent = biased - 1023;
    ln_reduced::ln_reduced(significand) + f64::from(exponent) * ln2
}

#[derive(Clone, Copy)]
struct Table {
    a: u64,
    row1: u64,
    col1: u64,
    n: u64,
}

fn contingency(n_g1: u32, n_g2: u32, total_g1: u32, total_g2: u32) -> Option<Table> {
    total_g1.checked_sub(n_g1)?;
    total_g2.checked_sub(n_g2)?;
    let a = u64::from(n_g1);
    let c = u64::from(n_g2);
    Some(Table {
        a,
        row1: a + c,
        col1: u64::from(total_g1),
        n: u64::from(total_g1) + u64::from(total_g2),
    })
}

fn support(table: &Table) -> (u64, u64) {
    let absent = table.n - table.row1;
    let min_a = table.col1.saturating_sub(absent);
    let max_a = table.row1.min(table.col1);
    (min_a, max_a)
}

fn step_factors(ai: u64, table: &Table) -> (u64, u64, u64, u64) {
    let group1_absent = table.col1 - ai;
    let present_remaining = table.row1 - ai;
    let absent_margin = table.n - table.row1;
    let group2_absent = absent_margin - group1_absent;
    (present_remaining, group1_absent, ai + 1, group2_absent + 1)
}

fn binom_u128(n: u64, k: u64) -> Option<u128> {
    if k > n {
        return Some(0);
    }
    let k = k.min(n - k);
    let mut acc = 1u128;
    for i in 0..k {
        let mut num = n - i;
        let mut den = i + 1;
        let g = gcd_u64(num, den);
        num /= g;
        den /= g;
        let g2 = gcd_u128(acc, u128::from(den)) as u64;
        acc /= u128::from(g2);
        den /= g2;
        acc = acc.checked_mul(u128::from(num))?;
        if den != 1 {
            if acc % u128::from(den) != 0 {
                return None;
            }
            acc /= u128::from(den);
        }
    }
    Some(acc)
}

fn step_u128(weight: u128, ai: u64, table: &Table) -> Option<u128> {
    let (mut n1, mut n2, mut d1, mut d2) = step_factors(ai, table);
    let mut weight = weight;
    cancel_u64(&mut n1, &mut d1);
    cancel_u64(&mut n1, &mut d2);
    cancel_u64(&mut n2, &mut d1);
    cancel_u64(&mut n2, &mut d2);
    cancel_weight(&mut weight, &mut d1);
    cancel_weight(&mut weight, &mut d2);
    weight = weight.checked_mul(u128::from(n1))?;
    weight = weight.checked_mul(u128::from(n2))?;
    if d1 != 1 {
        if weight % u128::from(d1) != 0 {
            return None;
        }
        weight /= u128::from(d1);
    }
    if d2 != 1 {
        if weight % u128::from(d2) != 0 {
            return None;
        }
        weight /= u128::from(d2);
    }
    Some(weight)
}

fn cancel_u64(numer: &mut u64, denom: &mut u64) {
    let g = gcd_u64(*numer, *denom);
    *numer /= g;
    *denom /= g;
}

fn cancel_weight(weight: &mut u128, denom: &mut u64) {
    if *denom == 1 {
        return;
    }
    let g = gcd_u128(*weight, u128::from(*denom)) as u64;
    *weight /= u128::from(g);
    *denom /= g;
}

fn try_fisher_u128(table: &Table) -> Option<f64> {
    let (min_a, max_a) = support(table);
    if min_a > max_a || table.a < min_a || table.a > max_a {
        return Some(1.0);
    }
    let absent = table.n - table.row1;
    let mut weight =
        binom_u128(table.row1, min_a)?.checked_mul(binom_u128(absent, table.col1 - min_a)?)?;
    let mut observed = None;
    let mut ai = min_a;
    loop {
        if ai == table.a {
            observed = Some(weight);
        }
        if ai == max_a {
            break;
        }
        weight = step_u128(weight, ai, table)?;
        ai += 1;
    }
    let observed = observed?;
    weight = binom_u128(table.row1, min_a)?.checked_mul(binom_u128(absent, table.col1 - min_a)?)?;
    let mut total = 0u128;
    let mut tail = 0u128;
    ai = min_a;
    loop {
        total = total.checked_add(weight)?;
        if weight <= observed {
            tail = tail.checked_add(weight)?;
        }
        if ai == max_a {
            break;
        }
        weight = step_u128(weight, ai, table)?;
        ai += 1;
    }
    if total == 0 {
        return Some(1.0);
    }
    Some(rational_u128_to_f64(tail, total))
}

fn binom_limbs(n: u64, k: u64) -> Limbs {
    if k > n {
        return Limbs::zero();
    }
    let k = k.min(n - k);
    let mut acc = Limbs::one();
    for i in 0..k {
        let mut num = n - i;
        let mut den = i + 1;
        let g = gcd_u64(num, den);
        num /= g;
        den /= g;
        let g2 = gcd_u64(acc.rem_u64(den), den);
        acc = acc.div_exact_u64(g2);
        den /= g2;
        acc = acc.mul_u64(num);
        if den != 1 {
            acc = acc.div_exact_u64(den);
        }
    }
    acc
}

fn step_limbs(weight: &Limbs, ai: u64, table: &Table) -> Limbs {
    let (n1, n2, d1, d2) = step_factors(ai, table);
    let wide = weight.mul_u64(n1).mul_u64(n2);
    let (after_first, rem1) = wide.div_u64(d1);
    debug_assert_eq!(rem1, 0);
    let (after_second, rem2) = after_first.div_u64(d2);
    debug_assert_eq!(rem2, 0);
    after_second
}

fn endpoint_limbs(table: &Table, ai: u64) -> Limbs {
    let absent = table.n - table.row1;
    binom_limbs(table.row1, ai).mul(&binom_limbs(absent, table.col1 - ai))
}

fn fisher_limbs(table: &Table) -> f64 {
    let (min_a, max_a) = support(table);
    if min_a > max_a || table.a < min_a || table.a > max_a {
        return 1.0;
    }
    let observed = endpoint_limbs(table, table.a);
    let mut weight = endpoint_limbs(table, min_a);
    let mut total = Limbs::zero();
    let mut tail = Limbs::zero();
    let mut ai = min_a;
    loop {
        total = total.add(&weight);
        if weight.cmp_limbs(&observed) != std::cmp::Ordering::Greater {
            tail = tail.add(&weight);
        }
        if ai == max_a {
            break;
        }
        weight = step_limbs(&weight, ai, table);
        ai += 1;
    }
    if total.is_zero() {
        return 1.0;
    }
    ratio_to_f64(&tail, &total)
}

/// Two-sided Fisher probability, as an exact ratio of hypergeometric weights
/// rounded once to `f64`.
///
/// Invalid margins (`count > total`) return `1.0`, matching [`crate::stats::fisher_exact`].
/// The probability is not clamped; [`fisher_p_value`] applies the `1e-16` floor.
#[must_use]
pub fn fisher_probability(n_g1: u32, n_g2: u32, total_g1: u32, total_g2: u32) -> f64 {
    let Some(table) = contingency(n_g1, n_g2, total_g1, total_g2) else {
        return 1.0;
    };
    if let Some(probability) = try_fisher_u128(&table) {
        return probability;
    }
    fisher_limbs(&table)
}

/// [`fisher_probability`] clamped to `[1e-16, 1]`, the range the association
/// commands report.
#[must_use]
pub fn fisher_p_value(n_g1: u32, n_g2: u32, total_g1: u32, total_g2: u32) -> f64 {
    fisher_probability(n_g1, n_g2, total_g1, total_g2).clamp(1e-16, 1.0)
}

struct Ratio {
    num: Limbs,
    den: Limbs,
}

impl Ratio {
    fn one() -> Self {
        Self {
            num: Limbs::one(),
            den: Limbs::one(),
        }
    }

    fn mul_u64(&mut self, mut factor: u64) {
        if factor == 0 {
            self.num = Limbs::zero();
            return;
        }
        if self.den.is_zero() {
            return;
        }
        let g = gcd_u64(factor, self.den.rem_u64(factor));
        factor /= g;
        self.den = self.den.div_exact_u64(g);
        self.num = self.num.mul_u64(factor);
    }

    fn div_u64(&mut self, mut divisor: u64) {
        if divisor == 0 || self.num.is_zero() {
            self.num = Limbs::zero();
            return;
        }
        let g = gcd_u64(divisor, self.num.rem_u64(divisor));
        divisor /= g;
        self.num = self.num.div_exact_u64(g);
        self.den = self.den.mul_u64(divisor);
    }

    fn mul_rising(&mut self, start: u64, len: u64) {
        for offset in 0..len {
            self.mul_u64(start + offset);
        }
    }

    fn div_rising(&mut self, start: u64, len: u64) {
        for offset in 0..len {
            self.div_u64(start + offset);
        }
    }

    fn to_f64(&self) -> f64 {
        ratio_to_f64(&self.num, &self.den)
    }
}

struct IntPriors {
    a1: u64,
    b1: u64,
    a2: u64,
    b2: u64,
    a0: u64,
    b0: u64,
}

fn integer_prior(value: f64) -> Option<u64> {
    if !value.is_finite() || value < 1.0 || value > MAX_INTEGER_PRIOR as f64 {
        return None;
    }
    let rounded = value as u64;
    if (rounded as f64) != value {
        return None;
    }
    Some(rounded)
}

fn integer_priors(model: &BayesFactorModel) -> Option<IntPriors> {
    Some(IntPriors {
        a1: integer_prior(model.alternative_group1.alpha)?,
        b1: integer_prior(model.alternative_group1.beta)?,
        a2: integer_prior(model.alternative_group2.alpha)?,
        b2: integer_prior(model.alternative_group2.beta)?,
        a0: integer_prior(model.null.alpha)?,
        b0: integer_prior(model.null.beta)?,
    })
}

fn bayes_uniform_u128(k1: u64, n1: u64, k2: u64, n2: u64) -> Option<f64> {
    let n = n1 + n2;
    let k = k1 + k2;
    let mut num = binom_u128(n, k)?;
    let c1 = binom_u128(n1, k1)?;
    let c2 = binom_u128(n2, k2)?;
    let mut den = c1.checked_mul(c2)?;
    let g = gcd_u128(num, den);
    num /= g;
    den /= g;
    num = num.checked_mul(u128::from(n) + 1)?;
    den = den.checked_mul(u128::from(n1) + 1)?;
    den = den.checked_mul(u128::from(n2) + 1)?;
    let g = gcd_u128(num, den);
    num /= g;
    den /= g;
    if den == 0 {
        return None;
    }
    Some(rational_u128_to_f64(num, den))
}

fn bayes_rising(k1: u64, n1: u64, k2: u64, n2: u64, priors: &IntPriors) -> Ratio {
    let mut ratio = Ratio::one();
    ratio.mul_rising(priors.a1, k1);
    ratio.mul_rising(priors.b1, n1 - k1);
    ratio.div_rising(priors.a1 + priors.b1, n1);
    ratio.mul_rising(priors.a2, k2);
    ratio.mul_rising(priors.b2, n2 - k2);
    ratio.div_rising(priors.a2 + priors.b2, n2);
    let k = k1 + k2;
    let n = n1 + n2;
    ratio.div_rising(priors.a0, k);
    ratio.div_rising(priors.b0, n - k);
    ratio.mul_rising(priors.a0 + priors.b0, n);
    ratio
}

fn bayes_ratio(
    n_g1: u32,
    n_g2: u32,
    total_g1: u32,
    total_g2: u32,
    model: &BayesFactorModel,
) -> Option<Ratio> {
    if n_g1 > total_g1 || n_g2 > total_g2 {
        return None;
    }
    let priors = integer_priors(model)?;
    let k1 = u64::from(n_g1);
    let n1 = u64::from(total_g1);
    let k2 = u64::from(n_g2);
    let n2 = u64::from(total_g2);
    Some(bayes_rising(k1, n1, k2, n2, &priors))
}

/// Bayes factor for positive integer Beta priors at most [`MAX_INTEGER_PRIOR`].
///
/// Returns [`None`] when a prior is not such an integer or a count exceeds its
/// margin. The uniform `Beta(1, 1)` case is the binomial form
/// `(N+1)/((n1+1)(n2+1)) * C(N, K)/(C(n1, k1) C(n2, k2))`, which is the
/// marginal the `lgamma` kernel computes (the binomial coefficient of the
/// likelihood is not part of that marginal).
#[must_use]
pub fn bayes_factor_integer_prior(
    n_g1: u32,
    n_g2: u32,
    total_g1: u32,
    total_g2: u32,
    model: &BayesFactorModel,
) -> Option<f64> {
    if n_g1 > total_g1 || n_g2 > total_g2 {
        return None;
    }
    let priors = integer_priors(model)?;
    let k1 = u64::from(n_g1);
    let n1 = u64::from(total_g1);
    let k2 = u64::from(n_g2);
    let n2 = u64::from(total_g2);
    let uniform = priors.a1 == 1
        && priors.b1 == 1
        && priors.a2 == 1
        && priors.b2 == 1
        && priors.a0 == 1
        && priors.b0 == 1;
    if uniform {
        if let Some(value) = bayes_uniform_u128(k1, n1, k2, n2) {
            return Some(value);
        }
    }
    Some(bayes_rising(k1, n1, k2, n2, &priors).to_f64())
}

/// Natural log of [`bayes_factor_integer_prior`], with the same [`None`] cases.
#[must_use]
pub fn ln_bayes_factor_integer_prior(
    n_g1: u32,
    n_g2: u32,
    total_g1: u32,
    total_g2: u32,
    model: &BayesFactorModel,
) -> Option<f64> {
    let ratio = bayes_ratio(n_g1, n_g2, total_g1, total_g2, model)?;
    Some(ln_ratio(&ratio.num, &ratio.den))
}

/// Yates chi-squared statistic as one correctly rounded rational.
///
/// A zero margin returns `NaN`, matching [`crate::stats::chi_squared_yates`].
/// For biological margins the cross-products fit in the `f64` mantissa and the
/// two functions are bit-identical; past that range this one is the rounded
/// exact value and the `f64` formula is not.
#[must_use]
pub fn chi_squared_yates_rounded(
    n_group1: u32,
    n_group2: u32,
    total_group1: u32,
    total_group2: u32,
) -> f64 {
    if n_group1 > total_group1 || n_group2 > total_group2 {
        return f64::NAN;
    }
    let t1 = u64::from(total_group1);
    let t2 = u64::from(total_group2);
    let present = u64::from(n_group1) + u64::from(n_group2);
    let n = t1 + t2;
    let absent = n - present;
    if t1 == 0 || t2 == 0 || present == 0 || absent == 0 {
        return f64::NAN;
    }
    let ad = u64::from(n_group1) * t2;
    let bc = u64::from(n_group2) * t1;
    let delta = ad.abs_diff(bc);
    if u128::from(delta) * 2 <= u128::from(n) {
        return 0.0;
    }
    let mut den = Limbs::from_u64(t1)
        .mul_u64(t2)
        .mul_u64(present)
        .mul_u64(absent);
    let num = if n % 2 == 0 {
        let y = delta - n / 2;
        Limbs::from_u64(y).mul_u64(y).mul_u64(n)
    } else {
        let base = u128::from(delta) * 2 - u128::from(n);
        den = den.mul_u64(4);
        Limbs::from_u128(base)
            .mul(&Limbs::from_u128(base))
            .mul_u64(n)
    };
    ratio_to_f64(&num, &den)
}

/// Pointwise sum of equal-length count shards.
///
/// Addition of bounded naturals is a commutative monoid on the values that fit
/// in `u64`. Any binary reduction tree agrees with the left fold, and a sum
/// that does not fit returns [`ShardMergeError::Overflow`] rather than wrapping.
pub fn merge_count_shards(shards: &[Vec<u64>]) -> Result<Vec<u64>, ShardMergeError> {
    let Some(first) = shards.first() else {
        return Ok(Vec::new());
    };
    let width = first.len();
    let mut acc = vec![0u64; width];
    for shard in shards {
        if shard.len() != width {
            return Err(ShardMergeError::Shape);
        }
        for (slot, &count) in acc.iter_mut().zip(shard.iter()) {
            *slot = slot.checked_add(count).ok_or(ShardMergeError::Overflow)?;
        }
    }
    Ok(acc)
}

/// Line-aligned slices of a byte buffer.
///
/// A chunk ends on a newline or at the end of the buffer, and the concatenation
/// of the chunks is the whole buffer. `target == 0` still terminates: the scan
/// for the next newline always consumes at least one byte.
#[must_use]
pub fn partition_line_chunks(data: &[u8], target: usize) -> Vec<std::ops::Range<usize>> {
    let mut chunks = Vec::new();
    let mut start = 0usize;
    while start < data.len() {
        let mut end = match start.checked_add(target) {
            Some(sum) => sum.min(data.len()),
            None => data.len(),
        };
        if end < data.len() {
            if let Some(rel) = memchr::memchr(b'\n', &data[end..]) {
                end = end + rel + 1;
            } else {
                end = data.len();
            }
        }
        if end > start {
            chunks.push(start..end);
            start = end;
        } else {
            // No progress is a bug in the newline scan. Consume one byte.
            chunks.push(start..start + 1);
            start += 1;
        }
    }
    chunks
}

/// Chinese-remainder reconstruction.
///
/// `residues[i]` is taken modulo `moduli[i]`. The moduli must be pairwise
/// coprime and nonzero. The result is the unique representative in
/// `[0, product)`.
#[must_use]
pub fn crt_reconstruct(residues: &[u64], moduli: &[u64]) -> Option<Vec<u64>> {
    if residues.len() != moduli.len() || moduli.is_empty() {
        return None;
    }
    if moduli.iter().any(|&modulus| modulus == 0) {
        return None;
    }
    if residues
        .iter()
        .zip(moduli.iter())
        .any(|(&residue, &modulus)| residue >= modulus)
    {
        return None;
    }
    for i in 0..moduli.len() {
        for j in 0..i {
            if gcd_u64(moduli[i], moduli[j]) != 1 {
                return None;
            }
        }
    }
    let mut value = Limbs::from_u64(residues[0]);
    let mut product = Limbs::from_u64(moduli[0]);
    for i in 1..moduli.len() {
        let modulus = moduli[i];
        let current = value.rem_u64(modulus);
        let diff = (residues[i] + modulus - current) % modulus;
        let prod_mod = product.rem_u64(modulus);
        let inverse = mod_inverse(prod_mod, modulus)?;
        let term = ((u128::from(diff) * u128::from(inverse)) % u128::from(modulus)) as u64;
        value = value.add(&product.mul_u64(term));
        product = product.mul_u64(modulus);
    }
    Some(value.d)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stats::{self, BayesFactorModel, BetaPrior};

    fn ulp_diff(left: f64, right: f64) -> u64 {
        left.to_bits().abs_diff(right.to_bits())
    }

    fn binom_ref(n: u64, k: u64) -> u128 {
        assert!(k <= n);
        let k = k.min(n - k);
        let mut acc = 1u128;
        for i in 0..k {
            acc = acc * u128::from(n - i) / u128::from(i + 1);
        }
        acc
    }

    #[test]
    fn limb_product_splits_inside_u128() {
        let factor = u64::MAX;
        let wide = u128::from(factor) * u128::from(factor) + u128::from(factor);
        assert!(wide < u128::MAX);
        assert_eq!(wide >> 64, u128::from(u64::MAX));
        assert_eq!(wide as u64, 0);
        let limbs = Limbs::from_u64(factor).mul_u64(factor).add_u64(factor);
        assert_eq!(limbs.d, vec![0, u64::MAX]);
    }

    #[test]
    fn ratio_rounds_ties_and_subnormals() {
        let cases = [
            (1u128, 1u128, 1.0f64),
            (1, 2, 0.5),
            (3, 2, 1.5),
            (1, 3, 1.0 / 3.0),
            (2, 3, 2.0 / 3.0),
        ];
        for (num, den, expected) in cases {
            assert_eq!(
                rational_u128_to_f64(num, den).to_bits(),
                expected.to_bits(),
                "{num}/{den}"
            );
        }
        let tie = (1u128 << 53) + 1;
        assert_eq!(
            rational_u128_to_f64(tie, 1).to_bits(),
            (2.0f64.powi(53)).to_bits()
        );
        let up = (1u128 << 53) + 3;
        assert_eq!(
            rational_u128_to_f64(up, 1).to_bits(),
            (2.0f64.powi(53) + 4.0).to_bits()
        );
        let min_sub = Limbs::one();
        let den = Limbs::one().shl(1074);
        assert_eq!(ratio_to_f64(&min_sub, &den).to_bits(), 1);
        let tie_under = Limbs::one().shl(1075);
        assert_eq!(ratio_to_f64(&min_sub, &tie_under).to_bits(), 0);
        assert_eq!(
            ratio_to_f64(&Limbs::one(), &Limbs::one().shl(1076)).to_bits(),
            0
        );
        assert!(ratio_to_f64(&Limbs::one().shl(1024), &Limbs::one()).is_infinite());
    }

    #[test]
    fn small_integers_match_ieee_division() {
        for num in 1u128..200 {
            for den in 1u128..200 {
                let got = rational_u128_to_f64(num, den);
                let expect = (num as f64) / (den as f64);
                assert_eq!(got.to_bits(), expect.to_bits(), "{num}/{den}");
            }
        }
    }

    #[test]
    fn fisher_matches_integer_weights_for_every_small_table() {
        for t1 in 0..9 {
            for t2 in 0..9 {
                for n1 in 0..=t1 {
                    for n2 in 0..=t2 {
                        let got = fisher_probability(n1, n2, t1, t2);
                        let Some(table) = contingency(n1, n2, t1, t2) else {
                            continue;
                        };
                        let (min_a, max_a) = support(&table);
                        let absent = table.n - table.row1;
                        let mut total = 0u128;
                        let mut tail = 0u128;
                        let observed = binom_ref(table.row1, table.a)
                            * binom_ref(absent, table.col1 - table.a);
                        for ai in min_a..=max_a {
                            let weight =
                                binom_ref(table.row1, ai) * binom_ref(absent, table.col1 - ai);
                            total += weight;
                            if weight <= observed {
                                tail += weight;
                            }
                        }
                        let expect = if total == 0 {
                            1.0
                        } else {
                            rational_u128_to_f64(tail, total)
                        };
                        assert_eq!(got.to_bits(), expect.to_bits(), "{n1},{n2},{t1},{t2}");
                    }
                }
            }
        }
    }

    #[test]
    fn fisher_u128_and_limbs_agree_past_the_fast_path() {
        // C(80, 40) fits in u128; the limb walk is the fallback's arithmetic.
        // Central weight C(50, 25)^2 still fits in u128; C(160, 80)^2 does not.
        let table = contingency(25, 25, 50, 50).unwrap();
        let fast = try_fisher_u128(&table).unwrap();
        let slow = fisher_limbs(&table);
        assert_eq!(fast.to_bits(), slow.to_bits());
        let wide = contingency(80, 80, 160, 160).unwrap();
        assert!(try_fisher_u128(&wide).is_none());
        let p = fisher_limbs(&wide);
        assert!(p > 0.0 && p <= 1.0);
    }

    #[test]
    fn fisher_invalid_and_floor() {
        assert_eq!(fisher_probability(3, 0, 2, 2), 1.0);
        let p = fisher_p_value(20, 0, 20, 20);
        assert!((1e-16..=1.0).contains(&p));
    }

    #[test]
    fn uniform_bayes_is_one_and_a_half_and_matches_rising_factorials() {
        let model = BayesFactorModel::uniform_v1();
        let got = bayes_factor_integer_prior(1, 0, 1, 1, &model).unwrap();
        assert_eq!(got.to_bits(), 1.5f64.to_bits());
        let forced = bayes_rising(1, 1, 0, 1, &integer_priors(&model).unwrap()).to_f64();
        assert_eq!(got.to_bits(), forced.to_bits());

        for n in [4u32, 10, 20, 30] {
            let bf = bayes_factor_integer_prior(n / 2, n / 2, n, n, &model).unwrap();
            let rising = bayes_rising(
                u64::from(n / 2),
                u64::from(n),
                u64::from(n / 2),
                u64::from(n),
                &integer_priors(&model).unwrap(),
            )
            .to_f64();
            assert_eq!(bf.to_bits(), rising.to_bits(), "n={n}");
            let lgamma = stats::bayes_factor_2x2(n / 2, n / 2, n, n);
            let rel = (bf - lgamma).abs() / bf;
            assert!(rel < 1e-12, "n={n} exact {bf} lgamma {lgamma}");
        }
    }

    #[test]
    fn integer_prior_matches_a_hand_reduced_fraction() {
        let prior = BetaPrior {
            alpha: 2.0,
            beta: 3.0,
        };
        let model = BayesFactorModel {
            alternative_group1: prior,
            alternative_group2: prior,
            null: prior,
        };
        // See the derivation in the design note: the ratio is 21/20.
        let got = bayes_factor_integer_prior(1, 0, 2, 1, &model).unwrap();
        assert_eq!(got.to_bits(), 1.05f64.to_bits());
        assert!(
            bayes_factor_integer_prior(
                1,
                0,
                1,
                1,
                &BayesFactorModel {
                    alternative_group1: BetaPrior {
                        alpha: 0.5,
                        beta: 0.5
                    },
                    alternative_group2: BetaPrior::uniform(),
                    null: BetaPrior::uniform(),
                }
            )
            .is_none()
        );
        assert!(bayes_factor_integer_prior(2, 0, 1, 1, &BayesFactorModel::uniform_v1()).is_none());
    }

    #[test]
    fn ln_bayes_tracks_the_rounded_factor() {
        let model = BayesFactorModel::uniform_v1();
        let bf = bayes_factor_integer_prior(8, 1, 12, 10, &model).unwrap();
        let ln_bf = ln_bayes_factor_integer_prior(8, 1, 12, 10, &model).unwrap();
        let rel = (ln_bf.exp() - bf).abs() / bf;
        assert!(rel < 1e-12, "{ln_bf} vs {bf}");
    }

    #[test]
    fn ln_panel_stays_near_libm_and_above_a_real_error() {
        assert_eq!(libm::log(1.0), ln_reduced::ln_reduced(1.0));
        let mut worst = 0.0f64;
        for i in 0..(32 * 64) {
            let m = 1.0 + (i as f64) / (32.0 * 64.0);
            if m >= 2.0 {
                continue;
            }
            let err = (ln_reduced::ln_reduced(m) - libm::log(m)).abs();
            worst = worst.max(err);
        }
        assert!(worst < 5.0e-16, "horner error {worst}");
        assert!(
            worst > 1.0e-18,
            "polynomial was not exercised, worst {worst}"
        );
        let _ = LN_POLY_ABS_BOUND;
    }

    #[test]
    fn chi_square_matches_f64_while_margins_are_small_and_not_after() {
        for t1 in 1..12u32 {
            for t2 in 1..12u32 {
                for n1 in 0..=t1 {
                    for n2 in 0..=t2 {
                        let left = stats::chi_squared_yates(n1, n2, t1, t2);
                        let right = chi_squared_yates_rounded(n1, n2, t1, t2);
                        if left.is_nan() {
                            assert!(right.is_nan());
                        } else {
                            assert_eq!(left.to_bits(), right.to_bits(), "{n1},{n2},{t1},{t2}");
                        }
                    }
                }
            }
        }
        // A huge but balanced cross-product can still round to the same f64.
        let tied = stats::chi_squared_yates(3_000_000_000, 1, 4_000_000_000, 4_000_000_000);
        let tied_exact = chi_squared_yates_rounded(3_000_000_000, 1, 4_000_000_000, 4_000_000_000);
        assert_eq!(tied.to_bits(), tied_exact.to_bits());
        // This unbalanced table does not: the f64 cast of |ad-bc| moves the statistic.
        let approx =
            stats::chi_squared_yates(9_045_414, 2_988_579_416, 3_273_968_006, 3_294_916_954);
        let exact =
            chi_squared_yates_rounded(9_045_414, 2_988_579_416, 3_273_968_006, 3_294_916_954);
        assert!(approx.is_finite() && exact.is_finite());
        assert_ne!(approx.to_bits(), exact.to_bits());
        assert_eq!(ulp_diff(approx, exact), 1);
    }

    #[test]
    fn crt_roundtrip_and_pairwise_rejection() {
        let moduli = [3u64, 5, 7, 11];
        let value = 100u64;
        let residues: Vec<u64> = moduli.iter().map(|m| value % m).collect();
        let back = crt_reconstruct(&residues, &moduli).unwrap();
        assert_eq!(back, vec![value]);
        assert!(crt_reconstruct(&[1, 1], &[4, 6]).is_none());
        assert!(crt_reconstruct(&[1], &[1]).is_none());
    }

    #[test]
    fn shards_are_associative_until_they_overflow() {
        let a = vec![1u64, 2, 3];
        let b = vec![4, 5, 6];
        let c = vec![7, 8, 9];
        let left = merge_count_shards(&[a.clone(), b.clone(), c.clone()]).unwrap();
        let ab = merge_count_shards(&[a.clone(), b.clone()]).unwrap();
        let right = merge_count_shards(&[ab, c.clone()]).unwrap();
        let bc = merge_count_shards(&[b, c]).unwrap();
        let other = merge_count_shards(&[a, bc]).unwrap();
        assert_eq!(left, right);
        assert_eq!(left, other);
        assert_eq!(left, vec![12, 15, 18]);
        assert!(matches!(
            merge_count_shards(&[vec![1, 2], vec![3]]),
            Err(ShardMergeError::Shape)
        ));
        assert!(matches!(
            merge_count_shards(&[vec![u64::MAX], vec![1]]),
            Err(ShardMergeError::Overflow)
        ));
        assert!(merge_count_shards(&[]).unwrap().is_empty());
    }

    #[test]
    fn line_chunks_cover_without_splitting_newlines() {
        for len in 1..10 {
            for mask in 0..(1u32 << len) {
                let mut data = vec![b'a'; len];
                for (i, byte) in data.iter_mut().enumerate() {
                    if (mask >> i) & 1 == 1 {
                        *byte = b'\n';
                    }
                }
                for target in 1..=len {
                    let chunks = partition_line_chunks(&data, target);
                    let mut cursor = 0;
                    for (index, range) in chunks.iter().enumerate() {
                        assert_eq!(range.start, cursor);
                        assert!(range.end > range.start);
                        cursor = range.end;
                        let last = data[range.end - 1];
                        let at_eof = index + 1 == chunks.len();
                        assert!(last == b'\n' || at_eof, "target {target} mask {mask}");
                    }
                    assert_eq!(cursor, data.len());
                }
            }
        }
        let data = b"a\nb\n";
        assert_eq!(partition_line_chunks(data, 1), vec![0..2, 2..4]);
    }
}
