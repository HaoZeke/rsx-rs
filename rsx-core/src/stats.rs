// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Statistical functions: chi-squared test with Yates correction,
//! Bonferroni multiple testing correction, and group bias.

use statrs::distribution::{ChiSquared, ContinuousCDF};
use std::fmt;

/// Format a float like C++ `operator<<` default: `%g` with 6 significant digits.
/// This matches the C++ radsex output format exactly.
pub struct Cg(pub f64);

impl fmt::Display for Cg {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // %g with 6 significant digits, strip trailing zeros
        let val = self.0;
        if val == 0.0 {
            return write!(f, "0");
        }
        let abs = val.abs();
        let exp = abs.log10().floor() as i32;
        if (-4..6).contains(&exp) {
            // Fixed notation
            let precision = (5 - exp).max(0) as usize;
            let formatted = format!("{:.*}", precision, val);
            // Strip trailing zeros after decimal point
            if formatted.contains('.') {
                let trimmed = formatted.trim_end_matches('0').trim_end_matches('.');
                write!(f, "{trimmed}")
            } else {
                write!(f, "{formatted}")
            }
        } else {
            // Scientific notation
            let mantissa = val / 10f64.powi(exp);
            let formatted = format!("{:.5}", mantissa);
            let trimmed = formatted.trim_end_matches('0').trim_end_matches('.');
            write!(f, "{trimmed}e{exp:+03}")
        }
    }
}

/// Chi-squared statistic with Yates continuity correction for a 2x2 table.
///
/// Implements the shortcut formula:
///   chi2 = N * (|ad - bc| - N/2)^2 / (a+b)(c+d)(a+c)(b+d)
///
/// where the contingency table is:
///   |           | marker present | marker absent |
///   |-----------|----------------|---------------|
///   | group1    | n_group1       | total1 - n1   |
///   | group2    | n_group2       | total2 - n2   |
pub fn chi_squared_yates(
    n_group1: u32,
    n_group2: u32,
    total_group1: u32,
    total_group2: u32,
) -> f64 {
    let n = (total_group1 + total_group2) as f64;
    let ns = total_group1 as f64;
    let nf = total_group2 as f64;
    let na = (n_group1 + n_group2) as f64;
    let nb = n - na;

    let ad_bc = ((n_group1 as i64) * (total_group2 as i64)
        - (n_group2 as i64) * (total_group1 as i64))
        .unsigned_abs() as f64;
    let yates = (ad_bc - n / 2.0).max(0.0);

    n * yates * yates / (ns * nf * na * nb)
}

/// P-value for a chi-squared statistic with df=1.
pub fn chi_squared_p(chi_sq: f64) -> f64 {
    if chi_sq.is_nan() || chi_sq <= 0.0 {
        return 1.0;
    }
    let dist = ChiSquared::new(1.0).unwrap();
    (1.0 - dist.cdf(chi_sq)).min(1.0)
}

/// Compute p-value of association with group using chi-squared test
/// with Yates correction. Matches C++ `get_p_association` exactly.
pub fn p_association(n_group1: u32, n_group2: u32, total_group1: u32, total_group2: u32) -> f64 {
    let chi_sq = chi_squared_yates(n_group1, n_group2, total_group1, total_group2);

    // NaN guard: if chi_sq is NaN (e.g. division by zero when marker
    // is present in all individuals), return p = 1.0
    let p = if chi_sq.is_nan() {
        1.0
    } else {
        chi_squared_p(chi_sq).min(1.0)
    };

    // Lower bound: 1e-16 (matching C++ behavior)
    p.max(1e-16)
}

/// Apply Bonferroni correction to a p-value.
pub fn bonferroni_correct(p: f64, n_markers: u64) -> f64 {
    (p * n_markers as f64).min(1.0)
}

/// Group bias: difference in marker frequency between two groups.
/// Ranges from -1.0 (only in group2) to +1.0 (only in group1).
pub fn group_bias(n_group1: u32, total_group1: u32, n_group2: u32, total_group2: u32) -> f64 {
    (n_group1 as f64 / total_group1 as f64) - (n_group2 as f64 / total_group2 as f64)
}

/// Find median of a mutable slice (partially sorts in-place).
pub fn find_median(data: &mut [u16]) -> u16 {
    let len = data.len();
    if len == 0 {
        return 0;
    }
    data.sort_unstable();
    if len % 2 == 0 {
        let a = data[len / 2 - 1] as u32;
        let b = data[len / 2] as u32;
        ((a + b) / 2) as u16
    } else {
        data[len / 2]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chi_squared_known_values() {
        // 10 males with marker, 0 females, out of 15 males and 10 females
        let chi_sq = chi_squared_yates(10, 0, 15, 10);
        assert!(chi_sq > 0.0);
        let p = chi_squared_p(chi_sq);
        assert!(p < 0.05);
    }

    #[test]
    fn test_equal_distribution_not_significant() {
        // 5 males, 5 females, out of 10 each -- no association
        let p = p_association(5, 5, 10, 10);
        assert!(p > 0.05);
    }

    #[test]
    fn test_p_association_floor() {
        // Extreme case: should be capped at 1e-16
        let p = p_association(20, 0, 20, 20);
        assert!(p >= 1e-16);
    }

    #[test]
    fn test_bonferroni() {
        assert_eq!(bonferroni_correct(0.01, 10), 0.1);
        assert_eq!(bonferroni_correct(0.5, 10), 1.0); // capped at 1.0
    }

    #[test]
    fn test_group_bias() {
        let bias = group_bias(10, 10, 0, 10);
        assert!((bias - 1.0).abs() < f64::EPSILON);
        let bias = group_bias(0, 10, 10, 10);
        assert!((bias - (-1.0)).abs() < f64::EPSILON);
    }

    #[test]
    fn test_find_median_odd() {
        let mut data = vec![3, 1, 2];
        assert_eq!(find_median(&mut data), 2);
    }

    #[test]
    fn test_find_median_even() {
        let mut data = vec![4, 1, 3, 2];
        assert_eq!(find_median(&mut data), 2); // (2+3)/2 = 2 (integer division)
    }
}
