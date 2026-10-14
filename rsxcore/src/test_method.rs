// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Statistical test and correction method selection.

/// Which statistical test to use for sex-marker association.
#[derive(Debug, Clone, Copy, Default)]
pub enum TestMethod {
    #[default]
    ChiSquared,
    Fisher,
    GTest,
}

/// Which multiple testing correction to apply.
#[derive(Debug, Clone, Copy, Default)]
pub enum CorrectionMethod {
    #[default]
    Bonferroni,
    Fdr,
    None,
}

impl TestMethod {
    pub fn parse_str(s: &str) -> Result<Self, String> {
        match s.to_lowercase().as_str() {
            "chisq" | "chi-squared" | "chi2" => Ok(TestMethod::ChiSquared),
            "fisher" | "fisher-exact" => Ok(TestMethod::Fisher),
            "gtest" | "g-test" | "g" => Ok(TestMethod::GTest),
            _ => Err(format!(
                "Unknown test method: {}. Options: chisq, fisher, gtest",
                s
            )),
        }
    }
}

impl CorrectionMethod {
    pub fn parse_str(s: &str) -> Result<Self, String> {
        match s.to_lowercase().as_str() {
            "bonferroni" | "bonf" => Ok(CorrectionMethod::Bonferroni),
            "fdr" | "bh" | "benjamini-hochberg" => Ok(CorrectionMethod::Fdr),
            "none" | "disabled" => Ok(CorrectionMethod::None),
            _ => Err(format!(
                "Unknown correction: {}. Options: bonferroni, fdr, none",
                s
            )),
        }
    }
}

/// Compute p-value using the selected test method.
pub fn compute_p(method: TestMethod, n_g1: u32, n_g2: u32, total_g1: u32, total_g2: u32) -> f64 {
    match method {
        TestMethod::ChiSquared => crate::stats::p_association(n_g1, n_g2, total_g1, total_g2),
        TestMethod::Fisher => crate::stats::fisher_exact(n_g1, n_g2, total_g1, total_g2),
        TestMethod::GTest => crate::stats::g_test(n_g1, n_g2, total_g1, total_g2),
    }
}

/// Precomputed association p-values for every `(n_g1, n_g2)` pair.
///
/// A marker's contingency table is fully determined by the two group counts.
/// With group sizes in the tens or hundreds there are far fewer distinct tables
/// than markers, so one lattice replaces a per-marker Fisher or chi-squared
/// evaluation. Values are exactly `compute_p` for the same inputs.
pub struct AssociationLattice {
    total_g1: u32,
    total_g2: u32,
    method: TestMethod,
    stride: usize,
    values: Vec<f64>,
}

impl AssociationLattice {
    /// Build the lattice when it is smaller than about a quarter million cells.
    /// Larger group counts fall back to `compute_p` inside [`Self::p`].
    pub fn build(method: TestMethod, total_g1: u32, total_g2: u32) -> Self {
        let rows = total_g1 as usize + 1;
        let cols = total_g2 as usize + 1;
        let cells = rows.saturating_mul(cols);
        if cells > 250_000 {
            return Self {
                total_g1,
                total_g2,
                method,
                stride: 0,
                values: Vec::new(),
            };
        }
        let mut values = Vec::with_capacity(cells);
        for g1 in 0..=total_g1 {
            for g2 in 0..=total_g2 {
                values.push(compute_p(method, g1, g2, total_g1, total_g2));
            }
        }
        Self {
            total_g1,
            total_g2,
            method,
            stride: cols,
            values,
        }
    }

    #[inline]
    pub fn p(&self, n_g1: u32, n_g2: u32) -> f64 {
        if self.stride == 0 || n_g1 > self.total_g1 || n_g2 > self.total_g2 {
            return compute_p(self.method, n_g1, n_g2, self.total_g1, self.total_g2);
        }
        self.values[n_g1 as usize * self.stride + n_g2 as usize]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lattice_matches_direct_p_values() {
        for method in [
            TestMethod::ChiSquared,
            TestMethod::Fisher,
            TestMethod::GTest,
        ] {
            let lattice = AssociationLattice::build(method, 6, 5);
            for g1 in 0..=6 {
                for g2 in 0..=5 {
                    let direct = compute_p(method, g1, g2, 6, 5);
                    assert_eq!(lattice.p(g1, g2).to_bits(), direct.to_bits());
                }
            }
        }
    }
}
