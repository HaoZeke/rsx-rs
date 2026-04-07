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
            _ => Err(format!("Unknown test method: {}. Options: chisq, fisher, gtest", s)),
        }
    }
}

impl CorrectionMethod {
    pub fn parse_str(s: &str) -> Result<Self, String> {
        match s.to_lowercase().as_str() {
            "bonferroni" | "bonf" => Ok(CorrectionMethod::Bonferroni),
            "fdr" | "bh" | "benjamini-hochberg" => Ok(CorrectionMethod::Fdr),
            "none" | "disabled" => Ok(CorrectionMethod::None),
            _ => Err(format!("Unknown correction: {}. Options: bonferroni, fdr, none", s)),
        }
    }
}

/// Compute p-value using the selected test method.
pub fn compute_p(
    method: TestMethod,
    n_g1: u32, n_g2: u32, total_g1: u32, total_g2: u32,
) -> f64 {
    match method {
        TestMethod::ChiSquared => crate::stats::p_association(n_g1, n_g2, total_g1, total_g2),
        TestMethod::Fisher => crate::stats::fisher_exact(n_g1, n_g2, total_g1, total_g2),
        TestMethod::GTest => crate::stats::g_test(n_g1, n_g2, total_g1, total_g2),
    }
}
