// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Statistical test and correction method selection.

/// Which statistical test to use for sex-marker association.
#[derive(Debug, Clone, Copy, Default)]
pub enum TestMethod {
    /// Yates-corrected chi-square on the 2x2 presence table. This is the default.
    #[default]
    ChiSquared,
    /// Two-sided Fisher exact test on the same 2x2 table.
    Fisher,
    /// G-test (likelihood ratio) on the same 2x2 table.
    GTest,
}

/// Which multiple testing correction to apply.
#[derive(Debug, Clone, Copy, Default)]
pub enum CorrectionMethod {
    /// Divide the p-value threshold by the number of markers. This is the default.
    #[default]
    Bonferroni,
    /// Benjamini-Hochberg false-discovery rate.
    Fdr,
    /// Compare the raw p-value with the threshold.
    None,
}

impl TestMethod {
    /// Parse a test name from the CLI or a run profile.
    ///
    /// Accepted values are `chisq` (aliases `chi-squared`, `chi2`), `fisher`
    /// (alias `fisher-exact`), and `gtest` (aliases `g-test`, `g`).
    pub fn parse_str(s: &str) -> Result<Self, String> {
        match s.to_lowercase().as_str() {
            "chisq" | "chi-squared" | "chi2" => Ok(TestMethod::ChiSquared),
            "fisher" | "fisher-exact" => Ok(TestMethod::Fisher),
            "gtest" | "g-test" | "g" => Ok(TestMethod::GTest),
            _ => Err(format!(
                "test '{s}' is not recognized. Accepted values: chisq (aliases: chi-squared, chi2), fisher (alias: fisher-exact), gtest (aliases: g-test, g)"
            )),
        }
    }
}

impl CorrectionMethod {
    /// Parse a correction name from the CLI or a run profile.
    ///
    /// Accepted values are `bonferroni` (alias `bonf`), `fdr` (aliases `bh`,
    /// `benjamini-hochberg`), and `none` (alias `disabled`).
    pub fn parse_str(s: &str) -> Result<Self, String> {
        match s.to_lowercase().as_str() {
            "bonferroni" | "bonf" => Ok(CorrectionMethod::Bonferroni),
            "fdr" | "bh" | "benjamini-hochberg" => Ok(CorrectionMethod::Fdr),
            "none" | "disabled" => Ok(CorrectionMethod::None),
            _ => Err(format!(
                "correction '{s}' is not recognized. Accepted values: bonferroni (alias: bonf), fdr (aliases: bh, benjamini-hochberg), none (alias: disabled)"
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
