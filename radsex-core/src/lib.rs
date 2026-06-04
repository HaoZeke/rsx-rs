// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

// Allow unwrap/expect in tests and internal trusted code (e.g. wire format deserializers,
// test fixtures). Production paths should prefer ? or explicit error mapping.
#![cfg_attr(test, allow(clippy::unwrap_used, clippy::expect_used))]

//! Core library for RADSex: sex-determination analysis from RAD-Sequencing data.
//!
//! This crate provides the computational core for analyzing RAD-seq data to
//! identify sex-linked markers. It exposes both a Rust API and a C-compatible
//! FFI layer (via cbindgen) for integration with C++, R, and other languages.
//!
//! # Key types
//!
//! - [`markers_table::MarkersTableStream`]: streaming / bounded-memory reader over
//!   marker depth tables (the central data structure for all commands).
//! - [`commands`]: implementations of `process`, `distrib`, `signif`, `pca`, `merge`, etc.
//! - [`stats`]: Bayesian and frequentist tests (Bayes factor, posterior, chi2, Fisher, G-test).
//!
//! All commands are designed to run with O(n_individuals) or bounded temporary
//! memory even on 50 GB+ tables.

pub mod bitset;
pub mod c_api;
pub mod commands;
pub mod io;
pub mod kmer;
pub mod marker;
pub mod markers_table;
pub mod popmap;
pub mod source;
pub mod stats;
pub mod status;
pub mod test_method;
