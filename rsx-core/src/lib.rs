// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! Core library for RADSex: sex-determination analysis from RAD-Sequencing data.
//!
//! This crate provides the computational core for analyzing RAD-seq data to
//! identify sex-linked markers. It exposes both a Rust API and a C-compatible
//! FFI layer (via cbindgen) for integration with C++, R, and other languages.

pub mod bitset;
pub mod c_api;
pub mod kmer;
pub mod commands;
pub mod io;
pub mod marker;
pub mod markers_table;
pub mod popmap;
pub mod stats;
pub mod status;
pub mod test_method;
