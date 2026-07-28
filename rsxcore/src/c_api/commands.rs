// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! C-compatible wrappers for RADSex commands.

use crate::status::{catch_unwind, rsx_status_t, set_last_error};
use crate::test_method::{CorrectionMethod, TestMethod};
use std::ffi::CStr;
use std::os::raw::c_char;

/// Parse a correction method from a C string, recording the error on failure.
unsafe fn parse_correction(ptr: *const c_char) -> Result<CorrectionMethod, rsx_status_t> {
    let s = unsafe { cstr_to_string(ptr, "correction") }?;
    CorrectionMethod::parse_str(&s).map_err(|e| {
        set_last_error(&e);
        rsx_status_t::RSX_INVALID_PARAMETER
    })
}

/// Parse a statistical test from a C string, recording the error on failure.
unsafe fn parse_test(ptr: *const c_char) -> Result<TestMethod, rsx_status_t> {
    let s = unsafe { cstr_to_string(ptr, "test") }?;
    TestMethod::parse_str(&s).map_err(|e| {
        set_last_error(&e);
        rsx_status_t::RSX_INVALID_PARAMETER
    })
}

/// Helper to convert a C string pointer to a Rust string, returning error on null
/// or invalid UTF-8 (never invents an empty string for bad encodings).
pub(crate) unsafe fn cstr_to_string(
    ptr: *const c_char,
    name: &str,
) -> Result<String, rsx_status_t> {
    if ptr.is_null() {
        set_last_error(&format!("null pointer for {name}"));
        return Err(rsx_status_t::RSX_INVALID_PARAMETER);
    }
    match unsafe { CStr::from_ptr(ptr) }.to_str() {
        Ok(s) => Ok(s.to_string()),
        Err(_) => {
            set_last_error(&format!("invalid UTF-8 in {name}"));
            Err(rsx_status_t::RSX_INVALID_PARAMETER)
        }
    }
}

/// Run the `process` command.
///
/// # Safety
/// All string pointers must be valid null-terminated C strings.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsx_process(
    input_dir: *const c_char,
    output_path: *const c_char,
    n_threads: u32,
    min_depth: u32,
) -> rsx_status_t {
    catch_unwind(|| {
        let input_dir = match unsafe { cstr_to_string(input_dir, "input_dir") } {
            Ok(s) => s,
            Err(e) => return e,
        };
        let output_path = match unsafe { cstr_to_string(output_path, "output_path") } {
            Ok(s) => s,
            Err(e) => return e,
        };

        let params = crate::commands::process::ProcessParams {
            input_dir_path: input_dir,
            output_file_path: output_path,
            n_threads,
            min_depth: min_depth as u16,
            kmer_dedup: None,
        };

        match crate::commands::process::run(&params) {
            Ok(()) => rsx_status_t::RSX_SUCCESS,
            Err(e) => {
                set_last_error(&format!("process failed: {e}"));
                rsx_status_t::RSX_INTERNAL_ERROR
            }
        }
    })
}

/// Run the `freq` command.
///
/// # Safety
/// All string pointers must be valid null-terminated C strings.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsx_freq(
    table_path: *const c_char,
    output_path: *const c_char,
    min_depth: u32,
) -> rsx_status_t {
    catch_unwind(|| {
        let table_path = match unsafe { cstr_to_string(table_path, "table_path") } {
            Ok(s) => s,
            Err(e) => return e,
        };
        let output_path = match unsafe { cstr_to_string(output_path, "output_path") } {
            Ok(s) => s,
            Err(e) => return e,
        };

        let params = crate::commands::freq::FreqParams {
            markers_table_path: table_path,
            output_file_path: output_path,
            min_depth: min_depth as u16,
        };

        match crate::commands::freq::run(&params) {
            Ok(()) => rsx_status_t::RSX_SUCCESS,
            Err(e) => {
                set_last_error(&format!("freq failed: {e}"));
                rsx_status_t::RSX_INTERNAL_ERROR
            }
        }
    })
}

/// Run the `distrib` command.
///
/// # Safety
/// All string pointers must be valid null-terminated C strings.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn rsx_distrib(
    table_path: *const c_char,
    popmap_path: *const c_char,
    output_path: *const c_char,
    min_depth: u32,
    signif_threshold: f32,
    group1: *const c_char,
    group2: *const c_char,
    correction: *const c_char,
    test: *const c_char,
    output_bayes: bool,
) -> rsx_status_t {
    catch_unwind(|| {
        let table_path = match unsafe { cstr_to_string(table_path, "table_path") } {
            Ok(s) => s,
            Err(e) => return e,
        };
        let popmap_path = match unsafe { cstr_to_string(popmap_path, "popmap_path") } {
            Ok(s) => s,
            Err(e) => return e,
        };
        let output_path = match unsafe { cstr_to_string(output_path, "output_path") } {
            Ok(s) => s,
            Err(e) => return e,
        };
        let group1 = match unsafe { cstr_to_string(group1, "group1") } {
            Ok(s) => s,
            Err(e) => return e,
        };
        let group2 = match unsafe { cstr_to_string(group2, "group2") } {
            Ok(s) => s,
            Err(e) => return e,
        };
        let correction = match unsafe { parse_correction(correction) } {
            Ok(c) => c,
            Err(e) => return e,
        };
        let test_method = match unsafe { parse_test(test) } {
            Ok(t) => t,
            Err(e) => return e,
        };

        let params = crate::commands::distrib::DistribParams {
            markers_table_path: table_path,
            popmap_file_path: popmap_path,
            output_file_path: output_path,
            min_depth: min_depth as u16,
            signif_threshold,
            correction,
            test_method,
            output_bayes,
            group1,
            group2,
        };

        match crate::commands::distrib::run(&params) {
            Ok(()) => rsx_status_t::RSX_SUCCESS,
            Err(e) => {
                set_last_error(&format!("distrib failed: {e}"));
                rsx_status_t::RSX_INTERNAL_ERROR
            }
        }
    })
}

/// Run the `signif` command.
///
/// # Safety
/// All string pointers must be valid null-terminated C strings.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn rsx_signif(
    table_path: *const c_char,
    popmap_path: *const c_char,
    output_path: *const c_char,
    min_depth: u32,
    signif_threshold: f32,
    group1: *const c_char,
    group2: *const c_char,
    correction: *const c_char,
    test: *const c_char,
    output_fasta: bool,
    output_bayes: bool,
) -> rsx_status_t {
    catch_unwind(|| {
        let table_path = match unsafe { cstr_to_string(table_path, "table_path") } {
            Ok(s) => s,
            Err(e) => return e,
        };
        let popmap_path = match unsafe { cstr_to_string(popmap_path, "popmap_path") } {
            Ok(s) => s,
            Err(e) => return e,
        };
        let output_path = match unsafe { cstr_to_string(output_path, "output_path") } {
            Ok(s) => s,
            Err(e) => return e,
        };
        let group1 = match unsafe { cstr_to_string(group1, "group1") } {
            Ok(s) => s,
            Err(e) => return e,
        };
        let group2 = match unsafe { cstr_to_string(group2, "group2") } {
            Ok(s) => s,
            Err(e) => return e,
        };
        let correction = match unsafe { parse_correction(correction) } {
            Ok(c) => c,
            Err(e) => return e,
        };
        let test_method = match unsafe { parse_test(test) } {
            Ok(t) => t,
            Err(e) => return e,
        };

        let params = crate::commands::signif::SignifParams {
            markers_table_path: table_path,
            popmap_file_path: popmap_path,
            output_file_path: output_path,
            min_depth: min_depth as u16,
            signif_threshold,
            correction,
            test_method,
            output_fasta,
            output_bayes,
            group1,
            group2,
        };

        match crate::commands::signif::run(&params) {
            Ok(()) => rsx_status_t::RSX_SUCCESS,
            Err(e) => {
                set_last_error(&format!("signif failed: {e}"));
                rsx_status_t::RSX_INTERNAL_ERROR
            }
        }
    })
}

/// Run the `triage` command.
///
/// # Safety
/// All string pointers must be valid null-terminated C strings.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn rsx_triage(
    table_path: *const c_char,
    popmap_path: *const c_char,
    output_path: *const c_char,
    min_depth: u32,
    signif_threshold: f32,
    posterior_threshold: f64,
    bayes_factor_threshold: f64,
    prior_probability: f64,
    linked_probability: f64,
    group1: *const c_char,
    group2: *const c_char,
) -> rsx_status_t {
    catch_unwind(|| {
        let table_path = match unsafe { cstr_to_string(table_path, "table_path") } {
            Ok(s) => s,
            Err(e) => return e,
        };
        let popmap_path = match unsafe { cstr_to_string(popmap_path, "popmap_path") } {
            Ok(s) => s,
            Err(e) => return e,
        };
        let output_path = match unsafe { cstr_to_string(output_path, "output_path") } {
            Ok(s) => s,
            Err(e) => return e,
        };
        let group1 = match unsafe { cstr_to_string(group1, "group1") } {
            Ok(s) => s,
            Err(e) => return e,
        };
        let group2 = match unsafe { cstr_to_string(group2, "group2") } {
            Ok(s) => s,
            Err(e) => return e,
        };

        let params = crate::commands::triage::TriageParams {
            markers_table_path: table_path,
            popmap_file_path: popmap_path,
            output_file_path: output_path,
            min_depth: min_depth as u16,
            signif_threshold,
            posterior_threshold,
            bayes_factor_threshold,
            prior_probability,
            linked_probability,
            group1,
            group2,
        };

        match crate::commands::triage::run(&params) {
            Ok(()) => rsx_status_t::RSX_SUCCESS,
            Err(e) => {
                set_last_error(&format!("triage failed: {e}"));
                rsx_status_t::RSX_INTERNAL_ERROR
            }
        }
    })
}

/// Run the `depth` command.
///
/// # Safety
/// All string pointers must be valid null-terminated C strings.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsx_depth(
    table_path: *const c_char,
    popmap_path: *const c_char,
    output_path: *const c_char,
    min_frequency: f32,
    streaming: bool,
) -> rsx_status_t {
    catch_unwind(|| {
        let table_path = match unsafe { cstr_to_string(table_path, "table_path") } {
            Ok(s) => s,
            Err(e) => return e,
        };
        let popmap_path = match unsafe { cstr_to_string(popmap_path, "popmap_path") } {
            Ok(s) => s,
            Err(e) => return e,
        };
        let output_path = match unsafe { cstr_to_string(output_path, "output_path") } {
            Ok(s) => s,
            Err(e) => return e,
        };

        let params = crate::commands::depth::DepthParams {
            markers_table_path: table_path,
            popmap_file_path: popmap_path,
            output_file_path: output_path,
            min_frequency,
            streaming,
        };

        match crate::commands::depth::run(&params) {
            Ok(()) => rsx_status_t::RSX_SUCCESS,
            Err(e) => {
                set_last_error(&format!("depth failed: {e}"));
                rsx_status_t::RSX_INTERNAL_ERROR
            }
        }
    })
}

/// Run the `merge` command.
///
/// `input_files` is an array of `n_files` null-terminated C strings.
/// `buffer_size` of 0 selects the default buffer.
///
/// # Safety
/// `input_files` must point to `n_files` valid C string pointers and the
/// remaining pointers must be valid null-terminated C strings.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsx_merge(
    input_files: *const *const c_char,
    n_files: usize,
    output_path: *const c_char,
    buffer_size: usize,
    output_parquet: bool,
) -> rsx_status_t {
    catch_unwind(|| {
        if input_files.is_null() {
            set_last_error("null pointer for input_files");
            return rsx_status_t::RSX_INVALID_PARAMETER;
        }
        let mut files = Vec::with_capacity(n_files);
        for i in 0..n_files {
            let ptr = unsafe { *input_files.add(i) };
            match unsafe { cstr_to_string(ptr, "input_file") } {
                Ok(s) => files.push(s),
                Err(e) => return e,
            }
        }
        let output_path = match unsafe { cstr_to_string(output_path, "output_path") } {
            Ok(s) => s,
            Err(e) => return e,
        };

        let params = crate::commands::merge::MergeParams {
            input_files: files,
            output_file_path: output_path,
            buffer_size: if buffer_size == 0 {
                None
            } else {
                Some(buffer_size)
            },
            output_parquet,
        };

        match crate::commands::merge::run(&params) {
            Ok(()) => rsx_status_t::RSX_SUCCESS,
            Err(e) => {
                set_last_error(&format!("merge failed: {e}"));
                rsx_status_t::RSX_INTERNAL_ERROR
            }
        }
    })
}

/// Run the `pca` command.
///
/// `n_components` of 0 selects the default number of components.
///
/// # Safety
/// All string pointers must be valid null-terminated C strings.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsx_pca(
    table_path: *const c_char,
    output_dir: *const c_char,
    min_depth: u32,
    n_components: usize,
) -> rsx_status_t {
    catch_unwind(|| {
        let table_path = match unsafe { cstr_to_string(table_path, "table_path") } {
            Ok(s) => s,
            Err(e) => return e,
        };
        let output_dir = match unsafe { cstr_to_string(output_dir, "output_dir") } {
            Ok(s) => s,
            Err(e) => return e,
        };

        let params = crate::commands::pca::PcaParams {
            markers_table_path: table_path,
            output_dir,
            min_depth: min_depth as u16,
            n_components: if n_components == 0 {
                None
            } else {
                Some(n_components)
            },
        };

        match crate::commands::pca::run(&params) {
            Ok(()) => rsx_status_t::RSX_SUCCESS,
            Err(e) => {
                set_last_error(&format!("pca failed: {e}"));
                rsx_status_t::RSX_INTERNAL_ERROR
            }
        }
    })
}

#[cfg(test)]
mod cstr_tests {
    use super::*;
    use std::ffi::CString;

    #[test]
    fn cstr_to_string_accepts_valid_utf8() {
        let s = CString::new("markers.tsv").unwrap();
        let got = unsafe { cstr_to_string(s.as_ptr(), "path") }.unwrap();
        assert_eq!(got, "markers.tsv");
    }

    #[test]
    fn cstr_to_string_rejects_null() {
        let err = unsafe { cstr_to_string(std::ptr::null(), "path") }.unwrap_err();
        assert_eq!(err, rsx_status_t::RSX_INVALID_PARAMETER);
    }

    #[test]
    fn cstr_to_string_rejects_invalid_utf8() {
        // C string with invalid UTF-8 (0xFF byte) and trailing NUL.
        let bytes = [b'a', 0xFFu8, b'b', 0];
        let ptr = bytes.as_ptr() as *const c_char;
        let err = unsafe { cstr_to_string(ptr, "path") }.unwrap_err();
        assert_eq!(err, rsx_status_t::RSX_INVALID_PARAMETER);
    }
}
