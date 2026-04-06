// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! C-compatible wrappers for RADSex commands.

use crate::status::{rsx_status_t, catch_unwind, set_last_error};
use std::ffi::CStr;
use std::os::raw::c_char;

/// Helper to convert a C string pointer to a Rust string, returning error on null.
unsafe fn cstr_to_string(ptr: *const c_char, name: &str) -> Result<String, rsx_status_t> {
    if ptr.is_null() {
        set_last_error(&format!("null pointer for {name}"));
        return Err(rsx_status_t::RSX_INVALID_PARAMETER);
    }
    Ok(unsafe { CStr::from_ptr(ptr) }
        .to_str()
        .unwrap_or("")
        .to_string())
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
