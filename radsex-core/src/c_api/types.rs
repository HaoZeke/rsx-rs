// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers

//! C-compatible type constructors and destructors.

use crate::status::{catch_unwind, rsx_status_t, set_last_error};
use std::ffi::CStr;
use std::os::raw::c_char;

/// Opaque handle to a loaded popmap.
#[allow(non_camel_case_types)]
pub struct rsx_popmap_t {
    pub(crate) inner: crate::popmap::Popmap,
}

/// Load a popmap from a file path.
///
/// # Safety
/// `path` must be a valid null-terminated C string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsx_popmap_load(
    path: *const c_char,
    out: *mut *mut rsx_popmap_t,
) -> rsx_status_t {
    catch_unwind(|| {
        if path.is_null() || out.is_null() {
            set_last_error("null pointer passed to rsx_popmap_load");
            return rsx_status_t::RSX_INVALID_PARAMETER;
        }

        let path_str = unsafe { CStr::from_ptr(path) }.to_str().unwrap_or("");
        match crate::popmap::Popmap::from_file(std::path::Path::new(path_str)) {
            Ok(popmap) => {
                let boxed = Box::new(rsx_popmap_t { inner: popmap });
                unsafe { *out = Box::into_raw(boxed) };
                rsx_status_t::RSX_SUCCESS
            }
            Err(e) => {
                set_last_error(&format!("Failed to load popmap: {e}"));
                rsx_status_t::RSX_IO_ERROR
            }
        }
    })
}

/// Free a popmap handle.
///
/// # Safety
/// `popmap` must have been created by `rsx_popmap_load`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsx_popmap_free(popmap: *mut rsx_popmap_t) {
    if !popmap.is_null() {
        unsafe { drop(Box::from_raw(popmap)) };
    }
}

/// Get the number of individuals in a popmap.
///
/// # Safety
/// `popmap` must be a valid handle from `rsx_popmap_load`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsx_popmap_n_individuals(popmap: *const rsx_popmap_t) -> u16 {
    if popmap.is_null() {
        return 0;
    }
    unsafe { (*popmap).inner.n_individuals }
}
