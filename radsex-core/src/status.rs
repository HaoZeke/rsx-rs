// GPL-3.0-or-later
// Copyright 2024--present radsex-rs developers

//! Error handling following the metatensor/rgpot pattern.
//!
//! Provides:
//! 1. [`radsex_status_t`] -- integer enum returned from every `extern "C"` function.
//! 2. Thread-local error message retrievable via [`radsex_last_error()`].
//! 3. [`catch_unwind`] wrapper to prevent panics crossing the FFI boundary.

use std::cell::RefCell;
use std::ffi::CString;
use std::os::raw::c_char;

/// Status codes returned by all C API functions.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(non_camel_case_types)]
pub enum radsex_status_t {
    /// Operation completed successfully.
    RADSEX_SUCCESS = 0,
    /// An invalid parameter was passed (null pointer, wrong size, etc.).
    RADSEX_INVALID_PARAMETER = 1,
    /// An internal error occurred (e.g. a Rust panic was caught).
    RADSEX_INTERNAL_ERROR = 2,
    /// An I/O error occurred (file not found, permission denied, etc.).
    RADSEX_IO_ERROR = 3,
    /// An alignment error occurred (index missing, alignment failed, etc.).
    RADSEX_ALIGNMENT_ERROR = 4,
}

thread_local! {
    static LAST_ERROR: RefCell<CString> = RefCell::new(CString::default());
}

/// Store an error message in the thread-local slot.
pub fn set_last_error(msg: &str) {
    LAST_ERROR.with(|cell| {
        let c = CString::new(msg).unwrap_or_else(|_| {
            CString::new("(error message contained interior NUL)").unwrap()
        });
        *cell.borrow_mut() = c;
    });
}

/// Retrieve a pointer to the last error message for the current thread.
///
/// The pointer is valid until the next call to any `radsex_*` function
/// on the same thread.
///
/// # Safety
/// This is intended to be called from C. The returned pointer must not
/// be freed by the caller.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn radsex_last_error() -> *const c_char {
    LAST_ERROR.with(|cell| cell.borrow().as_ptr())
}

/// Execute a closure, catching any panics and converting them to status codes.
pub fn catch_unwind<F>(f: F) -> radsex_status_t
where
    F: FnOnce() -> radsex_status_t + std::panic::UnwindSafe,
{
    match std::panic::catch_unwind(f) {
        Ok(status) => status,
        Err(e) => {
            let msg = if let Some(s) = e.downcast_ref::<&str>() {
                s.to_string()
            } else if let Some(s) = e.downcast_ref::<String>() {
                s.clone()
            } else {
                "unknown panic".to_string()
            };
            set_last_error(&msg);
            radsex_status_t::RADSEX_INTERNAL_ERROR
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_set_and_get_last_error() {
        set_last_error("test error");
        let ptr = unsafe { radsex_last_error() };
        let msg = unsafe { std::ffi::CStr::from_ptr(ptr) };
        assert_eq!(msg.to_str().unwrap(), "test error");
    }

    #[test]
    fn test_catch_unwind_success() {
        let status = catch_unwind(|| radsex_status_t::RADSEX_SUCCESS);
        assert_eq!(status, radsex_status_t::RADSEX_SUCCESS);
    }

    #[test]
    fn test_catch_unwind_panic() {
        let status = catch_unwind(|| panic!("boom"));
        assert_eq!(status, radsex_status_t::RADSEX_INTERNAL_ERROR);
        let ptr = unsafe { radsex_last_error() };
        let msg = unsafe { std::ffi::CStr::from_ptr(ptr) };
        assert_eq!(msg.to_str().unwrap(), "boom");
    }

    #[test]
    fn test_error_overwrite() {
        set_last_error("first");
        set_last_error("second");
        let ptr = unsafe { radsex_last_error() };
        let msg = unsafe { std::ffi::CStr::from_ptr(ptr) };
        assert_eq!(msg.to_str().unwrap(), "second");
    }

    #[test]
    fn test_interior_nul_is_handled() {
        set_last_error("has\0interior nul");
        let ptr = unsafe { radsex_last_error() };
        let msg = unsafe { std::ffi::CStr::from_ptr(ptr) };
        assert_eq!(
            msg.to_str().unwrap(),
            "(error message contained interior NUL)"
        );
    }

    #[test]
    fn test_catch_unwind_returns_callback_status() {
        let status = catch_unwind(|| radsex_status_t::RADSEX_IO_ERROR);
        assert_eq!(status, radsex_status_t::RADSEX_IO_ERROR);
    }
}
