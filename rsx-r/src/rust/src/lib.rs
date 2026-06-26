// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers
//
//! Static-library anchor for the rsxr R package.
//!
//! The package's C glue (`src/rsxr.c`) calls the rsx C API directly. Those
//! `#[no_mangle] extern "C"` symbols live in `rsx-core`; when `rsx-core` is
//! consumed as an rlib dependency and this crate is archived into a
//! `staticlib`, the linker would garbage-collect the C-API object files
//! because nothing in Rust references them. `rsxr_anchor` names every entry
//! point, which makes those object files reachable and keeps the symbols in
//! `librsxr.a` so the R package resolves them at load time.

use rsx_core::c_api::{commands, types};
use rsx_core::status;

/// References each rsx-core C entry point so its object file is archived into
/// `librsxr.a`. Exported (`#[no_mangle]`) so it is never dead-code-eliminated.
#[no_mangle]
pub extern "C" fn rsxr_anchor() -> usize {
    let entries: [*const (); 12] = [
        commands::rsx_process as *const (),
        commands::rsx_freq as *const (),
        commands::rsx_distrib as *const (),
        commands::rsx_signif as *const (),
        commands::rsx_triage as *const (),
        commands::rsx_depth as *const (),
        commands::rsx_merge as *const (),
        commands::rsx_pca as *const (),
        types::rsx_popmap_load as *const (),
        types::rsx_popmap_free as *const (),
        types::rsx_popmap_n_individuals as *const (),
        status::rsx_last_error as *const (),
    ];
    entries.len()
}
