#' rsxr: R bindings for the rsx RAD-seq sex-determination toolkit
#'
#' rsxr binds the rsx C API directly through R's native C interface. The
#' low-level functions ([rsx_process()], [rsx_freq()], [rsx_distrib()],
#' [rsx_signif()], [rsx_triage()], [rsx_depth()], [rsx_merge()],
#' [rsx_pca()]) mirror the rsx CLI one to one. The high-level
#' [marker_table()] workflow returns tibbles.
#'
#' @keywords internal
#' @useDynLib rsxr, .registration = TRUE
"_PACKAGE"
