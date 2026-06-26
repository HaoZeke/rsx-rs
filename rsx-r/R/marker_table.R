# Ergonomic high-level API. A marker_table wraps a marker depth table on disk;
# the analysis verbs run the corresponding rsx command into a temporary file and
# read the result back as a tibble.

#' Read an rsx core TSV output
#'
#' rsx commands prefix their TSV output with a `#Number of markers` comment
#' line; this reads past it and returns a tibble.
#' @param path TSV path.
#' @return A tibble.
#' @keywords internal
rsxr_read_core_tsv <- function(path) {
  tibble::as_tibble(
    readr::read_tsv(path, comment = "#", show_col_types = FALSE,
                    progress = FALSE)
  )
}

#' Create a marker table handle
#'
#' @param x A path to an rsx marker depth table (TSV), or a data frame in the
#'   rsx marker-table layout.
#' @return An object of class `marker_table`.
#' @export
#' @examples
#' tmp <- tempfile(fileext = ".tsv")
#' writeLines(c("#Number of markers: 0", "id\tsequence"), tmp)
#' mt <- marker_table(tmp)
#' print(mt)
marker_table <- function(x) {
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) {
      stop("marker_table: file does not exist: ", x, call. = FALSE)
    }
    return(structure(list(path = x), class = "marker_table"))
  }
  if (is.data.frame(x)) {
    tmp <- tempfile(fileext = ".tsv")
    readr::write_tsv(x, tmp)
    return(structure(list(path = tmp), class = "marker_table"))
  }
  stop("marker_table: 'x' must be a file path or a data frame", call. = FALSE)
}

#' @export
print.marker_table <- function(x, ...) {
  cat("<marker_table>\n  path:", x$path, "\n")
  invisible(x)
}

#' @export
format.marker_table <- function(x, ...) {
  paste0("<marker_table: ", x$path, ">")
}

# Verb generics (named to avoid clashing with base::signif and base::merge).

#' Bayesian sex-linkage triage
#' @param x A [marker_table].
#' @param popmap Population map path.
#' @param ... Passed to [rsx_triage()].
#' @return A tibble of triaged markers.
#' @export
#' @examples
#' \dontrun{
#' mt <- marker_table("markers.tsv")
#' triage(mt, popmap = "popmap.tsv", min_depth = 10L)
#' }
triage <- function(x, ...) UseMethod("triage")

#' @rdname triage
#' @export
triage.marker_table <- function(x, popmap, ...) {
  out <- tempfile(fileext = ".tsv")
  rsx_triage(x$path, popmap, out, ...)
  rsxr_read_core_tsv(out)
}

#' Significant sex-linked markers
#' @param x A [marker_table].
#' @param popmap Population map path.
#' @param ... Passed to [rsx_signif()].
#' @return A tibble of significant markers.
#' @export
#' @examples
#' \dontrun{
#' mt <- marker_table("markers.tsv")
#' signif_markers(mt, popmap = "popmap.tsv", test = "fisher")
#' }
signif_markers <- function(x, ...) UseMethod("signif_markers")

#' @rdname signif_markers
#' @export
signif_markers.marker_table <- function(x, popmap, ...) {
  out <- tempfile(fileext = ".tsv")
  rsx_signif(x$path, popmap, out, ...)
  rsxr_read_core_tsv(out)
}

#' Marker distribution across groups
#' @param x A [marker_table].
#' @param popmap Population map path.
#' @param ... Passed to [rsx_distrib()].
#' @return A tibble of the marker distribution.
#' @export
#' @examples
#' \dontrun{
#' mt <- marker_table("markers.tsv")
#' distrib(mt, popmap = "popmap.tsv", group1 = "M", group2 = "F")
#' }
distrib <- function(x, ...) UseMethod("distrib")

#' @rdname distrib
#' @export
distrib.marker_table <- function(x, popmap, ...) {
  out <- tempfile(fileext = ".tsv")
  rsx_distrib(x$path, popmap, out, ...)
  rsxr_read_core_tsv(out)
}

#' Per-individual marker depth
#' @param x A [marker_table].
#' @param popmap Population map path.
#' @param ... Passed to [rsx_depth()].
#' @return A tibble of per-individual depths.
#' @export
#' @examples
#' \dontrun{
#' mt <- marker_table("markers.tsv")
#' depth(mt, popmap = "popmap.tsv")
#' }
depth <- function(x, ...) UseMethod("depth")

#' @rdname depth
#' @export
depth.marker_table <- function(x, popmap, ...) {
  out <- tempfile(fileext = ".tsv")
  rsx_depth(x$path, popmap, out, ...)
  rsxr_read_core_tsv(out)
}

#' Per-marker allele frequencies
#' @param x A [marker_table].
#' @param ... Passed to [rsx_freq()].
#' @return A tibble of frequencies.
#' @export
#' @examples
#' \dontrun{
#' mt <- marker_table("markers.tsv")
#' frequencies(mt, min_depth = 5L)
#' }
frequencies <- function(x, ...) UseMethod("frequencies")

#' @rdname frequencies
#' @export
frequencies.marker_table <- function(x, ...) {
  out <- tempfile(fileext = ".tsv")
  rsx_freq(x$path, out, ...)
  rsxr_read_core_tsv(out)
}
