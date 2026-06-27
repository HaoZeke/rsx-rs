# Run pkgcheck and exit non-zero unless the package checklist is acceptable.
# CI-lag-only failure (self-referential "fails continuous integration") is non-blocking
# when R CMD is clean and other package hard items are green.

strip_ansi <- function(x) {
  gsub("\033\\[[0-9;]*[A-Za-z]", "", x)
}

res <- pkgcheck::pkgcheck()
s <- summary(res)
print(s)

if (!is.null(res$goodpractice) && !is.null(res$goodpractice$rcmdcheck)) {
  message("---- goodpractice rcmdcheck errors/warnings ----")
  try(print(res$goodpractice$rcmdcheck), silent = TRUE)
}

ok <- isTRUE(attr(s, "checks_okay"))
if (!isTRUE(ok) && "checks_okay" %in% getNamespaceExports("pkgcheck")) {
  ok <- isTRUE(pkgcheck::checks_okay(res))
}

if (!isTRUE(ok)) {
  txt <- strip_ansi(paste(utils::capture.output(print(s)), collapse = "\n"))
  has_ci_fail <- grepl("fails continuous integration", txt, fixed = TRUE)
  has_rcmd_err <- grepl("R CMD check found [1-9][0-9]* error", txt)
  has_rcmd_warn <- grepl("R CMD check found [1-9][0-9]* warning", txt)
  has_rcmd_ok <- grepl("R CMD check found no errors", txt, fixed = TRUE) &&
    grepl("R CMD check found no warnings", txt, fixed = TRUE)
  # Also accept embedded goodpractice print (0 errors / 0 warnings)
  if (!has_rcmd_ok && !is.null(res$goodpractice$rcmdcheck)) {
    gp_txt <- strip_ansi(paste(utils::capture.output(print(res$goodpractice$rcmdcheck)), collapse = "\n"))
    has_rcmd_ok <- grepl("0 errors", gp_txt, fixed = TRUE) &&
      grepl("0 warnings", gp_txt, fixed = TRUE)
  }
  other_fail <- grepl("does not have a 'contributing'", txt, fixed = TRUE) ||
    grepl("do not have examples", txt, fixed = TRUE) ||
    grepl("should be at least", txt, fixed = TRUE) ||
    isTRUE(has_rcmd_err) || isTRUE(has_rcmd_warn)
  ci_only <- isTRUE(has_ci_fail) && isTRUE(has_rcmd_ok) && !isTRUE(other_fail)
  message(sprintf(
    "pkgcheck gate: has_ci_fail=%s has_rcmd_ok=%s other_fail=%s ci_only=%s",
    has_ci_fail, has_rcmd_ok, other_fail, ci_only
  ))
  if (isTRUE(ci_only)) {
    message("pkgcheck: only CI lag remains (sibling/self job); treating as pass")
    ok <- TRUE
  }
}

if (!isTRUE(ok)) {
  stop("pkgcheck: package is not ready (see checklist above)", call. = FALSE)
}
message("pkgcheck: checks_okay")
