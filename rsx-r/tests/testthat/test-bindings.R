test_that("rsx_version reports a semantic version", {
  v <- rsx_version()
  expect_type(v, "character")
  expect_match(v, "^[0-9]+\\.[0-9]+\\.[0-9]+")
})

test_that("marker_table validates its input", {
  expect_error(marker_table("/does/not/exist.tsv"), "does not exist")
  expect_error(marker_table(42), "file path or a data frame")
})

test_that("a missing input surfaces the rsx error message", {
  # The C API returns a non-zero status; the glue turns it into an R error
  # carrying the thread-local rsx_last_error() text.
  expect_error(
    rsx_freq("/no/such/markers.tsv", tempfile(), min_depth = 1L),
    "rsx:"
  )
})

test_that("marker_table prints its path", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines("#Number of markers: 0\nid\tsequence", tmp)
  mt <- marker_table(tmp)
  expect_s3_class(mt, "marker_table")
  expect_output(print(mt), "marker_table")
})
