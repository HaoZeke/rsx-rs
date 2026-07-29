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
  expect_match(format(mt), "marker_table")
})

test_that("marker_table accepts a data frame", {
  df <- data.frame(id = "m1", sequence = "ACGT", check.names = FALSE)
  mt <- marker_table(df)
  expect_s3_class(mt, "marker_table")
  expect_true(file.exists(mt$path))
})

test_that("low-level commands error on missing paths", {
  miss <- tempfile()
  out <- tempfile()
  expect_error(rsx_process(miss, out), "rsx:")
  expect_error(rsx_distrib(miss, miss, out), "rsx:")
  expect_error(rsx_signif(miss, miss, out), "rsx:")
  expect_error(rsx_triage(miss, miss, out), "rsx:")
  expect_error(rsx_depth(miss, miss, out), "rsx:")
  expect_error(rsx_merge(miss, out), "rsx:")
  expect_error(rsx_pca(miss, tempdir()), "rsx:")
})

test_that("Bayesian R bindings accept the complete directional model", {
  miss <- tempfile()
  out <- tempfile()
  model <- list(
    prior_probability = 0.02,
    linked_probability = 0.85,
    null_prevalence = 0.4,
    group1_linked_weight = 0.7,
    bf_group1_alpha = 8.0,
    bf_group1_beta = 2.0,
    bf_group2_alpha = 2.0,
    bf_group2_beta = 8.0,
    bf_null_alpha = 10.0,
    bf_null_beta = 10.0,
    posterior_linked_family = "beta",
    posterior_linked_alpha = 9.0,
    posterior_linked_beta = 1.0,
    posterior_null_family = "beta",
    posterior_null_alpha = 5.0,
    posterior_null_beta = 5.0
  )
  expect_error(do.call(rsx_distrib, c(list(miss, miss, out, output_bayes = TRUE), model)), "rsx:")
  expect_error(do.call(rsx_signif, c(list(miss, miss, out, bayes = TRUE), model)), "rsx:")
  expect_error(do.call(rsx_triage, c(list(miss, miss, out), model)), "rsx:")
})

test_that("high-level verbs error when inputs are missing", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines("#Number of markers: 0\nid\tsequence", tmp)
  mt <- marker_table(tmp)
  miss <- tempfile()
  expect_error(triage(mt, popmap = miss), "rsx:")
  expect_error(signif_markers(mt, popmap = miss), "rsx:")
  expect_error(distrib(mt, popmap = miss), "rsx:")
  expect_error(depth(mt, popmap = miss), "rsx:")
  # frequencies only needs the table path; missing file still errors from C
  bad <- marker_table(tmp)
  # corrupt by pointing at a gone path: build handle then delete
  path <- bad$path
  file.remove(path)
  # object still holds path; frequencies should fail via rsx
  expect_error(frequencies(bad), "rsx:")
})
