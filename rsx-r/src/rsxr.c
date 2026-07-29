// GPL-3.0-or-later
// Copyright 2024--present rsx-rs developers
//
// Raw C <-> R glue over the rsx C API. Each entry point converts R SEXP
// arguments to the C types expected by rsx.h, calls into the static
// librsx_core, and turns a non-zero rsx_status_t into an R error carrying
// the thread-local rsx_last_error() message.

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdint.h>
#include <stddef.h>

#include "rsx.h"

static void rsxr_check(rsx_status_t st) {
  if (st != RSX_SUCCESS) {
    const char *msg = rsx_last_error();
    Rf_error("rsx: %s (status %d)", (msg != NULL) ? msg : "unknown error", (int) st);
  }
}

static const char *rsxr_str(SEXP x, const char *name) {
  if (TYPEOF(x) != STRSXP || LENGTH(x) < 1 || STRING_ELT(x, 0) == NA_STRING) {
    Rf_error("rsx: argument '%s' must be a non-NA character scalar", name);
  }
  return CHAR(STRING_ELT(x, 0));
}

SEXP C_rsx_process(SEXP input_dir, SEXP output_file, SEXP threads, SEXP min_depth) {
  rsx_status_t st = rsx_process(
      rsxr_str(input_dir, "input_dir"),
      rsxr_str(output_file, "output_file"),
      (uint32_t) Rf_asInteger(threads),
      (uint32_t) Rf_asInteger(min_depth));
  rsxr_check(st);
  return R_NilValue;
}

SEXP C_rsx_freq(SEXP table_path, SEXP output_file, SEXP min_depth) {
  rsx_status_t st = rsx_freq(
      rsxr_str(table_path, "table_path"),
      rsxr_str(output_file, "output_file"),
      (uint32_t) Rf_asInteger(min_depth));
  rsxr_check(st);
  return R_NilValue;
}

SEXP C_rsx_distrib(SEXP table_path, SEXP popmap_path, SEXP output_file,
                   SEXP min_depth, SEXP signif_threshold, SEXP group1,
                   SEXP group2, SEXP correction, SEXP test, SEXP output_bayes,
                   SEXP prior_probability, SEXP linked_probability,
                   SEXP null_prevalence, SEXP group1_linked_weight,
                   SEXP bf_group1_alpha, SEXP bf_group1_beta,
                   SEXP bf_group2_alpha, SEXP bf_group2_beta,
                   SEXP bf_null_alpha, SEXP bf_null_beta,
                   SEXP posterior_linked_family, SEXP posterior_linked_probability,
                   SEXP posterior_linked_alpha,
                   SEXP posterior_linked_beta, SEXP posterior_null_family,
                   SEXP posterior_null_probability, SEXP posterior_null_alpha, SEXP posterior_null_beta) {
  rsx_status_t st = rsx_distrib(
      rsxr_str(table_path, "table_path"),
      rsxr_str(popmap_path, "popmap_path"),
      rsxr_str(output_file, "output_file"),
      (uint32_t) Rf_asInteger(min_depth),
      (float) Rf_asReal(signif_threshold),
      rsxr_str(group1, "group1"),
      rsxr_str(group2, "group2"),
      rsxr_str(correction, "correction"),
      rsxr_str(test, "test"),
      (bool) (Rf_asLogical(output_bayes) == TRUE),
      (double) Rf_asReal(prior_probability),
      (double) Rf_asReal(linked_probability),
      (double) Rf_asReal(null_prevalence),
      (double) Rf_asReal(group1_linked_weight),
      (double) Rf_asReal(bf_group1_alpha),
      (double) Rf_asReal(bf_group1_beta),
      (double) Rf_asReal(bf_group2_alpha),
      (double) Rf_asReal(bf_group2_beta),
      (double) Rf_asReal(bf_null_alpha),
      (double) Rf_asReal(bf_null_beta),
      rsxr_str(posterior_linked_family, "posterior_linked_family"),
      (double) Rf_asReal(posterior_linked_probability),
      (double) Rf_asReal(posterior_linked_alpha),
      (double) Rf_asReal(posterior_linked_beta),
      rsxr_str(posterior_null_family, "posterior_null_family"),
      (double) Rf_asReal(posterior_null_probability),
      (double) Rf_asReal(posterior_null_alpha),
      (double) Rf_asReal(posterior_null_beta));
  rsxr_check(st);
  return R_NilValue;
}

SEXP C_rsx_signif(SEXP table_path, SEXP popmap_path, SEXP output_file,
                  SEXP min_depth, SEXP signif_threshold, SEXP group1,
                  SEXP group2, SEXP correction, SEXP test, SEXP output_fasta,
                  SEXP output_bayes, SEXP prior_probability,
                  SEXP linked_probability, SEXP null_prevalence,
                  SEXP group1_linked_weight, SEXP bf_group1_alpha,
                  SEXP bf_group1_beta, SEXP bf_group2_alpha,
                  SEXP bf_group2_beta, SEXP bf_null_alpha,
                  SEXP bf_null_beta, SEXP posterior_linked_family,
                  SEXP posterior_linked_probability,
                  SEXP posterior_linked_alpha, SEXP posterior_linked_beta,
                  SEXP posterior_null_family, SEXP posterior_null_probability,
                  SEXP posterior_null_alpha, SEXP posterior_null_beta) {
  rsx_status_t st = rsx_signif(
      rsxr_str(table_path, "table_path"),
      rsxr_str(popmap_path, "popmap_path"),
      rsxr_str(output_file, "output_file"),
      (uint32_t) Rf_asInteger(min_depth),
      (float) Rf_asReal(signif_threshold),
      rsxr_str(group1, "group1"),
      rsxr_str(group2, "group2"),
      rsxr_str(correction, "correction"),
      rsxr_str(test, "test"),
      (bool) (Rf_asLogical(output_fasta) == TRUE),
      (bool) (Rf_asLogical(output_bayes) == TRUE),
      (double) Rf_asReal(prior_probability),
      (double) Rf_asReal(linked_probability),
      (double) Rf_asReal(null_prevalence),
      (double) Rf_asReal(group1_linked_weight),
      (double) Rf_asReal(bf_group1_alpha),
      (double) Rf_asReal(bf_group1_beta),
      (double) Rf_asReal(bf_group2_alpha),
      (double) Rf_asReal(bf_group2_beta),
      (double) Rf_asReal(bf_null_alpha),
      (double) Rf_asReal(bf_null_beta),
      rsxr_str(posterior_linked_family, "posterior_linked_family"),
      (double) Rf_asReal(posterior_linked_probability),
      (double) Rf_asReal(posterior_linked_alpha),
      (double) Rf_asReal(posterior_linked_beta),
      rsxr_str(posterior_null_family, "posterior_null_family"),
      (double) Rf_asReal(posterior_null_probability),
      (double) Rf_asReal(posterior_null_alpha),
      (double) Rf_asReal(posterior_null_beta));
  rsxr_check(st);
  return R_NilValue;
}

SEXP C_rsx_triage(SEXP table_path, SEXP popmap_path, SEXP output_file,
                  SEXP min_depth, SEXP signif_threshold, SEXP posterior_threshold,
                  SEXP bayes_factor_threshold, SEXP prior_probability,
                  SEXP linked_probability, SEXP null_prevalence,
                  SEXP group1_linked_weight, SEXP bf_group1_alpha,
                  SEXP bf_group1_beta, SEXP bf_group2_alpha,
                  SEXP bf_group2_beta, SEXP bf_null_alpha,
                  SEXP bf_null_beta, SEXP group1, SEXP group2,
                  SEXP posterior_linked_family, SEXP posterior_linked_probability,
                   SEXP posterior_linked_alpha,
                  SEXP posterior_linked_beta, SEXP posterior_null_family,
                  SEXP posterior_null_probability, SEXP posterior_null_alpha, SEXP posterior_null_beta) {
  rsx_status_t st = rsx_triage(
      rsxr_str(table_path, "table_path"),
      rsxr_str(popmap_path, "popmap_path"),
      rsxr_str(output_file, "output_file"),
      (uint32_t) Rf_asInteger(min_depth),
      (float) Rf_asReal(signif_threshold),
      (double) Rf_asReal(posterior_threshold),
      (double) Rf_asReal(bayes_factor_threshold),
      (double) Rf_asReal(prior_probability),
      (double) Rf_asReal(linked_probability),
      (double) Rf_asReal(null_prevalence),
      (double) Rf_asReal(group1_linked_weight),
      (double) Rf_asReal(bf_group1_alpha),
      (double) Rf_asReal(bf_group1_beta),
      (double) Rf_asReal(bf_group2_alpha),
      (double) Rf_asReal(bf_group2_beta),
      (double) Rf_asReal(bf_null_alpha),
      (double) Rf_asReal(bf_null_beta),
      rsxr_str(group1, "group1"),
      rsxr_str(group2, "group2"),
      rsxr_str(posterior_linked_family, "posterior_linked_family"),
      (double) Rf_asReal(posterior_linked_probability),
      (double) Rf_asReal(posterior_linked_alpha),
      (double) Rf_asReal(posterior_linked_beta),
      rsxr_str(posterior_null_family, "posterior_null_family"),
      (double) Rf_asReal(posterior_null_probability),
      (double) Rf_asReal(posterior_null_alpha),
      (double) Rf_asReal(posterior_null_beta));
  rsxr_check(st);
  return R_NilValue;
}

SEXP C_rsx_depth(SEXP table_path, SEXP popmap_path, SEXP output_file,
                 SEXP min_frequency, SEXP streaming) {
  rsx_status_t st = rsx_depth(
      rsxr_str(table_path, "table_path"),
      rsxr_str(popmap_path, "popmap_path"),
      rsxr_str(output_file, "output_file"),
      (float) Rf_asReal(min_frequency),
      (bool) (Rf_asLogical(streaming) == TRUE));
  rsxr_check(st);
  return R_NilValue;
}

SEXP C_rsx_merge(SEXP input_files, SEXP output_file, SEXP buffer_size,
                 SEXP output_parquet) {
  if (TYPEOF(input_files) != STRSXP || LENGTH(input_files) < 1) {
    Rf_error("rsx: 'input_files' must be a non-empty character vector");
  }
  R_xlen_t n = XLENGTH(input_files);
  const char **arr = (const char **) R_alloc(n, sizeof(char *));
  for (R_xlen_t i = 0; i < n; i++) {
    if (STRING_ELT(input_files, i) == NA_STRING) {
      Rf_error("rsx: 'input_files' contains NA at position %lld", (long long) (i + 1));
    }
    arr[i] = CHAR(STRING_ELT(input_files, i));
  }
  rsx_status_t st = rsx_merge(
      arr, (size_t) n,
      rsxr_str(output_file, "output_file"),
      (size_t) Rf_asInteger(buffer_size),
      (bool) (Rf_asLogical(output_parquet) == TRUE));
  rsxr_check(st);
  return R_NilValue;
}

SEXP C_rsx_pca(SEXP table_path, SEXP output_dir, SEXP min_depth,
               SEXP n_components) {
  rsx_status_t st = rsx_pca(
      rsxr_str(table_path, "table_path"),
      rsxr_str(output_dir, "output_dir"),
      (uint32_t) Rf_asInteger(min_depth),
      (size_t) Rf_asInteger(n_components));
  rsxr_check(st);
  return R_NilValue;
}

SEXP C_rsx_version(void) {
  return Rf_mkString(RSX_VERSION);
}

static const R_CallMethodDef CallEntries[] = {
    {"C_rsx_process", (DL_FUNC) &C_rsx_process, 4},
    {"C_rsx_freq",    (DL_FUNC) &C_rsx_freq,    3},
    {"C_rsx_distrib", (DL_FUNC) &C_rsx_distrib, 28},
    {"C_rsx_signif",  (DL_FUNC) &C_rsx_signif,  29},
    {"C_rsx_triage",  (DL_FUNC) &C_rsx_triage,  27},
    {"C_rsx_depth",   (DL_FUNC) &C_rsx_depth,   5},
    {"C_rsx_merge",   (DL_FUNC) &C_rsx_merge,   4},
    {"C_rsx_pca",     (DL_FUNC) &C_rsx_pca,     4},
    {"C_rsx_version", (DL_FUNC) &C_rsx_version, 0},
    {NULL, NULL, 0}};

void R_init_rsxr(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
