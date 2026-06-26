# Low-level path-based wrappers. Each maps one rsx CLI command to one C call
# and returns the output path invisibly. These mirror the pyrsx surface.

#' Process demultiplexed reads into a marker depth table
#'
#' @param input_dir Directory of demultiplexed FASTQ files.
#' @param output_file Path to write the marker depth table (TSV).
#' @param threads Number of worker threads.
#' @param min_depth Minimum depth for a marker to be retained.
#' @return The output path, invisibly.
#' @export
rsx_process <- function(input_dir, output_file, threads = 1L, min_depth = 1L) {
  .Call(C_rsx_process, input_dir, output_file,
        as.integer(threads), as.integer(min_depth))
  invisible(output_file)
}

#' Compute per-marker allele frequencies
#'
#' @param table_path Marker depth table (TSV).
#' @param output_file Output path (TSV).
#' @param min_depth Minimum depth threshold.
#' @return The output path, invisibly.
#' @export
rsx_freq <- function(table_path, output_file, min_depth = 1L) {
  .Call(C_rsx_freq, table_path, output_file, as.integer(min_depth))
  invisible(output_file)
}

#' Marker distribution across groups
#'
#' @param table_path Marker depth table (TSV).
#' @param popmap_path Population map (TSV).
#' @param output_file Output path (TSV).
#' @param min_depth Minimum depth threshold.
#' @param signif_threshold Significance threshold.
#' @param group1,group2 Group labels to contrast (e.g. "M", "F").
#' @param correction Multiple-testing correction: "bonferroni", "fdr", "none".
#' @param test Statistical test: "chisq", "fisher", "gtest".
#' @param output_bayes Also emit Bayesian columns.
#' @return The output path, invisibly.
#' @export
rsx_distrib <- function(table_path, popmap_path, output_file,
                        min_depth = 1L, signif_threshold = 0.05,
                        group1 = "", group2 = "",
                        correction = "bonferroni", test = "chisq",
                        output_bayes = FALSE) {
  .Call(C_rsx_distrib, table_path, popmap_path, output_file,
        as.integer(min_depth), as.numeric(signif_threshold),
        group1, group2, correction, test, as.logical(output_bayes))
  invisible(output_file)
}

#' Extract significant sex-linked markers
#'
#' @inheritParams rsx_distrib
#' @param output_fasta Also write a FASTA of significant markers.
#' @param bayes Also emit Bayesian columns.
#' @return The output path, invisibly.
#' @export
rsx_signif <- function(table_path, popmap_path, output_file,
                       min_depth = 1L, signif_threshold = 0.05,
                       group1 = "", group2 = "",
                       correction = "bonferroni", test = "chisq",
                       output_fasta = FALSE, bayes = FALSE) {
  .Call(C_rsx_signif, table_path, popmap_path, output_file,
        as.integer(min_depth), as.numeric(signif_threshold),
        group1, group2, correction, test,
        as.logical(output_fasta), as.logical(bayes))
  invisible(output_file)
}

#' Bayesian sex-linkage triage of markers
#'
#' @inheritParams rsx_distrib
#' @param posterior_threshold Posterior probability cutoff.
#' @param bayes_factor_threshold Bayes factor cutoff.
#' @param prior_probability Prior probability of sex linkage.
#' @param linked_probability Probability model for a linked marker.
#' @return The output path, invisibly.
#' @export
rsx_triage <- function(table_path, popmap_path, output_file,
                       min_depth = 10L, signif_threshold = 0.05,
                       posterior_threshold = 0.9,
                       bayes_factor_threshold = 10.0,
                       prior_probability = 0.01,
                       linked_probability = 0.9,
                       group1 = "M", group2 = "F") {
  .Call(C_rsx_triage, table_path, popmap_path, output_file,
        as.integer(min_depth), as.numeric(signif_threshold),
        as.numeric(posterior_threshold), as.numeric(bayes_factor_threshold),
        as.numeric(prior_probability), as.numeric(linked_probability),
        group1, group2)
  invisible(output_file)
}

#' Marker depth per individual
#'
#' @param table_path Marker depth table (TSV).
#' @param popmap_path Population map (TSV).
#' @param output_file Output path (TSV).
#' @param min_frequency Minimum group presence frequency.
#' @param streaming Use the streaming (bounded-memory) path.
#' @return The output path, invisibly.
#' @export
rsx_depth <- function(table_path, popmap_path, output_file,
                      min_frequency = 0.5, streaming = FALSE) {
  .Call(C_rsx_depth, table_path, popmap_path, output_file,
        as.numeric(min_frequency), as.logical(streaming))
  invisible(output_file)
}

#' Merge marker depth tables (bounded memory)
#'
#' @param input_files Character vector of marker table paths.
#' @param output_file Output path.
#' @param buffer_size External-sort buffer size; 0 uses the default.
#' @param output_parquet Write Parquet instead of TSV.
#' @return The output path, invisibly.
#' @export
rsx_merge <- function(input_files, output_file, buffer_size = 0L,
                      output_parquet = FALSE) {
  .Call(C_rsx_merge, as.character(input_files), output_file,
        as.integer(buffer_size), as.logical(output_parquet))
  invisible(output_file)
}

#' Streaming PCA of the marker matrix
#'
#' @param table_path Marker depth table (TSV).
#' @param output_dir Directory for PCA outputs.
#' @param min_depth Minimum depth threshold.
#' @param n_components Number of components; 0 uses the default.
#' @return The output directory, invisibly.
#' @export
rsx_pca <- function(table_path, output_dir, min_depth = 1L, n_components = 0L) {
  .Call(C_rsx_pca, table_path, output_dir,
        as.integer(min_depth), as.integer(n_components))
  invisible(output_dir)
}

#' rsx C library version
#'
#' @return The rsx version string the bindings were compiled against.
#' @export
#' @examples
#' rsx_version()
rsx_version <- function() {
  .Call(C_rsx_version)
}
