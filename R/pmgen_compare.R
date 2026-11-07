#' Compare Two pmgen_result Objects
#'
#' Produces a side-by-side data frame comparing the summary statistics from two
#' `pmgen_result` objects at their shared N values, for example the results of
#' two [pmgen_simulate()] runs with different methods (such as ridge and lasso).
#'
#' @param result_1 A `pmgen_result` object.
#' @param result_2 A `pmgen_result` object.
#' @param stat Character string. Which summary statistic to compare: `"mean"`,
#'   `"sd"`, `"q025"`, `"q975"`, or `"degradation_mean"`. Default `"mean"`.
#'
#' @return A data frame with columns: `N`, `metric`, `<method_1>`, `<method_2>`,
#'   `difference` (result_1 minus result_2).
#'
#' @examples
#' \dontrun{
#' result_ridge <- pmgen_simulate(ref, tp, N = 456, method = "ridge", seed = 1)
#' result_lasso <- pmgen_simulate(ref, tp, N = 456, method = "lasso", seed = 1)
#' pmgen_compare(result_ridge, result_lasso)
#' }
#' @export
pmgen_compare <- function(result_1, result_2, stat = "mean") {

  if (!inherits(result_1, "pmgen_result")) {
    rlang::abort("'result_1' must be a pmgen_result object.")
  }
  if (!inherits(result_2, "pmgen_result")) {
    rlang::abort("'result_2' must be a pmgen_result object.")
  }

  stat <- match.arg(stat, c("mean", "sd", "q025", "q975", "degradation_mean"))

  shared_N <- intersect(result_1$N, result_2$N)
  if (length(shared_N) == 0L) {
    rlang::abort("The two results share no N values. Cannot compare.")
  }

  shared_metrics <- intersect(result_1$metrics, result_2$metrics)
  if (length(shared_metrics) == 0L) {
    rlang::abort("The two results share no metrics. Cannot compare.")
  }

  rows <- lapply(shared_N, function(n_val) {
    key <- as.character(n_val)
    s1  <- result_1$summaries[[key]]
    s2  <- result_2$summaries[[key]]

    # Filter to shared metrics
    s1 <- s1[s1$metric %in% shared_metrics, ]
    s2 <- s2[s2$metric %in% shared_metrics, ]

    # Align by metric
    s1 <- s1[order(s1$metric), ]
    s2 <- s2[order(s2$metric), ]

    data.frame(
      N          = n_val,
      metric     = s1$metric,
      v1         = s1[[stat]],
      v2         = s2[[stat]],
      difference = s1[[stat]] - s2[[stat]],
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  # Rename columns to reflect method labels
  names(out)[3L] <- result_1$method
  names(out)[4L] <- result_2$method

  out
}
