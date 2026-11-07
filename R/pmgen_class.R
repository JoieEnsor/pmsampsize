# S3 class constructors and print/summary/plot methods

# ---- pmgen_ref_model ----------------------------------------------------------

#' Construct a pmgen_ref_model object
#' @noRd
new_pmgen_ref_model <- function(coef, intercept, prevalence_achieved,
                                cstat_achieved, predictors, joint_data,
                                unit_info_matrix, inv_unit_info,
                                call) {
  structure(
    list(
      coef                = coef,
      intercept           = intercept,
      prevalence_achieved = prevalence_achieved,
      cstat_achieved      = cstat_achieved,
      predictors          = predictors,
      joint_data          = joint_data,
      unit_info_matrix    = unit_info_matrix,
      inv_unit_info       = inv_unit_info,
      call                = call
    ),
    class = "pmgen_ref_model"
  )
}

#' @export
print.pmgen_ref_model <- function(x, ...) {
  cat("pmgen reference model\n")
  cat("---------------------\n")
  cat(sprintf("  Predictors        : %s\n", paste(x$predictors, collapse = ", ")))
  cat(sprintf("  Intercept         : %.4f\n", x$intercept))
  cat(sprintf("  Prevalence (achieved): %.4f\n", x$prevalence_achieved))
  cat(sprintf("  C-statistic (achieved): %.4f\n", x$cstat_achieved))
  cat(sprintf("  Coefficients      :\n"))
  print(round(x$coef, 4))
  cat("  Unit info matrix  : cached (", nrow(x$unit_info_matrix), "x",
      ncol(x$unit_info_matrix), ")\n", sep = "")
  invisible(x)
}


# ---- pmgen_target_pop --------------------------------------------------------

#' Construct a pmgen_target_pop object
#' @noRd
new_pmgen_target_pop <- function(data, dev_pool, ref_model, n, seed) {
  structure(
    list(
      data      = data,
      dev_pool  = dev_pool,
      ref_model = ref_model,
      n         = n,
      seed      = seed
    ),
    class = "pmgen_target_pop"
  )
}

#' @export
print.pmgen_target_pop <- function(x, ...) {
  cat("pmgen target population\n")
  cat("-----------------------\n")
  cat(sprintf("  Target pop N    : %d\n", x$n))
  cat(sprintf("  Dev pool N      : %d\n", nrow(x$dev_pool)))
  cat(sprintf("  Seed            : %s\n", if (is.null(x$seed)) "NULL" else x$seed))
  cat(sprintf("  Prevalence (observed): %.4f\n", mean(x$data$outcome)))
  cat(sprintf("  Mean true risk  : %.4f\n", mean(x$data$p_true)))
  cat(sprintf("  Predictors      : %s\n",
              paste(x$ref_model$predictors, collapse = ", ")))
  invisible(x)
}


# ---- pmgen_result ------------------------------------------------------------

#' Construct a pmgen_result object
#' @noRd
new_pmgen_result <- function(call, ref_model, N, method, metrics,
                             posterior, summaries, pred_matrices = NULL,
                             p_true = NULL, fit_info = NULL,
                             dropped_predictors = NULL) {
  structure(
    list(
      call               = call,
      ref_model          = ref_model,
      N                  = N,
      method             = method,
      metrics            = metrics,
      posterior          = posterior,
      summaries          = summaries,
      pred_matrices      = pred_matrices,
      p_true             = p_true,
      fit_info           = fit_info,
      dropped_predictors = dropped_predictors
    ),
    class = "pmgen_result"
  )
}

#' @export
print.pmgen_result <- function(x, ...) {
  cat("pmgen result\n")
  cat("------------\n")
  cat(sprintf("  Method              : %s\n", x$method))
  cat(sprintf("  N values            : %s\n", paste(x$N, collapse = ", ")))
  cat(sprintf("  Metrics             : %s\n",
              paste(x$metrics, collapse = ", ")))
  cat(sprintf("  Draws/N             : %d\n", nrow(x$posterior[[1]])))
  cat(sprintf("  Predictions stored  : %s\n",
              if (!is.null(x$pred_matrices)) "yes" else "no"))
  fi <- x$fit_info
  if (!is.null(fi) && any(fi$n_failed > 0)) {
    f <- fi[fi$n_failed > 0, ]
    cat(sprintf("  Failed fits         : %s\n",
                paste0("N = ", f$N, ": ", f$n_failed, "/", f$n_sims,
                       collapse = "; ")))
  }
  if (!is.null(fi) && any(fi$n_dropped > 0)) {
    cat("  Perfect prediction  : fits with predictors dropped\n")
    for (i in which(fi$n_dropped > 0)) {
      d <- x$dropped_predictors[i, ]
      d <- d[d > 0]
      cat(sprintf("    N = %d: %d/%d (%s)\n", fi$N[i], fi$n_dropped[i],
                  fi$n_sims[i], paste(names(d), d, collapse = ", ")))
    }
  }
  cat("\nSummary (means across simulations):\n")
  print(summary(x, ...), row.names = FALSE)
  invisible(x)
}

#' Plot a pmgen_result object
#'
#' Plots the results of [pmgen_simulate()]. Requires the \pkg{ggplot2} and
#' \pkg{scales} packages.
#'
#' @param x A `pmgen_result` object.
#' @param type Character. One of `"degradation"`, `"assurance"`,
#'   `"instability"`, `"calibration"`. The `"instability"` and `"calibration"`
#'   plots need the predictions stored (`store_predictions = TRUE` in
#'   [pmgen_simulate()]).
#' @param N Integer. For `"instability"` and `"calibration"`, which N value
#'   to plot. Defaults to the first N in the result.
#' @param ... Additional arguments passed to the specific plot function.
#' @return A ggplot object.
#' @export
plot.pmgen_result <- function(x, type = "degradation", N = NULL, ...) {
  rlang::check_installed(c("ggplot2", "scales"), reason = "to plot pmgen results.")
  type <- match.arg(type,
                    c("degradation", "assurance", "instability", "calibration"))
  switch(type,
    degradation  = .plot_degradation(x, ...),
    assurance    = .plot_assurance(x, ...),
    instability  = .plot_instability(x, N = N, ...),
    calibration  = .plot_calibration(x, N = N, ...)
  )
}
