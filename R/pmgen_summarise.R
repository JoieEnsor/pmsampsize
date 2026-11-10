# Summarisation and minimum-N identification

# ---- summary.pmgen_result ----------------------------------------------------

#' Summarise a pmgen_result object
#'
#' Returns a formatted data frame of posterior summary statistics across all
#' evaluated N values, one row per metric per N. Columns: N, metric, mean, sd,
#' q025 (2.5th percentile), q975 (97.5th percentile), degradation_mean.
#'
#' @param object A `pmgen_result` object from [pmgen_simulate()].
#' @param ... Ignored.
#' @return A data frame with one row per (N, metric) combination.
#' @export
summary.pmgen_result <- function(object, ...) {
  out <- do.call(rbind, lapply(object$N, function(n) {
    s <- object$summaries[[as.character(n)]]
    cbind(N = n, s, stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}


# ---- pmgen_min_n ---------------------------------------------------------------

#' Identify the Minimum Sample Size Meeting Specified Criteria
#'
#' Scans the posterior summaries and raw draws in a `pmgen_result` object and
#' returns the smallest N at which **all** user-specified criteria are
#' simultaneously met.
#'
#' @param result A `pmgen_result` object from [pmgen_simulate()].
#' @param criteria Named list of criterion specifications. The name must match a
#'   metric computed in `result` (e.g. `"cstat"`, `"cal_slope"`, `"mape"`,
#'   `"net_benefit"`, `"ervsi"`). Each element is itself a named list choosing
#'   **one** criterion type:
#'   \describe{
#'     \item{`range` + `min_prob`}{P(metric in `range`) >= `min_prob`. Use for
#'       calibration slope: `list(range = c(0.9, 1.1), min_prob = 0.90)`.}
#'     \item{`max_mean`}{mean(metric) <= `max_mean` (for lower-is-better metrics
#'       like MAPE, or absolute degradation). E.g. `list(max_mean = 0.05)`.}
#'     \item{`min_mean`}{mean(metric) >= `min_mean`. E.g. `list(min_mean = 0.70)`
#'       for c-statistic floor.}
#'     \item{`max_abs_deg_mean`}{|mean degradation| <= threshold. Uses
#'       `degradation_mean` from summaries. E.g. `list(max_abs_deg_mean = 0.025)`.}
#'     \item{`min_prob_above`}{P(metric >= `threshold`) >= `min_prob`. E.g. for
#'       ERVSI: `list(min_prob_above = 0.90, min_prob = 0.90)`.}
#'     \item{`min_prob_below`}{P(metric <= `threshold`) >= `min_prob`. E.g.
#'       `list(min_prob_below = 0.05, min_prob = 0.90)` for MAPE.}
#'   }
#' @param ... Reserved for future arguments.
#'
#' @return An object of class `pmgen_min_n`, a list with:
#'   \describe{
#'     \item{`min_n`}{The smallest N meeting all criteria, or `NA` if none do.}
#'     \item{`table`}{Data frame: one row per N, one column per criterion, plus
#'       `all_met` (logical) and `N`.}
#'     \item{`criteria`}{The criteria list as supplied.}
#'   }
#'
#' @examples
#' # Small settings so the example runs quickly; see pmgen_simulate()
#' ref <- pmgen_reference(
#'   joint_data = dvt_data, predictors = setdiff(names(dvt_data), "dvt"),
#'   prevalence = 0.21, cstatistic = 0.78,
#'   n_bisect_pop = 1e5, n_synth = 1e5, seed = 42
#' )
#' tp <- pmgen_target_population(ref, n = 5000, seed = 42)
#' result <- pmgen_simulate(ref, tp, N = c(250, 1000), n_sims = 20,
#'                          metrics = c("cstat", "cal_slope", "mape"),
#'                          seed = 1, verbose = FALSE)
#'
#' pmgen_min_n(result, criteria = list(
#'   cal_slope = list(range = c(0.8, 1.2), min_prob = 0.80),
#'   mape      = list(max_mean = 0.05)
#' ))
#' @references
#' Riley RD, Whittle R, Sadatsafavi M, Martin GP, Pate A, Collins GS, Ensor J. A general sample size framework
#' for developing or updating a predictive algorithm: with application to
#' clinical prediction models. *BMC Medical Research Methodology*.
#' 2026;26:121. \doi{10.1186/s12874-026-02856-7}
#' @export
pmgen_min_n <- function(result, criteria = list(), ...) {
  if (!inherits(result, "pmgen_result")) {
    rlang::abort("'result' must be a pmgen_result object.")
  }
  if (length(criteria) == 0) {
    rlang::abort("'criteria' must be a non-empty named list.")
  }
  if (is.null(names(criteria)) || any(names(criteria) == "")) {
    rlang::abort("Every element of 'criteria' must be named (the metric name).")
  }

  # Validate that all named metrics were computed
  bad <- setdiff(names(criteria), result$metrics)
  if (length(bad) > 0) {
    rlang::abort(paste0(
      "Criteria reference metric(s) not in result: ",
      paste(bad, collapse = ", "),
      ". Available: ", paste(result$metrics, collapse = ", ")
    ))
  }

  # Evaluate each criterion at each N
  n_vals      <- result$N
  crit_names  <- names(criteria)
  compliance  <- matrix(NA, nrow = length(n_vals), ncol = length(crit_names),
                        dimnames = list(as.character(n_vals), crit_names))

  for (n_chr in as.character(n_vals)) {
    summ  <- result$summaries[[n_chr]]   # data frame: metric, mean, ..., degradation_mean
    draws <- result$posterior[[n_chr]]   # data frame: n_sims x metrics

    for (crit_nm in crit_names) {
      spec <- criteria[[crit_nm]]
      compliance[n_chr, crit_nm] <- .evaluate_criterion(crit_nm, spec, summ, draws)
    }
  }

  all_met <- apply(compliance, 1, all)
  passing <- n_vals[all_met]
  min_n   <- if (length(passing) > 0) min(passing) else NA_integer_

  tbl <- as.data.frame(compliance, stringsAsFactors = FALSE)
  tbl <- cbind(N = n_vals, tbl, all_met = all_met)
  rownames(tbl) <- NULL

  structure(
    list(min_n = min_n, table = tbl, criteria = criteria, result = result),
    class = "pmgen_min_n"
  )
}


# ---- print.pmgen_min_n -------------------------------------------------------

#' @export
print.pmgen_min_n <- function(x, ...) {
  cat("pmgen minimum sample size\n")
  cat("-------------------------\n")
  if (is.na(x$min_n)) {
    cat("  Minimum N : none of the evaluated N values meets all criteria.\n")
    cat("  Consider expanding the N grid (larger values).\n")
  } else {
    cat(sprintf("  Minimum N : %d\n", x$min_n))
  }
  cat("\nCriterion compliance by N:\n")
  tbl <- x$table
  # Format logicals as tick/cross for readability
  for (col in names(x$criteria)) {
    tbl[[col]] <- ifelse(tbl[[col]], "YES", "no ")
  }
  tbl$all_met <- ifelse(x$table$all_met, "YES", "no ")
  print(tbl, row.names = FALSE)
  invisible(x)
}


# ---- Internal: evaluate one criterion ----------------------------------------

#' Evaluate a single criterion specification for one N value
#'
#' @param metric Character. Metric name.
#' @param spec Named list. One of the criterion forms documented in pmgen_min_n.
#' @param summ Data frame. Summary for this N (from summarise_posterior).
#' @param draws Data frame. Raw draws for this N (n_sims rows).
#' @return Logical scalar.
#' @noRd
.evaluate_criterion <- function(metric, spec, summ, draws) {
  row      <- summ[summ$metric == metric, , drop = FALSE]
  raw_vals <- draws[[metric]]

  if (nrow(row) == 0) {
    rlang::warn(paste0("Metric '", metric, "' not found in summaries; criterion set to FALSE."))
    return(FALSE)
  }

  # --- range + min_prob -------------------------------------------------------
  if (!is.null(spec$range) && !is.null(spec$min_prob)) {
    lo  <- spec$range[1]
    hi  <- spec$range[2]
    prob <- mean(raw_vals >= lo & raw_vals <= hi, na.rm = TRUE)
    return(prob >= spec$min_prob)
  }

  # --- max_mean ---------------------------------------------------------------
  if (!is.null(spec$max_mean)) {
    return(row$mean <= spec$max_mean)
  }

  # --- min_mean ---------------------------------------------------------------
  if (!is.null(spec$min_mean)) {
    return(row$mean >= spec$min_mean)
  }

  # --- max_abs_deg_mean -------------------------------------------------------
  if (!is.null(spec$max_abs_deg_mean)) {
    return(abs(row$degradation_mean) <= spec$max_abs_deg_mean)
  }

  # --- min_prob_above ---------------------------------------------------------
  if (!is.null(spec$min_prob_above) && !is.null(spec$min_prob)) {
    prob <- mean(raw_vals >= spec$min_prob_above, na.rm = TRUE)
    return(prob >= spec$min_prob)
  }

  # --- min_prob_below ---------------------------------------------------------
  if (!is.null(spec$min_prob_below) && !is.null(spec$min_prob)) {
    prob <- mean(raw_vals <= spec$min_prob_below, na.rm = TRUE)
    return(prob >= spec$min_prob)
  }

  rlang::abort(paste0(
    "Unrecognised criterion specification for metric '", metric, "'. ",
    "Must include one of: 'range'+'min_prob', 'max_mean', 'min_mean', ",
    "'max_abs_deg_mean', 'min_prob_above'+'min_prob', 'min_prob_below'+'min_prob'."
  ))
}
