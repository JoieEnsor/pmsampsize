# Internal performance metric functions
# All functions take vectors and return a scalar unless noted.
# None are exported; they are called from the simulation engine.

# ---- Discrimination ----------------------------------------------------------

#' C-statistic (AUC) via fast Wilcoxon rank-sum
#' @param p_hat Numeric vector of predicted probabilities.
#' @param outcome Integer vector of binary outcomes (0/1).
#' @noRd
compute_cstat <- function(p_hat, outcome) {
  .fast_auc(p_hat, outcome)
}


# ---- Calibration -------------------------------------------------------------

#' Calibration slope (logistic recalibration slope)
#'
#' Fits outcome ~ lp in a logistic model and returns the slope coefficient. A
#' slope of 1 indicates perfect calibration.
#' @param lp Numeric vector. The model's linear predictor (log-odds) on the
#'   target population; see [compute_metrics_one_draw()].
#' @param outcome Integer vector of binary outcomes (0/1).
#' @noRd
compute_cal_slope <- function(lp, outcome) {
  # Guard against perfect separation producing infinite coefficients
  tryCatch({
    fit <- glm(outcome ~ lp, family = binomial())
    unname(coef(fit)[2])
  }, warning = function(w) {
    # glm warns on boundary fits — still return the estimate
    fit <- suppressWarnings(glm(outcome ~ lp, family = binomial()))
    unname(coef(fit)[2])
  }, error = function(e) NA_real_)
}

#' Calibration intercept (calibration-in-the-large)
#' @inheritParams compute_cal_slope
#' @noRd
compute_cal_intercept <- function(lp, outcome) {
  tryCatch({
    fit <- glm(outcome ~ 1, offset = lp, family = binomial())
    unname(coef(fit)[1])
  }, error = function(e) NA_real_)
}

#' O:E ratio (observed / expected events)
#' @noRd
compute_oe_ratio <- function(p_hat, outcome) {
  sum(outcome) / sum(p_hat)
}


# ---- Prediction error --------------------------------------------------------

#' Mean Absolute Prediction Error vs. true risks
#'
#' Requires access to the true risk `p_true` from the target population.
#' @param p_hat Numeric vector of model predictions on target population.
#' @param p_true Numeric vector of true risks from reference model.
#' @noRd
compute_mape <- function(p_hat, p_true) {
  mean(abs(p_hat - p_true))
}

#' Root Mean Squared Prediction Error vs. true risks
#' @noRd
compute_rmspe <- function(p_hat, p_true) {
  sqrt(mean((p_hat - p_true)^2))
}

#' Brier score
#' @noRd
compute_brier <- function(p_hat, outcome) {
  mean((outcome - p_hat)^2)
}


# ---- Utility -----------------------------------------------------------------

#' Net benefit at a clinical decision threshold
#'
#' @param p_hat Numeric vector of predicted probabilities.
#' @param outcome Integer 0/1 vector.
#' @param threshold Numeric. Clinical probability threshold (0, 1).
#' @noRd
compute_net_benefit <- function(p_hat, outcome, threshold) {
  n    <- length(outcome)
  tp   <- sum((p_hat >= threshold) & (outcome == 1)) / n
  fp   <- sum((p_hat >= threshold) & (outcome == 0)) / n
  tp - fp * (threshold / (1 - threshold))
}

#' Expected Relative Value of Sample Information (ERVSI)
#'
#' ERVSI = net_benefit_model / net_benefit_perfect_model.
#' The perfect model's NB = prevalence - prevalence * threshold / (1 - threshold).
#' @noRd
compute_ervsi <- function(p_hat, outcome, threshold) {
  nb_model   <- compute_net_benefit(p_hat, outcome, threshold)
  prev       <- mean(outcome)
  nb_perfect <- prev - prev * threshold / (1 - threshold)
  if (nb_perfect <= 0) return(NA_real_)
  nb_model / nb_perfect
}


# ---- Degradation and assurance -----------------------------------------------

#' Degradation of a metric relative to its reference value
#' @param metric_val Numeric. Achieved metric value.
#' @param ref_val Numeric. Reference (ideal) metric value.
#' @noRd
compute_degradation <- function(metric_val, ref_val) {
  metric_val - ref_val
}

#' Assurance probability: P(criterion met) across simulation draws
#'
#' @param vals Numeric vector of per-simulation metric values.
#' @param criterion A function that takes a scalar and returns TRUE/FALSE.
#' @noRd
compute_assurance <- function(vals, criterion) {
  mean(vapply(vals, criterion, logical(1)), na.rm = TRUE)
}


# ---- Compute all requested metrics for one simulation draw -------------------

#' Compute all requested metrics from one simulation draw
#'
#' @param p_hat Numeric vector. Predictions on the target population.
#' @param target_pop A `pmgen_target_pop` object.
#' @param metrics Character vector of metric names to compute.
#' @param threshold Numeric. Clinical threshold for NB / ERVSI.
#' @param lp Numeric vector or `NULL`. The model's linear predictor on the
#'   target population, used for the calibration slope and intercept. If
#'   `NULL` (e.g. a custom method that returns only probabilities), it is
#'   recovered as `qlogis(p_hat)`, with `p_hat` first clipped to
#'   \[1e-15, 1 - 1e-15\] so that predictions of exactly 0 or 1 stay finite.
#' @return Named numeric vector.
#' @noRd
compute_metrics_one_draw <- function(p_hat, target_pop, metrics, threshold = 0.5,
                                     lp = NULL) {

  outcome <- target_pop$data$outcome
  p_true  <- target_pop$data$p_true

  if (is.null(lp) && any(c("cal_slope", "cal_intercept") %in% metrics)) {
    eps <- 1e-15
    lp  <- stats::qlogis(pmin(pmax(p_hat, eps), 1 - eps))
  }

  result <- numeric(0)

  if ("cstat" %in% metrics) {
    result["cstat"] <- compute_cstat(p_hat, outcome)
  }
  if ("cal_slope" %in% metrics) {
    result["cal_slope"] <- compute_cal_slope(lp, outcome)
  }
  if ("cal_intercept" %in% metrics) {
    result["cal_intercept"] <- compute_cal_intercept(lp, outcome)
  }
  if ("oe_ratio" %in% metrics) {
    result["oe_ratio"] <- compute_oe_ratio(p_hat, outcome)
  }
  if ("mape" %in% metrics) {
    result["mape"] <- compute_mape(p_hat, p_true)
  }
  if ("rmspe" %in% metrics) {
    result["rmspe"] <- compute_rmspe(p_hat, p_true)
  }
  if ("brier" %in% metrics) {
    result["brier"] <- compute_brier(p_hat, outcome)
  }
  if ("net_benefit" %in% metrics) {
    result["net_benefit"] <- compute_net_benefit(p_hat, outcome, threshold)
  }
  if ("ervsi" %in% metrics) {
    result["ervsi"] <- compute_ervsi(p_hat, outcome, threshold)
  }

  result
}


#' Summarise posterior draws for one N value
#'
#' @param draws Data frame of n_sims rows x metrics columns.
#' @param ref_metrics Named numeric vector of reference model metric values
#'   (computed once on the full target population).
#' @return Data frame of summary statistics: mean, sd, q025, q975, degradation.
#' @noRd
summarise_posterior <- function(draws, ref_metrics = NULL) {
  metrics <- names(draws)

  do.call(rbind, lapply(metrics, function(m) {
    vals <- draws[[m]]
    row <- data.frame(
      metric  = m,
      mean    = mean(vals, na.rm = TRUE),
      sd      = sd(vals,   na.rm = TRUE),
      q025    = quantile(vals, 0.025, na.rm = TRUE),
      q975    = quantile(vals, 0.975, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    if (!is.null(ref_metrics) && m %in% names(ref_metrics)) {
      row$degradation_mean <- mean(vals - ref_metrics[m], na.rm = TRUE)
    } else {
      row$degradation_mean <- NA_real_
    }
    row
  }))
}
