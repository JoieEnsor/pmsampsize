# Visualisation functions for pmgen_result objects.
# Dispatched from plot.pmgen_result() in pmgen_class.R.

#' @importFrom rlang .data
NULL

# ---- Shared helpers ----------------------------------------------------------

.metric_labels <- c(
  cstat          = "C-statistic",
  cal_slope      = "Calibration slope",
  cal_intercept  = "Calibration intercept",
  oe_ratio       = "O:E ratio",
  mape           = "MAPE",
  rmspe          = "RMSPE",
  brier          = "Brier score",
  net_benefit    = "Net benefit",
  ervsi          = "ERVSI",
  interval_width = "95% UI width"
)

.label_metric <- function(m) {
  ifelse(m %in% names(.metric_labels), .metric_labels[m], m)
}

.check_pred_matrices <- function(result, type) {
  if (is.null(result$pred_matrices)) {
    rlang::abort(paste0(
      "Prediction matrices not found. Re-run pmgen_simulate() ",
      "with store_predictions = TRUE to enable '", type, "' plots."
    ))
  }
}

.resolve_n <- function(result, N) {
  if (is.null(N)) {
    n_val <- result$N[1]
    rlang::inform(sprintf(
      "No N specified; using N = %d (first value in result).", n_val
    ))
    return(n_val)
  }
  N <- as.integer(N)
  if (!N %in% result$N) {
    rlang::abort(sprintf(
      "N = %d not found in result. Available: %s.",
      N, paste(result$N, collapse = ", ")
    ))
  }
  N
}


# ---- degradation plot --------------------------------------------------------

#' @noRd
.plot_degradation <- function(result, ...) {
  summ <- summary(result)

  has_deg <- !is.na(summ$degradation_mean)
  if (!any(has_deg)) {
    rlang::abort(paste0(
      "No degradation values available. Ensure the result was produced with ",
      "a reference model so that ref_metrics are computed."
    ))
  }

  df <- summ[has_deg, ]
  df$metric_label <- .label_metric(df$metric)
  df$ref_val      <- df$mean - df$degradation_mean
  df$deg_lo       <- df$q025 - df$ref_val
  df$deg_hi       <- df$q975 - df$ref_val

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x      = .data$N,
      y      = .data$degradation_mean,
      colour = .data$metric_label,
      fill   = .data$metric_label
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "grey50") +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$deg_lo, ymax = .data$deg_hi),
      alpha = 0.15, colour = NA
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~ metric_label, scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = unique(df$N)) +
    ggplot2::labs(
      title  = paste0("Performance degradation vs. reference model (",
                      result$method, ")"),
      x      = "Development sample size (N)",
      y      = "Degradation (metric \u2212 reference)",
      colour = NULL,
      fill   = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
}


# ---- assurance plot ----------------------------------------------------------

.default_assurance_fns <- list(
  cstat = function(vals, ref_val) {
    mean(vals >= ref_val - 0.025, na.rm = TRUE)
  },
  cal_slope = function(vals, ref_val) {
    mean(vals >= 0.9 & vals <= 1.1, na.rm = TRUE)
  },
  cal_intercept = function(vals, ref_val) {
    mean(abs(vals) <= 0.5, na.rm = TRUE)
  },
  oe_ratio = function(vals, ref_val) {
    mean(vals >= 0.9 & vals <= 1.1, na.rm = TRUE)
  },
  mape = function(vals, ref_val) {
    mean(vals <= 0.10, na.rm = TRUE)
  },
  rmspe = function(vals, ref_val) {
    mean(vals <= 0.10, na.rm = TRUE)
  },
  brier = function(vals, ref_val) {
    mean(vals <= ref_val + 0.05, na.rm = TRUE)
  },
  net_benefit = function(vals, ref_val) {
    mean(vals > 0, na.rm = TRUE)
  },
  ervsi = function(vals, ref_val) {
    mean(vals >= 0.9, na.rm = TRUE)
  }
)

#' @noRd
.plot_assurance <- function(result, ref_line = 0.9, ...) {

  metrics    <- result$metrics
  n_vals     <- result$N
  first_summ <- result$summaries[[as.character(n_vals[1])]]

  rows <- do.call(rbind, lapply(metrics, function(m) {
    if (!m %in% names(.default_assurance_fns)) return(NULL)
    crit_fn <- .default_assurance_fns[[m]]
    ref_row <- first_summ[first_summ$metric == m, , drop = FALSE]
    ref_val <- if (nrow(ref_row) > 0 && !is.na(ref_row$degradation_mean)) {
      ref_row$mean - ref_row$degradation_mean
    } else {
      NA_real_
    }

    do.call(rbind, lapply(n_vals, function(n) {
      vals  <- result$posterior[[as.character(n)]][[m]]
      assur <- crit_fn(vals, ref_val)
      data.frame(
        N            = n,
        metric       = m,
        metric_label = .label_metric(m),
        assurance    = assur,
        stringsAsFactors = FALSE
      )
    }))
  }))

  if (is.null(rows) || nrow(rows) == 0) {
    rlang::abort("No metrics with default assurance criteria found in result.")
  }

  ggplot2::ggplot(
    rows,
    ggplot2::aes(
      x      = .data$N,
      y      = .data$assurance,
      colour = .data$metric_label
    )
  ) +
    ggplot2::geom_hline(yintercept = ref_line, linetype = "dashed",
                        colour = "grey40") +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~ metric_label) +
    ggplot2::scale_x_continuous(breaks = unique(rows$N)) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::labs(
      title  = paste0("Assurance probabilities by N (", result$method, ")"),
      x      = "Development sample size (N)",
      y      = "Assurance probability",
      colour = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
}


# ---- instability plot --------------------------------------------------------

#' @noRd
.plot_instability <- function(result, N = NULL, n_sample = 500L, ...) {

  .check_pred_matrices(result, "instability")

  n_val  <- .resolve_n(result, N)
  pmat   <- result$pred_matrices[[as.character(n_val)]]
  p_true <- result$p_true

  p_lo <- apply(pmat, 1, quantile, 0.025, na.rm = TRUE)
  p_hi <- apply(pmat, 1, quantile, 0.975, na.rm = TRUE)
  iw   <- p_hi - p_lo

  df    <- data.frame(p_true = p_true, iw = iw)
  n_pop <- nrow(df)
  if (n_pop > n_sample) {
    df <- df[sample.int(n_pop, n_sample), ]
  }

  ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$p_true, y = .data$iw, colour = .data$iw)
  ) +
    ggplot2::geom_point(alpha = 0.5, size = 1.2) +
    ggplot2::scale_colour_viridis_c(option = "plasma", direction = -1,
                                    name = "95% UI\nwidth") +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::scale_y_continuous(limits = c(0, NA)) +
    ggplot2::labs(
      title    = sprintf("Prediction instability at N = %d (%s)",
                         n_val, result$method),
      subtitle = sprintf("%d individuals shown (subsampled from %d)",
                         nrow(df), n_pop),
      x = "True risk (reference model)",
      y = "95% interval width across simulations"
    ) +
    ggplot2::theme_bw()
}


# ---- calibration fan plot ----------------------------------------------------

#' @noRd
.plot_calibration <- function(result, N = NULL, n_curves = 100L,
                              n_groups = 10L, ...) {

  .check_pred_matrices(result, "calibration")

  n_val         <- .resolve_n(result, N)
  pmat          <- result$pred_matrices[[as.character(n_val)]]
  p_true        <- result$p_true
  n_sims        <- ncol(pmat)
  n_curves_plot <- min(as.integer(n_curves), n_sims)
  curve_idx     <- sample.int(n_sims, n_curves_plot)

  curve_dfs <- lapply(curve_idx, function(s) {
    p_hat  <- pmat[, s]
    breaks <- quantile(p_hat,
                       probs = seq(0, 1, length.out = n_groups + 1L),
                       na.rm = TRUE)
    grp <- cut(p_hat, breaks = breaks, include.lowest = TRUE,
               labels = FALSE)
    agg <- stats::aggregate(
      cbind(p_hat, p_true) ~ grp,
      data = data.frame(p_hat = p_hat, p_true = p_true, grp = grp),
      FUN  = mean, na.rm = TRUE
    )
    agg$sim <- s
    agg
  })

  df <- do.call(rbind, curve_dfs)

  ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$p_hat, y = .data$p_true, group = .data$sim)
  ) +
    ggplot2::geom_abline(intercept = 0, slope = 1, colour = "black",
                         linewidth = 0.8, linetype = "dashed") +
    ggplot2::geom_line(alpha = 0.15, colour = "#2166ac", linewidth = 0.4) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::labs(
      title    = sprintf("Calibration instability at N = %d (%s)",
                         n_val, result$method),
      subtitle = sprintf("%d calibration curves overlaid", n_curves_plot),
      x        = "Mean predicted risk (decile bins)",
      y        = "Mean true risk (decile bins)"
    ) +
    ggplot2::theme_bw()
}
