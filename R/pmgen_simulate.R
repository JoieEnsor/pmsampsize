#' Run the Fully Simulation-Based Sample Size Pathway
#'
#' For each value in `N`, draws `n_sims` development datasets of size N from
#' the joint predictor distribution, fits the chosen model on each, evaluates
#' it on the fixed target population, and computes requested performance metrics.
#' Results are summarised into a `pmgen_result` object.
#'
#' @param ref_model A `pmgen_ref_model` object from [pmgen_reference()].
#' @param target_pop A `pmgen_target_pop` object from [pmgen_target_population()].
#' @param N Integer vector of development sample sizes to evaluate.
#' @param n_sims Integer. Number of simulation draws per N value. Default `1000`.
#' @param method Character string naming a built-in method (`"logistic"`,
#'   `"ridge"`, `"lasso"`, `"elasticnet"`, `"firth"`), or a user-supplied
#'   function with signature `(train_data, predictors, outcome_col, method_args)`
#'   that returns a prediction closure. See Details.
#' @param method_args Named list of method-specific tuning arguments passed
#'   through to the fit function. For `"logistic"`: `shrinkage = TRUE` applies
#'   van Houwelingen uniform shrinkage; `separation = "drop"` (default) or
#'   `"keep"` sets how perfect prediction is handled (see Details). For
#'   `"ridge"`, `"lasso"`, `"elasticnet"`:
#'   `nfolds`, `type.measure`, `lambda`. For `"elasticnet"`: `alpha` (required).
#'   For `"firth"`: `maxit`, `maxstep`.
#' @param threshold Numeric in (0, 1). Clinical probability threshold used for
#'   net benefit and ERVSI calculations. Default `0.5`.
#' @param metrics Character vector of metrics to compute. Any of: `"cstat"`,
#'   `"cal_slope"`, `"cal_intercept"`, `"oe_ratio"`, `"mape"`, `"rmspe"`,
#'   `"brier"`, `"net_benefit"`, `"ervsi"`. Default:
#'   `c("cstat", "cal_slope", "mape", "net_benefit")`.
#' @param seed Integer or NULL. Sets the global RNG seed before the simulation
#'   loop for reproducibility. Default `NULL`.
#' @param store_predictions Logical. If `TRUE`, the full matrix of predicted
#'   probabilities on the target population is stored for every simulation draw
#'   (a `n_pop x n_sims` matrix per N value). Required for
#'   `plot(result, type = "instability")` and `plot(result, type = "calibration")`.
#'   Memory use: approximately `n_pop * n_sims * 8` bytes per N value (~750 MB
#'   for n_pop = 100,000 and n_sims = 1,000). Default `FALSE`.
#' @param verbose Logical. Print progress messages per N value. Default `TRUE`.
#'
#' @details
#' **Parallel execution:** If the \pkg{furrr} package is installed, the
#' simulations are run with `furrr::future_map()`, so they follow whatever
#' \pkg{future} plan is set. Set a parallel backend before calling
#' `pmgen_simulate()` for speed:
#' ```r
#' future::plan(future::multisession, workers = 4)
#' result <- pmgen_simulate(...)
#' future::plan(future::sequential)  # reset afterwards
#' ```
#' Without a parallel plan, or without \pkg{furrr}, execution is sequential.
#' Each simulation's seed is drawn before the loop, so results for a given
#' `seed` are identical whether run sequentially or in parallel.
#'
#' The `"ridge"`, `"lasso"` and `"elasticnet"` methods need the \pkg{glmnet}
#' package, and `"firth"` needs \pkg{logistf}.
#'
#' **Custom methods:** Supply a function of the form
#' ```r
#' my_method <- function(train_data, predictors, outcome_col, method_args) {
#'   # fit some model on train_data
#'   function(new_data) {
#'     # return vector of predicted probabilities
#'   }
#' }
#' result <- pmgen_simulate(..., method = my_method)
#' ```
#'
#' **Development data:** At each simulation iteration, N rows are sampled with
#' replacement from `ref_model$joint_data`. Binary outcomes are generated from
#' the reference model's true risks (`Y ~ Bernoulli(plogis(Xβ))`). A warning
#' is issued if N exceeds the number of rows in `joint_data`.
#'
#' **Perfect prediction (`method = "logistic"`):** In small development
#' samples, one level of a binary predictor can have the same outcome for
#' everyone (e.g. all treated patients have the event). Maximum likelihood then
#' gives that predictor an effectively infinite coefficient, which distorts the
#' model's predictions. With `method_args = list(separation = "drop")` (the
#' default), before each fit any two-valued predictor with a level whose outcome
#' is constant is dropped, together with the observations at that level; this
#' is repeated until no such predictor remains, and the remaining predictors
#' are fitted. Dropped predictors get a coefficient of 0 when predicting. If a
#' predictor with more than two values separates the outcome completely, the
#' fit fails and that simulation's metrics are `NA`. Use `separation = "keep"`
#' to fit [stats::glm()] to all predictors and observations as they are. The
#' penalised methods (`"ridge"`, `"lasso"`, `"elasticnet"`, `"firth"`) handle
#' separation themselves and ignore this argument.
#'
#' **Failed fits:** If fitting or prediction fails in a simulation, its metrics
#' are `NA` and are left out of the summaries. A warning gives the number of
#' failed fits for each N and the first error message.
#'
#' @return A `pmgen_result` object. Use [summary()], [print()], and
#'   [plot()] methods to inspect results, and [pmgen_min_n()] to identify the
#'   minimum sample size meeting specified criteria. Element `fit_info` is a
#'   data frame giving, for each N, the number of failed fits (`n_failed`) and
#'   of fits in which predictors were dropped for perfect prediction
#'   (`n_dropped`); `dropped_predictors` is a matrix (one row per N) of how
#'   often each predictor was dropped.
#'
#' @examples
#' \dontrun{
#' ref <- pmgen_reference(joint_data = my_data, prevalence = 0.3,
#'                      cstatistic = 0.75, seed = 1)
#' tp  <- pmgen_target_population(ref, n = 50000, seed = 1)
#'
#' future::plan(future::multisession, workers = 4)
#' result <- pmgen_simulate(ref, tp, N = c(200, 500, 1000), n_sims = 500,
#'                        method = "lasso", seed = 42)
#' future::plan(future::sequential)
#'
#' summary(result)
#' }
#' @export
pmgen_simulate <- function(ref_model,
                          target_pop,
                          N,
                          n_sims            = 1000L,
                          method            = "logistic",
                          method_args       = list(),
                          threshold         = 0.5,
                          metrics           = c("cstat", "cal_slope", "mape",
                                                "net_benefit"),
                          store_predictions = FALSE,
                          seed              = NULL,
                          verbose           = TRUE) {

  call <- match.call()

  # ---- Input validation -------------------------------------------------------
  if (!inherits(ref_model, "pmgen_ref_model")) {
    rlang::abort("'ref_model' must be a pmgen_ref_model object.")
  }
  if (!inherits(target_pop, "pmgen_target_pop")) {
    rlang::abort("'target_pop' must be a pmgen_target_pop object.")
  }
  N <- as.integer(N)
  if (any(N <= 0)) rlang::abort("All values of N must be positive integers.")

  n_sims <- as.integer(n_sims)
  if (n_sims < 1L) rlang::abort("'n_sims' must be >= 1.")

  valid_metrics <- c("cstat", "cal_slope", "cal_intercept", "oe_ratio",
                     "mape", "rmspe", "brier", "net_benefit", "ervsi")
  bad_metrics <- setdiff(metrics, valid_metrics)
  if (length(bad_metrics) > 0) {
    rlang::abort(paste0("Unknown metric(s): ", paste(bad_metrics, collapse = ", "),
                        ". Valid options: ", paste(valid_metrics, collapse = ", ")))
  }

  if (threshold <= 0 || threshold >= 1) {
    rlang::abort("'threshold' must be strictly between 0 and 1.")
  }

  # Warn if any N exceeds the development pool size
  n_pool <- nrow(target_pop$dev_pool)
  large_N <- N[N > n_pool]
  if (length(large_N) > 0) {
    rlang::warn(paste0(
      "Some N values (", paste(large_N, collapse = ", "), ") exceed the ",
      "number of rows in the development pool (", n_pool, "). Development ",
      "samples will be drawn with replacement, which may inflate apparent ",
      "information at small pool sizes."
    ))
  }

  # ---- Resolve method ---------------------------------------------------------
  fit_fn <- resolve_method(method)
  method_label <- if (is.function(method)) "custom" else method

  # Fit errors are caught per simulation, so check method_args$separation here
  if (!is.null(method_args$separation)) {
    if (!identical(method_label, "logistic")) {
      rlang::warn("method_args$separation applies only to method = 'logistic' and is ignored.")
    } else if (!(is.character(method_args$separation) &&
                 length(method_args$separation) == 1L &&
                 method_args$separation %in% c("drop", "keep"))) {
      rlang::abort("method_args$separation must be \"drop\" or \"keep\".")
    }
  }

  # ---- Reference metrics (evaluate reference model once on target_pop) --------
  ref_phat <- target_pop$data$p_true  # reference model predictions = true risks
  ref_metrics <- compute_metrics_one_draw(ref_phat, target_pop, metrics, threshold,
                                          lp = target_pop$data$lp)

  # ---- Seed -------------------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)

  # ---- Main loop over N values ------------------------------------------------
  posterior_list <- vector("list", length(N))
  summary_list   <- vector("list", length(N))
  pred_mat_list  <- vector("list", length(N))
  fit_info_list  <- vector("list", length(N))
  names(posterior_list) <- names(summary_list) <-
    names(pred_mat_list) <- names(fit_info_list) <- as.character(N)

  predictors  <- ref_model$predictors
  outcome_col <- ".outcome"

  for (n_val in N) {
    if (verbose) {
      rlang::inform(
        sprintf("Simulating N = %d (%d iterations)...", n_val, n_sims)
      )
    }

    # Seeds per draw: reproducible regardless of parallel backend
    sim_seeds <- sample.int(.Machine$integer.max, n_sims)

    one_sim <- function(i) {
      run_single_sim(
        n_val             = n_val,
        target_pop        = target_pop,
        ref_model         = ref_model,
        fit_fn            = fit_fn,
        method_args       = method_args,
        metrics           = metrics,
        threshold         = threshold,
        predictors        = predictors,
        outcome_col       = outcome_col,
        sim_seed          = sim_seeds[i],
        store_predictions = store_predictions
      )
    }
    # Each simulation sets its own seed, so both routes give identical results
    draws_list <- if (.use_furrr()) {
      furrr::future_map(seq_len(n_sims), one_sim,
                        .options = furrr::furrr_options(seed = NULL))
    } else {
      .lapply_keep_rng(seq_len(n_sims), one_sim)
    }

    # Each element is list(metrics, p_hat, error, dropped); see run_single_sim()
    draws_df <- as.data.frame(
      do.call(rbind, lapply(draws_list, `[[`, "metrics"))
    )
    if (store_predictions) {
      pred_mat_list[[as.character(n_val)]] <- do.call(
        cbind, lapply(draws_list, `[[`, "p_hat")
      )  # n_pop x n_sims
    }

    posterior_list[[as.character(n_val)]] <- draws_df
    summary_list[[as.character(n_val)]]   <-
      summarise_posterior(draws_df, ref_metrics)

    info <- .fit_info_one_n(draws_list, predictors)
    fit_info_list[[as.character(n_val)]] <- info
    if (verbose && info$n_dropped > 0) {
      d <- info$dropped[info$dropped > 0]
      rlang::inform(sprintf(
        "  Perfect prediction in %d of %d fits; predictors dropped: %s.",
        info$n_dropped, n_sims, paste(names(d), d, collapse = ", ")
      ))
    }
  }

  if (verbose) rlang::inform("Done.")

  fit_info <- data.frame(
    N         = N,
    n_sims    = n_sims,
    n_failed  = vapply(fit_info_list, `[[`, integer(1), "n_failed"),
    n_dropped = vapply(fit_info_list, `[[`, integer(1), "n_dropped"),
    row.names = NULL
  )
  # N x predictors matrix: how often each predictor was dropped
  dropped_predictors <- do.call(rbind, lapply(fit_info_list, `[[`, "dropped"))
  rownames(dropped_predictors) <- as.character(N)

  if (any(fit_info$n_failed > 0)) {
    f <- fit_info[fit_info$n_failed > 0, ]
    first_error <- Filter(Negate(is.na),
                          lapply(fit_info_list, `[[`, "first_error"))[[1]]
    rlang::warn(paste0(
      "Model fitting failed in some simulations (",
      paste0("N = ", f$N, ": ", f$n_failed, " of ", n_sims, collapse = "; "),
      "). Their metrics are NA and are excluded from the summaries. ",
      "First error: ", first_error
    ))
  }

  new_pmgen_result(
    call          = call,
    ref_model     = ref_model,
    N             = N,
    method        = method_label,
    metrics       = metrics,
    posterior     = posterior_list,
    summaries     = summary_list,
    pred_matrices = if (store_predictions) pred_mat_list else NULL,
    p_true        = if (store_predictions) target_pop$data$p_true else NULL,
    fit_info           = fit_info,
    dropped_predictors = dropped_predictors
  )
}


# ---- Internal simulation engine -----------------------------------------------

#' Run one simulation draw
#'
#' Samples a development dataset of size n_val, fits a model, evaluates on the
#' target population, and computes the requested metrics.
#'
#' @return A list with `metrics` (named numeric vector, `NA` if the fit or
#'   prediction failed), `p_hat` (predictions on the target population if
#'   `store_predictions = TRUE`, otherwise `NULL`), `error` (the error message
#'   if the fit or prediction failed, otherwise `NULL`) and `dropped` (names of
#'   predictors dropped by the fit function for perfect prediction).
#' @noRd
run_single_sim <- function(n_val, target_pop, ref_model, fit_fn, method_args,
                            metrics, threshold, predictors, outcome_col,
                            sim_seed, store_predictions = FALSE) {

  set.seed(sim_seed)

  dev_pool <- target_pop$dev_pool
  n_pool   <- nrow(dev_pool)

  # 1. Sample N rows of predictors from the development pool (with replacement)
  idx      <- sample.int(n_pool, size = n_val, replace = TRUE)
  dev_pred <- dev_pool[idx, predictors, drop = FALSE]
  rownames(dev_pred) <- NULL

  # 2. Generate outcomes from the reference model's true risks
  lp      <- ref_model$intercept + as.vector(as.matrix(dev_pred) %*% ref_model$coef)
  p_true  <- stats::plogis(lp)
  outcome <- stats::rbinom(n_val, 1L, p_true)

  dev_pred[[outcome_col]] <- outcome

  # Returned when fitting or prediction fails: NA metrics plus the error message
  failed <- function(e) {
    list(metrics = stats::setNames(rep(NA_real_, length(metrics)), metrics),
         p_hat   = if (store_predictions) rep(NA_real_, nrow(target_pop$data)),
         error   = conditionMessage(e),
         dropped = character(0))
  }

  # 3. Fit the model; catch failures gracefully
  pred_fn <- tryCatch(
    fit_fn(dev_pred, predictors, outcome_col, method_args),
    error = function(e) e
  )
  if (inherits(pred_fn, "error")) return(failed(pred_fn))

  # 4. Predict on target population
  p_hat <- tryCatch(
    pred_fn(target_pop$data),
    error = function(e) e
  )
  if (inherits(p_hat, "error")) return(failed(p_hat))

  # Built-in methods attach the linear predictor; custom methods may not
  lp    <- attr(p_hat, "lp")
  p_hat <- as.vector(p_hat)

  # 5. Compute metrics
  metric_vals <- compute_metrics_one_draw(p_hat, target_pop, metrics, threshold,
                                          lp = lp)

  list(metrics = metric_vals,
       p_hat   = if (store_predictions) p_hat,
       error   = NULL,
       dropped = attr(pred_fn, "dropped") %||% character(0))
}


#' Run the simulations through furrr (respecting any future plan) if it is
#' installed; otherwise pmgen_simulate() falls back to lapply()
#' @noRd
.use_furrr <- function() {
  requireNamespace("furrr", quietly = TRUE)
}

#' lapply() that leaves the global RNG state as it found it
#'
#' Each simulation calls set.seed(), which would otherwise move the session's
#' random number stream and so change the seeds drawn for the next N. furrr
#' evaluates the simulations without touching the session's stream, so
#' restoring it here makes the sequential fallback give identical results.
#' @noRd
.lapply_keep_rng <- function(X, FUN) {
  seed_env <- globalenv()
  old_seed <- get0(".Random.seed", envir = seed_env, inherits = FALSE)
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = seed_env, inherits = FALSE)) {
        rm(".Random.seed", envir = seed_env)
      }
    } else {
      assign(".Random.seed", old_seed, envir = seed_env)
    }
  }, add = TRUE)
  lapply(X, FUN)
}


#' Summarise failed fits and perfect-prediction drops for one N value
#'
#' @param draws_list List of [run_single_sim()] results.
#' @param predictors Character vector of predictor names.
#' @return A list with `n_failed`, `n_dropped` (number of fits in which at
#'   least one predictor was dropped), `dropped` (named integer vector: how
#'   often each predictor was dropped) and `first_error` (the first error
#'   message, or `NA`).
#' @noRd
.fit_info_one_n <- function(draws_list, predictors) {
  errors  <- unlist(lapply(draws_list, `[[`, "error"))
  dropped <- lapply(draws_list, `[[`, "dropped")
  list(
    n_failed    = length(errors),
    n_dropped   = sum(lengths(dropped) > 0),
    dropped     = vapply(predictors, function(p) {
      sum(vapply(dropped, function(d) p %in% d, logical(1)))
    }, integer(1)),
    first_error = if (length(errors) > 0) errors[1] else NA_character_
  )
}
