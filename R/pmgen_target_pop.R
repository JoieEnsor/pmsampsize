#' Generate a Target Population for Simulation
#'
#' Partitions the joint predictor distribution stored in a `pmgen_ref_model`
#' object into two non-overlapping sets: a fixed target population of size `n`
#' (used to evaluate each simulated model) and a development pool containing
#' all remaining rows (used to draw development datasets in [pmgen_simulate()]).
#'
#' The target population rows are sampled **without replacement** from
#' `joint_data`, so the two sets are guaranteed to be disjoint. The reference
#' model is applied to the target population rows to compute true predicted
#' risks and simulated binary outcomes.
#'
#' @param ref_model A `pmgen_ref_model` object from [pmgen_reference()].
#' @param n Integer. Number of individuals in the target population. Must be
#'   strictly less than `nrow(ref_model$joint_data)`. Default `100000`.
#' @param dev_pool_n Integer or NULL. If `NULL` (default), the development pool
#'   contains all rows of `joint_data` not selected for the target population.
#'   If an integer, only that many rows are retained in the pool (sampled from
#'   the remainder). Useful when memory is limited or when intentionally
#'   bounding the available development data.
#' @param seed Integer or NULL. Random seed.
#'
#' @return A `pmgen_target_pop` object with slots:
#'   - `$data`: target population data frame (predictors, `lp`, `p_true`,
#'     `outcome`)
#'   - `$dev_pool`: development pool data frame (predictors only; outcomes are
#'     generated fresh each simulation draw in [pmgen_simulate()])
#'   - `$ref_model`, `$n`, `$seed`
#'
#' @examples
#' \dontrun{
#' ref <- pmgen_reference(
#'   joint_data = my_data, prevalence = 0.68,
#'   cstatistic = 0.76, seed = 42
#' )
#' tp  <- pmgen_target_population(ref, n = 100000, seed = 42)
#' print(tp)
#' }
#' @export
pmgen_target_population <- function(ref_model,
                                  n          = 100000L,
                                  dev_pool_n = NULL,
                                  seed       = NULL) {

  if (!inherits(ref_model, "pmgen_ref_model")) {
    rlang::abort(paste0(
      "'ref_model' must be a pmgen_ref_model object ",
      "from pmgen_reference()."
    ))
  }

  n <- as.integer(n)

  joint_data <- ref_model$joint_data
  n_pool     <- nrow(joint_data)

  if (n >= n_pool) {
    rlang::abort(paste0(
      "'n' (", n, ") must be strictly less than nrow(joint_data) (",
      n_pool, "). The target population and development pool must be ",
      "partitioned without overlap. If using distributional assumptions, ",
      "increase 'n_synth' in pmgen_reference()."
    ))
  }

  if (!is.null(seed)) set.seed(seed)

  # Partition joint_data without replacement
  target_idx <- sample.int(n_pool, size = n, replace = FALSE)
  pool_idx   <- setdiff(seq_len(n_pool), target_idx)

  samp <- joint_data[target_idx, , drop = FALSE]
  rownames(samp) <- NULL

  dev_pool <- joint_data[pool_idx, , drop = FALSE]
  rownames(dev_pool) <- NULL

  # Optionally cap the development pool size
  if (!is.null(dev_pool_n)) {
    dev_pool_n <- as.integer(dev_pool_n)
    if (dev_pool_n >= nrow(dev_pool)) {
      rlang::warn(paste0(
        "'dev_pool_n' (", dev_pool_n, ") is >= the number of remaining rows ",
        "(", nrow(dev_pool), "). Using all remaining rows."
      ))
    } else {
      pool_keep <- sample.int(nrow(dev_pool), size = dev_pool_n,
                              replace = FALSE)
      dev_pool <- dev_pool[pool_keep, , drop = FALSE]
      rownames(dev_pool) <- NULL
    }
  }

  # Apply reference model to target population
  x_pred  <- as.matrix(samp[, ref_model$predictors, drop = FALSE])
  lp      <- ref_model$intercept + as.vector(x_pred %*% ref_model$coef)
  p_true  <- plogis(lp)
  outcome <- rbinom(n, 1L, p_true)

  samp$lp      <- lp
  samp$p_true  <- p_true
  samp$outcome <- outcome

  new_pmgen_target_pop(
    data      = samp,
    dev_pool  = dev_pool,
    ref_model = ref_model,
    n         = n,
    seed      = seed
  )
}
