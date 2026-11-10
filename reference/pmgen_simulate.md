# Run the Fully Simulation-Based Sample Size Pathway

For each value in `N`, draws `n_sims` development datasets of size N
from the joint predictor distribution, fits the chosen model on each,
evaluates it on the fixed target population, and computes requested
performance metrics. Results are summarised into a `pmgen_result`
object.

## Usage

``` r
pmgen_simulate(
  ref_model,
  target_pop,
  N,
  n_sims = 1000L,
  method = "logistic",
  method_args = list(),
  threshold = 0.5,
  metrics = c("cstat", "cal_slope", "mape", "net_benefit"),
  store_predictions = FALSE,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- ref_model:

  A `pmgen_ref_model` object from
  [`pmgen_reference()`](https://joieensor.github.io/pmsampsize/reference/pmgen_reference.md).

- target_pop:

  A `pmgen_target_pop` object from
  [`pmgen_target_population()`](https://joieensor.github.io/pmsampsize/reference/pmgen_target_population.md).

- N:

  Integer vector of development sample sizes to evaluate.

- n_sims:

  Integer. Number of simulation draws per N value. Default `1000`.

- method:

  Character string naming a built-in method (`"logistic"`, `"ridge"`,
  `"lasso"`, `"elasticnet"`, `"firth"`), or a user-supplied function
  with signature `(train_data, predictors, outcome_col, method_args)`
  that returns a prediction closure. See Details.

- method_args:

  Named list of method-specific tuning arguments passed through to the
  fit function. For `"logistic"`: `shrinkage = TRUE` applies van
  Houwelingen uniform shrinkage; `separation = "drop"` (default) or
  `"keep"` sets how perfect prediction is handled (see Details). For
  `"ridge"`, `"lasso"`, `"elasticnet"`: `nfolds`, `type.measure`,
  `lambda`. For `"elasticnet"`: `alpha` (required). For `"firth"`:
  `maxit`, `maxstep`.

- threshold:

  Numeric in (0, 1). Clinical probability threshold used for net benefit
  and ERVSI calculations. Default `0.5`.

- metrics:

  Character vector of metrics to compute. Any of: `"cstat"`,
  `"cal_slope"`, `"cal_intercept"`, `"oe_ratio"`, `"mape"`, `"rmspe"`,
  `"brier"`, `"net_benefit"`, `"ervsi"`. Default:
  `c("cstat", "cal_slope", "mape", "net_benefit")`.

- store_predictions:

  Logical. If `TRUE`, the full matrix of predicted probabilities on the
  target population is stored for every simulation draw (a
  `n_pop x n_sims` matrix per N value). Required for
  `plot(result, type = "instability")` and
  `plot(result, type = "calibration")`. Memory use: approximately
  `n_pop * n_sims * 8` bytes per N value (~750 MB for n_pop = 100,000
  and n_sims = 1,000). Default `FALSE`.

- seed:

  Integer or NULL. Seed for reproducibility (see *Reproducibility* in
  Details). Default `NULL`, which uses the session's random number
  stream.

- verbose:

  Logical. Print progress messages per N value. Default `TRUE`.

## Value

A `pmgen_result` object. Use
[`summary()`](https://rdrr.io/r/base/summary.html),
[`print()`](https://rdrr.io/r/base/print.html), and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods to
inspect results, and
[`pmgen_min_n()`](https://joieensor.github.io/pmsampsize/reference/pmgen_min_n.md)
to identify the minimum sample size meeting specified criteria. Element
`fit_info` is a data frame giving, for each N, the number of failed fits
(`n_failed`) and of fits in which predictors were dropped for perfect
prediction (`n_dropped`); `dropped_predictors` is a matrix (one row per
N) of how often each predictor was dropped.

## Details

**Parallel execution:** If the furrr package is installed, the
simulations are run with
[`furrr::future_map()`](https://furrr.futureverse.org/reference/future_map.html),
so they follow whatever future plan is set. Set a parallel backend
before calling `pmgen_simulate()` for speed:

    future::plan(future::multisession, workers = 4)
    result <- pmgen_simulate(...)
    future::plan(future::sequential)  # reset afterwards

Without a parallel plan, or without furrr, execution is sequential.

**Reproducibility:** The simulations for each value of N use seeds that
depend only on `seed`, that N value and `n_sims`. So, for a given
`seed`:

- results are identical whether run sequentially or in parallel, with
  any future plan and any number of workers;

- results for an N value do not depend on the other values in `N`: N =
  400 gives the same results in `N = c(200, 400)` as in
  `N = c(400, 600)`, so a grid can be extended without re-running the
  values already done;

- increasing `n_sims` keeps the earlier simulations: the first 500
  simulations with `n_sims = 1000` are the simulations with
  `n_sims = 500`.

Each simulation uses the session's random number generator settings
([`base::RNGkind()`](https://rdrr.io/r/base/Random.html)), including on
parallel workers. Results for a given `seed` can differ between versions
of R if the generators change. `pmgen_simulate()` changes the session's
random number stream only by calling `set.seed(seed)` (if `seed` is
given) and making one draw.

The `"ridge"`, `"lasso"` and `"elasticnet"` methods need the glmnet
package, and `"firth"` needs logistf.

**Custom methods:** Supply a function of the form

    my_method <- function(train_data, predictors, outcome_col, method_args) {
      # fit some model on train_data
      function(new_data) {
        # return vector of predicted probabilities
      }
    }
    result <- pmgen_simulate(..., method = my_method)

**Development data:** At each simulation iteration, N rows are sampled
with replacement from `ref_model$joint_data`. Binary outcomes are
generated from the reference model's true risks
(`Y ~ Bernoulli(plogis(X %*% beta))`). A warning is issued if N exceeds
the number of rows in `joint_data`.

**Perfect prediction (`method = "logistic"`):** In small development
samples, one level of a binary predictor can have the same outcome for
everyone (e.g. all treated patients have the event). Maximum likelihood
then gives that predictor an effectively infinite coefficient, which
distorts the model's predictions. With
`method_args = list(separation = "drop")` (the default), before each fit
any two-valued predictor with a level whose outcome is constant is
dropped, together with the observations at that level; this is repeated
until no such predictor remains, and the remaining predictors are
fitted. Dropped predictors get a coefficient of 0 when predicting. If a
predictor with more than two values separates the outcome completely,
the fit fails and that simulation's metrics are `NA`. Use
`separation = "keep"` to fit
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) to all predictors and
observations as they are. The penalised methods (`"ridge"`, `"lasso"`,
`"elasticnet"`, `"firth"`) handle separation themselves and ignore this
argument.

**Failed fits:** If fitting or prediction fails in a simulation, its
metrics are `NA` and are left out of the summaries. A warning gives the
number of failed fits for each N and the first error message.

## References

Riley RD, Whittle R, Sadatsafavi M, Martin GP, Pate A, Collins GS, Ensor
J. A general sample size framework for developing or updating a
predictive algorithm: with application to clinical prediction models.
*BMC Medical Research Methodology*. 2026;26:121.
[doi:10.1186/s12874-026-02856-7](https://doi.org/10.1186/s12874-026-02856-7)

## Examples

``` r
# Small settings so the example runs quickly. In practice use the
# pmgen_reference() defaults, n_sims = 1000 or more, and a finer N grid.
ref <- pmgen_reference(
  joint_data = dvt_data, predictors = setdiff(names(dvt_data), "dvt"),
  prevalence = 0.21, cstatistic = 0.78,
  n_bisect_pop = 1e5, n_synth = 1e5, seed = 42
)
#> Running bisection to calibrate reference model...
#> Bisection converged: prevalence=0.2100 (target=0.2100), c-stat=0.7791 (target=0.7800)
#> Unit information matrix computed and cached.
tp <- pmgen_target_population(ref, n = 5000, seed = 42)

result <- pmgen_simulate(ref, tp, N = c(250, 1000), n_sims = 20,
                         method = "logistic",
                         metrics = c("cstat", "cal_slope", "mape"),
                         seed = 1, verbose = FALSE)
result
#> pmgen result
#> ------------
#>   Method              : logistic
#>   N values            : 250, 1000
#>   Metrics             : cstat, cal_slope, mape
#>   Draws/N             : 20
#>   Predictions stored  : no
#> 
#> Summary (means across simulations):
#>     N    metric       mean          sd       q025       q975 degradation_mean
#>   250     cstat 0.77490109 0.011578114 0.75100563 0.79234150     -0.016328173
#>   250 cal_slope 0.74127221 0.152276595 0.50331683 1.07611004     -0.271161746
#>   250      mape 0.06850693 0.014447641 0.04981431 0.09782688      0.068506927
#>  1000     cstat 0.78977302 0.003275223 0.78456655 0.79526384     -0.001456249
#>  1000 cal_slope 0.96798090 0.092676816 0.81243983 1.12431618     -0.044453053
#>  1000      mape 0.03055330 0.007499552 0.01577312 0.04280917      0.030553299
summary(result)
#>      N    metric       mean          sd       q025       q975 degradation_mean
#> 1  250     cstat 0.77490109 0.011578114 0.75100563 0.79234150     -0.016328173
#> 2  250 cal_slope 0.74127221 0.152276595 0.50331683 1.07611004     -0.271161746
#> 3  250      mape 0.06850693 0.014447641 0.04981431 0.09782688      0.068506927
#> 4 1000     cstat 0.78977302 0.003275223 0.78456655 0.79526384     -0.001456249
#> 5 1000 cal_slope 0.96798090 0.092676816 0.81243983 1.12431618     -0.044453053
#> 6 1000      mape 0.03055330 0.007499552 0.01577312 0.04280917      0.030553299

# To run in parallel, set a future plan first (needs furrr and future);
# results are identical to a sequential run with the same seed:
# future::plan(future::multisession, workers = 4)
# result <- pmgen_simulate(ref, tp, N = seq(500, 3000, by = 500),
#                          n_sims = 1000, method = "lasso", seed = 1)
# future::plan(future::sequential)
```
