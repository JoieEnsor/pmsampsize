# Identify the Minimum Sample Size Meeting Specified Criteria

Scans the posterior summaries and raw draws in a `pmgen_result` object
and returns the smallest N at which **all** user-specified criteria are
simultaneously met.

## Usage

``` r
pmgen_min_n(result, criteria = list(), ...)
```

## Arguments

- result:

  A `pmgen_result` object from
  [`pmgen_simulate()`](https://joieensor.github.io/pmsampsize/reference/pmgen_simulate.md).

- criteria:

  Named list of criterion specifications. The name must match a metric
  computed in `result` (e.g. `"cstat"`, `"cal_slope"`, `"mape"`,
  `"net_benefit"`, `"ervsi"`). Each element is itself a named list
  choosing **one** criterion type:

  `range` + `min_prob`

  :   P(metric in `range`) \>= `min_prob`. Use for calibration slope:
      `list(range = c(0.9, 1.1), min_prob = 0.90)`.

  `max_mean`

  :   mean(metric) \<= `max_mean` (for lower-is-better metrics like
      MAPE, or absolute degradation). E.g. `list(max_mean = 0.05)`.

  `min_mean`

  :   mean(metric) \>= `min_mean`. E.g. `list(min_mean = 0.70)` for
      c-statistic floor.

  `max_abs_deg_mean`

  :   \|mean degradation\| \<= threshold. Uses `degradation_mean` from
      summaries. E.g. `list(max_abs_deg_mean = 0.025)`.

  `min_prob_above`

  :   P(metric \>= `threshold`) \>= `min_prob`. E.g. for ERVSI:
      `list(min_prob_above = 0.90, min_prob = 0.90)`.

  `min_prob_below`

  :   P(metric \<= `threshold`) \>= `min_prob`. E.g.
      `list(min_prob_below = 0.05, min_prob = 0.90)` for MAPE.

- ...:

  Reserved for future arguments.

## Value

An object of class `pmgen_min_n`, a list with:

- `min_n`:

  The smallest N meeting all criteria, or `NA` if none do.

- `table`:

  Data frame: one row per N, one column per criterion, plus `all_met`
  (logical) and `N`.

- `criteria`:

  The criteria list as supplied.

## References

Riley RD, Whittle R, Sadatsafavi M, Martin GP, Pate A, Collins GS, Ensor
J. A general sample size framework for developing or updating a
predictive algorithm: with application to clinical prediction models.
*BMC Medical Research Methodology*. 2026;26:121.
[doi:10.1186/s12874-026-02856-7](https://doi.org/10.1186/s12874-026-02856-7)

## Examples

``` r
# Small settings so the example runs quickly; see pmgen_simulate()
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
                         metrics = c("cstat", "cal_slope", "mape"),
                         seed = 1, verbose = FALSE)

pmgen_min_n(result, criteria = list(
  cal_slope = list(range = c(0.8, 1.2), min_prob = 0.80),
  mape      = list(max_mean = 0.05)
))
#> pmgen minimum sample size
#> -------------------------
#>   Minimum N : 1000
#> 
#> Criterion compliance by N:
#>     N cal_slope mape all_met
#>   250       no   no      no 
#>  1000       YES  YES     YES
```
