# Compare Two pmgen_result Objects

Produces a side-by-side data frame comparing the summary statistics from
two `pmgen_result` objects at their shared N values, for example the
results of two
[`pmgen_simulate()`](https://joieensor.github.io/pmsampsize/reference/pmgen_simulate.md)
runs with different methods (such as ridge and lasso).

## Usage

``` r
pmgen_compare(result_1, result_2, stat = "mean")
```

## Arguments

- result_1:

  A `pmgen_result` object.

- result_2:

  A `pmgen_result` object.

- stat:

  Character string. Which summary statistic to compare: `"mean"`,
  `"sd"`, `"q025"`, `"q975"`, or `"degradation_mean"`. Default `"mean"`.

## Value

A data frame with columns: `N`, `metric`, `<method_1>`, `<method_2>`,
`difference` (result_1 minus result_2).

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

# Logistic regression vs Firth's penalised likelihood (needs logistf)
if (requireNamespace("logistf", quietly = TRUE)) {
  result_ml    <- pmgen_simulate(ref, tp, N = 250, n_sims = 10,
                                 method = "logistic", seed = 1,
                                 verbose = FALSE)
  result_firth <- pmgen_simulate(ref, tp, N = 250, n_sims = 10,
                                 method = "firth", seed = 1,
                                 verbose = FALSE)
  pmgen_compare(result_ml, result_firth)
}
#>     N      metric   logistic      firth   difference
#> 1 250   cal_slope 0.71451777 0.77700096 -0.062483192
#> 2 250       cstat 0.77197516 0.77340620 -0.001431043
#> 3 250        mape 0.06847273 0.06422492  0.004247811
#> 4 250 net_benefit 0.01626000 0.01618000  0.000080000
```
