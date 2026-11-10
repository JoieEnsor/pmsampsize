# Generate a Target Population for Simulation

Partitions the joint predictor distribution stored in a
`pmgen_ref_model` object into two non-overlapping sets: a fixed target
population of size `n` (used to evaluate each simulated model) and a
development pool containing all remaining rows (used to draw development
datasets in
[`pmgen_simulate()`](https://joieensor.github.io/pmsampsize/reference/pmgen_simulate.md)).

## Usage

``` r
pmgen_target_population(ref_model, n = 100000L, dev_pool_n = NULL, seed = NULL)
```

## Arguments

- ref_model:

  A `pmgen_ref_model` object from
  [`pmgen_reference()`](https://joieensor.github.io/pmsampsize/reference/pmgen_reference.md).

- n:

  Integer. Number of individuals in the target population. Must be
  strictly less than `nrow(ref_model$joint_data)`. Default `100000`.

- dev_pool_n:

  Integer or NULL. If `NULL` (default), the development pool contains
  all rows of `joint_data` not selected for the target population. If an
  integer, only that many rows are retained in the pool (sampled from
  the remainder). Useful when memory is limited or when intentionally
  bounding the available development data.

- seed:

  Integer or NULL. Random seed.

## Value

A `pmgen_target_pop` object with slots:

- `$data`: target population data frame (predictors, `lp`, `p_true`,
  `outcome`)

- `$dev_pool`: development pool data frame (predictors only; outcomes
  are generated fresh each simulation draw in
  [`pmgen_simulate()`](https://joieensor.github.io/pmsampsize/reference/pmgen_simulate.md))

- `$ref_model`, `$n`, `$seed`

## Details

The target population rows are sampled **without replacement** from
`joint_data`, so the two sets are guaranteed to be disjoint. The
reference model is applied to the target population rows to compute true
predicted risks and simulated binary outcomes.

## References

Riley RD, Whittle R, Sadatsafavi M, Martin GP, Pate A, Collins GS, Ensor
J. A general sample size framework for developing or updating a
predictive algorithm: with application to clinical prediction models.
*BMC Medical Research Methodology*. 2026;26:121.
[doi:10.1186/s12874-026-02856-7](https://doi.org/10.1186/s12874-026-02856-7)

## Examples

``` r
# Small settings so the example runs quickly; see pmgen_reference()
ref <- pmgen_reference(
  joint_data = dvt_data, predictors = setdiff(names(dvt_data), "dvt"),
  prevalence = 0.21, cstatistic = 0.78,
  n_bisect_pop = 1e5, n_synth = 1e5, seed = 42
)
#> Running bisection to calibrate reference model...
#> Bisection converged: prevalence=0.2100 (target=0.2100), c-stat=0.7791 (target=0.7800)
#> Unit information matrix computed and cached.
tp <- pmgen_target_population(ref, n = 5000, seed = 42)
tp
#> pmgen target population
#> -----------------------
#>   Target pop N    : 5000
#>   Dev pool N      : 4805
#>   Seed            : 42
#>   Prevalence (observed): 0.2088
#>   Mean true risk  : 0.2150
#>   Predictors      : male, OCuse, malignancy, surgery, ab_leg_trauma, vein_distension, calf_diff, ab_ddimer, pregnancy, immobile, clotting
```
