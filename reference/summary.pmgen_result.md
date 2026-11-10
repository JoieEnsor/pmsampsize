# Summarise a pmgen_result object

Returns a formatted data frame of posterior summary statistics across
all evaluated N values, one row per metric per N. Columns: N, metric,
mean, sd, q025 (2.5th percentile), q975 (97.5th percentile),
degradation_mean.

## Usage

``` r
# S3 method for class 'pmgen_result'
summary(object, ...)
```

## Arguments

- object:

  A `pmgen_result` object from
  [`pmgen_simulate()`](https://joieensor.github.io/pmsampsize/reference/pmgen_simulate.md).

- ...:

  Ignored.

## Value

A data frame with one row per (N, metric) combination.
