# Plot a pmgen_result object

Plots the results of
[`pmgen_simulate()`](https://joieensor.github.io/pmsampsize/reference/pmgen_simulate.md).
Requires the ggplot2 and scales packages.

## Usage

``` r
# S3 method for class 'pmgen_result'
plot(x, type = "degradation", N = NULL, ...)
```

## Arguments

- x:

  A `pmgen_result` object.

- type:

  Character. One of `"degradation"`, `"assurance"`, `"instability"`,
  `"calibration"`. The `"instability"` and `"calibration"` plots need
  the predictions stored (`store_predictions = TRUE` in
  [`pmgen_simulate()`](https://joieensor.github.io/pmsampsize/reference/pmgen_simulate.md)).

- N:

  Integer. For `"instability"` and `"calibration"`, which N value to
  plot. Defaults to the first N in the result.

- ...:

  Additional arguments passed to the specific plot function.

## Value

A ggplot object.
