# pmvalsampsize - Sample Size for External Validation of a Prediction Model

`pmvalsampsize` computes the minimum sample size required for the
external validation of an existing multivariable prediction model using
the criteria proposed by Archer *et al*. and Riley *et al*.
`pmvalsampsize` can currently be used to calculate the minimum sample
size for the external validation of models with binary outcomes.  
  
*Continuous and survival (time-to-event) outcome model calculations are
a work in progress*.  
  

## Usage

``` r
pmvalsampsize(
  type,
  cslope = 1,
  csciwidth = 0.2,
  oe = 1,
  oeciwidth = 0.2,
  cstatistic = NA,
  cstatciwidth = 0.1,
  simobs = 1e+06,
  lpnormal = NA,
  lpbeta = NA,
  lpcstat = FALSE,
  tolerance = NULL,
  increment = NULL,
  oeseincrement = 1e-04,
  graph = FALSE,
  trace = NULL,
  prevalence = NA,
  seed = 123456,
  sensitivity = NA,
  specificity = NA,
  threshold = NA,
  nbciwidth = 0.2,
  nbseincrement = 1e-04
)
```

## Arguments

- type:

  specifies the type of analysis for which sample size is being
  calculated

  - `"b"` specifies sample size calculation for a prediction model with
    a binary outcome

- cslope:

  specifies the anticipated c-slope performance in the validation
  sample. Default conservatively assumes perfect c-slope=1. The value
  could alternatively be based on a previous validation study for
  example. For **binary outcomes** the c-slope calculation requires the
  user to specify a distribution for the assumed LP in the validation
  sample (or alternatively the distribution of predicted probabilities
  in the validation sample). See `lp*()` options below.

- csciwidth:

  specifies the target CI width (acceptable precision) for the c-slope
  performance. Default assumes CI width=0.2.

- oe:

  specifies the anticipated O/E performance in the validation sample.
  Default conservatively assumes perfect O/E=1.

- oeciwidth:

  specifies the target CI width (acceptable precision) for the E/O
  performance. Default assumes CI width=0.2. The choice of CI width is
  context specific, and depends on the event probability in the
  population. See Riley et al. for further details.

- cstatistic:

  specifies the anticipated c-statistic performance in the validation
  sample. This is a **required input**. May be based on the
  optimism-adjusted c-statistic reported in the development study for
  the existing prediction model. Ideally, this should be an
  optimism-adjusted c-statistic. NB: This input is also used when using
  the `lpcstat()` option.

- cstatciwidth:

  specifies the target CI width (acceptable precision) for the
  c-statistic performance. Default assumes CI width=0.1.

- simobs:

  specifies the number of observations to use when simulating the LP
  distribution for c-slope calculation in criteria 2. Default
  observations=1,000,000. Higher `simobs()` values will reduce random
  variation further.

- lpnormal:

  defines parameters to simulate the LP distribution for criteria 2 from
  a normal distribution. The user must specify the mean and standard
  deviation (in this order) of the LP distribution.

- lpbeta:

  defines parameters to simulate the distribution of predicted
  probabilities for criteria 2 from a beta distribution. The user must
  specify the alpha and beta parameters (in this order) of the
  probability distribution. The LP distribution is then generated
  internally using this probability distribution.

- lpcstat:

  logical. If `TRUE`, the LP distribution for criteria 2 is simulated
  assuming that the distributions of the LP in events and non-events are
  normal with a common variance, derived from the anticipated
  `cstatistic` (which must be greater than 0.5). The mean for the
  non-events distribution is calculated so that the expected outcome
  proportion matches the anticipated `prevalence` in the validation
  sample, and the mean for the events distribution follows from this.
  Default `FALSE`. NB: this approach makes strong assumptions of
  normality and equal variances in each outcome group, which may be
  unrealistic in most situations. In pmvalsampsize 0.1.0 `lpcstat` took
  a starting value for the non-events mean; a numeric value is still
  accepted, with a warning, and treated as `TRUE`.

- tolerance, increment, trace:

  Retired in pmsampsize 2.0.0, and ignored with a warning if supplied.
  These controlled the iterative search for the non-events mean used by
  `lpcstat` in pmvalsampsize 0.1.0, which is now calculated exactly.

- oeseincrement:

  sets the increment by which to iterate when identifying the SE(ln(OE))
  to meet the target CI width specified for OE. The default iteration
  increment=0.0001. In the majority of cases this will be suitably small
  to ensure a precise SE is identified. The user should check the output
  table to ensure that the target CI width has been attained and adjust
  the increment if necessary.

- graph:

  specifies that a histogram of the simulated LP distribution for
  criteria 2 is produced. The graph also details summary statistics for
  the simulated distribution. Useful option for checking the simulated
  LP distribution against the source of input parameters. Also useful
  for reporting at publication.

- prevalence:

  specifies the overall outcome proportion (for a prognostic model) or
  overall prevalence (for a diagnostic model) expected within the model
  validation sample. This is a **required input**. This should be
  derived based on previous studies in the same population or directly
  from the validation sample if to hand.

- seed:

  specifies the initial value of the random-number seed used by the
  random-number functions when simulating data to approximate the LP
  distribution for criteria 2.

- sensitivity:

  specifies the anticipated sensitivity performance in the validation
  sample at the chosen risk threshold (specified using `threshold`). If
  sensitivity and specificity are not provided then `pmvalsampsize` uses
  the simulated LP distribution from criteria 2 and the user-specified
  risk threshold to estimate the anticipated sensitivity and specificity
  to be used in calculation of net benefit. NB: net benefit criteria is
  not calculated if either i) `sensitivity`, `specificity` and
  `threshold` or ii) `threshold` option are not specified.

- specificity:

  specifies the anticipated specificity performance in the validation
  sample at the chosen risk threshold (specified using `threshold`). If
  sensitivity and specificity are not provided then `pmvalsampsize` uses
  the simulated LP distribution from criteria 2 and the user-specified
  risk threshold to estimate the anticipated sensitivity and specificity
  to be used in calculation of net benefit. NB: net benefit criteria is
  not calculated if either i) `sensitivity`, `specificity` and
  `threshold` or ii) `threshold` option are not specified.

- threshold:

  specifies the risk threshold to be used for calculation of net benefit
  performance of the model in the validation sample. If sensitivity and
  specificity are not provided then `threshold` must be given in order
  for `pmvalsampsize` to assess sample size requirements for net
  benefit. NB: net benefit criteria is not calculated if either i)
  `sensitivity`, `specificity` and `threshold` or ii) `threshold` option
  are not specified.

- nbciwidth:

  specifies the target CI width (acceptable precision) for the
  standardised net benefit performance. Default assumes CI width=0.2.
  The choice of CI width is context specific. See Riley et al. for
  further details.

- nbseincrement:

  sets the increment by which to iterate when identifying the
  SE(standardised net benefit) to meet the target CI width specified for
  standardised net benefit. The default iteration increment=0.0001. In
  the majority of cases this will be suitably small to ensure a precise
  SE is identified. The user should check the output table to ensure
  that the target CI width has been attained and adjust the increment if
  necessary.

## Value

A list including a matrix of calculated sample size requirements for
each criteria defined under *'Details'*, and a series of vectors of
parameters used in the calculations as well as the the final recommended
minimum sample size and number of events required for external
validation.

## Details

A series of criteria define the sample size needed to ensure precise
estimation of key measures of prediction model performance, allowing
conclusions to be drawn about whether the model is potentially accurate
and useful in a given population of interest.  
  

For **binary outcomes**, there are three criteria to calculate the
sample size (N) needed for:  
i) precise estimation of the Observed/Expected (O/E) statistic,  
ii) precise estimation of the calibration slope (c-slope), and  
iii) precise estimation of the c-statistic.

The sample size calculation requires the user to pre-specify the
following:

- the outcome event proportion

- the anticipated O/E performance

- the target precision (CI width) for the O/E

- the anticipated c-slope

- the target precision (CI width) for the c-slope

- the anticipated c-statistic performance

- the target precision (CI width) for the c-statistic

- the distribution of estimated probabilities from the model, ideally
  specified on the log-odds scale - AKA the Linear Predictor (LP)

With thanks to Richard Riley for helpful feedback

## References

Archer L, Snell K, Ensor J, Hudda MT, Collins GS, Riley RD. Minimum
sample size for external validation of a clinical prediction model with
a continuous outcome. *Statistics in Medicine*. 2021;40(1):133-146.
[doi:10.1002/sim.8766](https://doi.org/10.1002/sim.8766)

Riley RD, Debray TPA, Collins G, Archer L, Ensor J, van Smeden M, Snell
KIE. Minimum sample size for external validation of a clinical
prediction model with a binary outcome. *Statistics in Medicine*.
2021;40(19):4230-4251.
[doi:10.1002/sim.9025](https://doi.org/10.1002/sim.9025)

## See also

[`pmsampsize()`](https://joieensor.github.io/pmsampsize/reference/pmsampsize.md)
and the `pmgen_*()` functions (e.g.
[`pmgen_simulate()`](https://joieensor.github.io/pmsampsize/reference/pmgen_simulate.md))
for model development.

## Author

Joie Ensor (University of Birmingham, j.ensor@bham.ac.uk)

## Examples

``` r
# Examples based on those included in the papers referenced below by
# Riley et al. & Archer et al. published in Statistics in Medicine.

# Note that the examples below use a very low simulation sample for criteria
# 2 for brevity. Sizes matching the default or higher are recommended to 
# minimise uncertainty in the calculated sample size requirements.

# Binary outcomes (Logistic prediction models)
# Use pmvalsampsize to calculate the minimum sample size required to
# externally validate an existing multivariable prediction model for a
# binary outcome (e.g. mechanical heart valve failure).  Based on previous
# evidence, the outcome prevalence is anticipated to be 0.018 (1.8%) and the
# reported c-statistic was 0.8.  The LP distribution was published and
# appeared normally distributed with mean(SD) of -5(2.5).  We target default
# CI widths for all but O/E CI width=1 (see Riley et al. for details). We can
# use the graph option to check the simulated distribution is appropriate.

pmvalsampsize(type = "b", prevalence = 0.018, cstatistic = 0.8,
lpnorm = c(-5,2.5), oeciwidth = 1, simobs = 1000)
#> Normal LP distribution with parameters - mean=-5, sd=2.5
#>                          Samp_size Perf    SE CI width
#> Criteria 1 - O/E               905  1.0 0.246      1.0
#> Criteria 2 - C-slope          4501  1.0 0.051      0.2
#> Criteria 3 - C statistic      4252  0.8 0.026      0.1
#> Final SS                      4501  1.0 0.051      0.2
#>  
#>  Minimum sample size required for model validation based on user inputs = 4501, 
#>  with 82 events (assuming an outcome prevalence = 0.018) 
#>  
#>  


# Alternatively, lets assume that the authors provided a distribution of
# predicted probabilities (e.g. as part of a calibration plot).  We can use
# this to specify parameters for a beta distribution to simulate the LP
# distribution as below:
pmvalsampsize(type = "b", prevalence = 0.018, cstatistic = 0.8,
lpbeta = c(0.5,0.5), oeciwidth = 1, simobs = 1000)
#> Beta P distribution with parameters - alpha=0.5, beta=0.5
#>                          Samp_size Perf    SE CI width
#> Criteria 1 - O/E               905  1.0 0.246      1.0
#> Criteria 2 - C-slope          1641  1.0 0.051      0.2
#> Criteria 3 - C statistic      4252  0.8 0.026      0.1
#> Final SS                      4252  0.8 0.026      0.1
#>  
#>  Minimum sample size required for model validation based on user inputs = 4252, 
#>  with 77 events (assuming an outcome prevalence = 0.018) 
#>  
#>  

# Finally, we can use the anticipated c-statistic to simulate the event and
# non-event distributions assuming normality and common variances.  The mean
# of the non-events distribution is calculated to match the prevalence:

pmvalsampsize(type = "b", prevalence = 0.018, cstatistic = 0.8,
lpcstat = TRUE, oeciwidth = 1, seed = 1234, simobs = 10000)
#> LP distribution from c-statistic (normal events & non-events, common variance) - mean in non-event group=-4.70755
#> Proportion of outcome events under simulation=0.0183, target prevalence=0.018
#>                          Samp_size Perf    SE CI width
#> Criteria 1 - O/E               905  1.0 0.246      1.0
#> Criteria 2 - C-slope         16270  1.0 0.051      0.2
#> Criteria 3 - C statistic      4252  0.8 0.026      0.1
#> Final SS                     16270  1.0 0.051      0.2
#>  
#>  Minimum sample size required for model validation based on user inputs = 16270, 
#>  with 293 events (assuming an outcome prevalence = 0.018) 
#>  
#>  


```
