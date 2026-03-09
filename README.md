# pmsampsize 

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/pmsampsize)](https://CRAN.R-project.org/package=pmsampsize)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/pmsampsize)](https://CRAN.R-project.org/package=pmsampsize)
[![R-CMD-check](https://github.com/JoieEnsor/pmsampsize/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JoieEnsor/pmsampsize/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->


## Overview

`pmsampsize` computes the **minimum sample size required for the development of a new multivariable prediction model**, using the criteria proposed by Riley *et al.* (2018, 2020). It supports models with:

- **Continuous** outcomes (linear prediction models)
- **Binary** outcomes (logistic prediction models)
- **Survival / time-to-event** outcomes (Cox prediction models)

The criteria implemented aim to minimise overfitting and ensure precise estimation of key parameters in the prediction model.

## Installation

Install the released version from CRAN:

```r
install.packages("pmsampsize")
```

Or install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("JoieEnsor/pmsampsize")
```

## Usage

### Binary outcome

Calculate the minimum sample size for a logistic prediction model with 24 candidate predictor parameters, an anticipated outcome prevalence of 17.4%, and a Cox-Snell R² of 0.288 (taken from an existing model):

```r
library(pmsampsize)

pmsampsize(type = "b", csrsquared = 0.288, parameters = 24, prevalence = 0.174)
```

If only a C-statistic (rather than R²) is available from an existing model, use the `cstatistic` argument instead:

```r
pmsampsize(type = "b", cstatistic = 0.89, parameters = 24, prevalence = 0.174)
```

### Survival outcome

Calculate the minimum sample size for a Cox prediction model with 30 candidate predictors, using information from an existing study (R² = 0.051, mean follow-up = 2.07 years, event rate = 0.065), with a prediction timepoint of 2 years:

```r
pmsampsize(type = "s", csrsquared = 0.051, parameters = 30, rate = 0.065,
           timepoint = 2, meanfup = 2.07)
```

### Continuous outcome

Calculate the minimum sample size for a linear prediction model with 25 candidate predictors, R² = 0.2, mean outcome value of 1.9, and SD of 0.6:

```r
pmsampsize(type = "c", rsquared = 0.2, parameters = 25, intercept = 1.9, sd = 0.6)
```

## Criteria implemented

### Continuous outcomes — four criteria:
1. Small overfitting: expected shrinkage of predictor effects ≤ 10%
2. Small absolute difference (≤ 0.05) between apparent and adjusted R²
3. Precise estimation of the residual standard deviation
4. Precise estimation of the average outcome value

### Binary and survival outcomes — three criteria:
1. Small overfitting: expected shrinkage of predictor effects ≤ 10%
2. Small absolute difference (≤ 0.05) between apparent and adjusted Nagelkerke's R²
3. Precise estimation (within ± 0.05) of the average outcome risk at the timepoint of interest

## Key arguments

| Argument | Type(s) | Description |
|---|---|---|
| `type` | all | `"c"` continuous, `"b"` binary, `"s"` survival |
| `parameters` | all | Number of candidate predictor parameters |
| `shrinkage` | all | Target shrinkage (default = 0.9) |
| `rsquared` | `"c"` | Anticipated R² of new model |
| `csrsquared` | `"b"`, `"s"` | Anticipated Cox-Snell R² |
| `nagrsquared` | `"b"`, `"s"` | Anticipated Nagelkerke's R² |
| `cstatistic` | `"b"` | C-statistic from existing model (used to approximate Cox-Snell R²) |
| `prevalence` | `"b"` | Anticipated outcome proportion |
| `rate` | `"s"` | Overall event rate in the population |
| `timepoint` | `"s"` | Timepoint of interest for prediction |
| `meanfup` | `"s"` | Anticipated mean follow-up time |
| `intercept` | `"c"` | Average outcome value in the population |
| `sd` | `"c"` | Standard deviation of outcome values |
| `mmoe` | `"c"` | Multiplicative margin of error for intercept (default = 1.1) |

## Citation

If you use `pmsampsize` in your research, please cite the following papers:

Riley RD, Ensor J, Snell KIE, Harrell FE, Martin GP, Reitsma JB, et al. Calculating the sample size required for developing a clinical prediction model. *BMJ*. 2020. <https://doi.org/10.1136/bmj.m441>

Riley RD, Snell KIE, Ensor J, Burke DL, Harrell FE Jr., Moons KG, Collins GS. Minimum sample size required for developing a multivariable prediction model: Part I — continuous outcomes. *Statistics in Medicine*. 2018. <https://doi.org/10.1002/sim.7993>

Riley RD, Snell KIE, Ensor J, Burke DL, Harrell FE Jr., Moons KG, Collins GS. Minimum sample size required for developing a multivariable prediction model: Part II — binary and time-to-event outcomes. *Statistics in Medicine*. 2018. <https://doi.org/10.1002/sim.7992>

Riley RD, Van Calster B, Collins GS. A note on estimating the Cox-Snell R² from a reported C statistic (AUROC) to inform sample size calculations for developing a prediction model with a binary outcome. *Statistics in Medicine*. 2020.

## Author

**Joie Ensor** (University of Birmingham, [j.ensor@bham.ac.uk](mailto:j.ensor@bham.ac.uk))

With thanks to Richard D. Riley, Emma C. Martin, Gary Collins, Glen Martin & Kym Snell for helpful input and feedback.

## License

GPL (≥ 3)
