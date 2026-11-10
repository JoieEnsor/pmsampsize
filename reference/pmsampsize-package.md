# pmsampsize: Sample Size Calculations for Clinical Prediction Models

Computes the minimum sample size required to develop or externally
validate a clinical prediction model. The package has three families of
functions.

## Model development (closed-form criteria)

[`pmsampsize()`](https://joieensor.github.io/pmsampsize/reference/pmsampsize.md)
gives the minimum sample size for developing a model with a continuous,
binary or time-to-event (survival) outcome, using the criteria of Riley
et al. (2019, 2020). The criteria aim to minimise overfitting and to
ensure precise estimation of key parameters in the model. For binary
outcomes, the anticipated Cox-Snell R-squared can be approximated from a
reported C-statistic (Riley et al. 2021a).

## External validation

[`pmvalsampsize()`](https://joieensor.github.io/pmsampsize/reference/pmvalsampsize.md)
gives the minimum sample size for the external validation of an existing
model with a binary outcome, using the criteria of Riley et al. (2021b)
and Archer et al. (2021): precise estimation of the observed/expected
ratio, the calibration slope and the C-statistic, and optionally net
benefit.

## Model development (generalised, simulation-based approach)

The `pmgen_*()` functions implement the general sample size framework of
Riley et al. (2026) for a model with a binary outcome developed with any
fitting method, by simulation:

- [`pmgen_reference()`](https://joieensor.github.io/pmsampsize/reference/pmgen_reference.md)
  specifies a reference (data-generating) model, calibrated to a target
  prevalence and C-statistic, from data or from distributional
  assumptions;

- [`pmgen_target_population()`](https://joieensor.github.io/pmsampsize/reference/pmgen_target_population.md)
  generates the population in which the developed models are evaluated;

- [`pmgen_simulate()`](https://joieensor.github.io/pmsampsize/reference/pmgen_simulate.md)
  simulates development datasets at each candidate sample size, fits the
  model (built in: logistic, ridge, lasso, elastic net, Firth; or a
  user-supplied function) and evaluates its performance;

- [`pmgen_min_n()`](https://joieensor.github.io/pmsampsize/reference/pmgen_min_n.md)
  gives the smallest sample size meeting user-specified criteria, and
  [`pmgen_compare()`](https://joieensor.github.io/pmsampsize/reference/pmgen_compare.md)
  compares two sets of results;

- [summary()](https://joieensor.github.io/pmsampsize/reference/summary.pmgen_result.md)
  and
  [plot()](https://joieensor.github.io/pmsampsize/reference/plot.pmgen_result.md)
  summarise and display the results.

See
[`vignette("pmgen", package = "pmsampsize")`](https://joieensor.github.io/pmsampsize/articles/pmgen.md)
for a worked example.

## Citation

Please cite the package and the methods papers for the functions used;
see `citation("pmsampsize")`.

## References

Riley RD, Snell KIE, Ensor J, Burke DL, Harrell FE Jr, Moons KGM,
Collins GS. Minimum sample size for developing a multivariable
prediction model: Part I - Continuous outcomes. *Statistics in
Medicine*. 2019;38(7):1262-1275.
[doi:10.1002/sim.7993](https://doi.org/10.1002/sim.7993)

Riley RD, Snell KIE, Ensor J, Burke DL, Harrell FE Jr, Moons KGM,
Collins GS. Minimum sample size for developing a multivariable
prediction model: Part II - Binary and time-to-event outcomes.
*Statistics in Medicine*. 2019;38(7):1276-1296.
[doi:10.1002/sim.7992](https://doi.org/10.1002/sim.7992)

Riley RD, Ensor J, Snell KIE, Harrell FE Jr, Martin GP, Reitsma JB,
Moons KGM, Collins G, van Smeden M. Calculating the sample size required
for developing a clinical prediction model. *BMJ*. 2020;368:m441.
[doi:10.1136/bmj.m441](https://doi.org/10.1136/bmj.m441)

Riley RD, Van Calster B, Collins GS. A note on estimating the Cox-Snell
R2 from a reported C statistic (AUROC) to inform sample size
calculations for developing a prediction model with a binary outcome.
*Statistics in Medicine*. 2021;40(4):859-864.
[doi:10.1002/sim.8806](https://doi.org/10.1002/sim.8806) (2021a)

Riley RD, Debray TPA, Collins GS, Archer L, Ensor J, van Smeden M, Snell
KIE. Minimum sample size for external validation of a clinical
prediction model with a binary outcome. *Statistics in Medicine*.
2021;40(19):4230-4251.
[doi:10.1002/sim.9025](https://doi.org/10.1002/sim.9025) (2021b)

Archer L, Snell KIE, Ensor J, Hudda MT, Collins GS, Riley RD. Minimum
sample size for external validation of a clinical prediction model with
a continuous outcome. *Statistics in Medicine*. 2021;40(1):133-146.
[doi:10.1002/sim.8766](https://doi.org/10.1002/sim.8766)

Riley RD, Whittle R, Sadatsafavi M, Martin GP, Pate A, Collins GS, Ensor
J. A general sample size framework for developing or updating a
predictive algorithm: with application to clinical prediction models.
*BMC Medical Research Methodology*. 2026;26:121.
[doi:10.1186/s12874-026-02856-7](https://doi.org/10.1186/s12874-026-02856-7)

## See also

Useful links:

- <https://github.com/JoieEnsor/pmsampsize>

- Report bugs at <https://github.com/JoieEnsor/pmsampsize/issues>

## Author

**Maintainer**: Joie Ensor <j.ensor@bham.ac.uk>
