#' pmsampsize: Sample Size Calculations for Clinical Prediction Models
#'
#' Computes the minimum sample size required to develop or externally
#' validate a clinical prediction model. The package has three families of
#' functions.
#'
#' @section Model development (closed-form criteria):
#' [pmsampsize()] gives the minimum sample size for developing a model with a
#' continuous, binary or time-to-event (survival) outcome, using the criteria
#' of Riley et al. (2019, 2020). The criteria aim to minimise overfitting and
#' to ensure precise estimation of key parameters in the model. For binary
#' outcomes, the anticipated Cox-Snell R-squared can be approximated from a
#' reported C-statistic (Riley et al. 2021a).
#'
#' @section External validation:
#' [pmvalsampsize()] gives the minimum sample size for the external validation
#' of an existing model with a binary outcome, using the criteria of Riley et
#' al. (2021b) and Archer et al. (2021): precise estimation of the
#' observed/expected ratio, the calibration slope and the C-statistic, and
#' optionally net benefit.
#'
#' @section Model development (generalised, simulation-based approach):
#' The `pmgen_*()` functions implement the general sample size framework of
#' Riley et al. (2026) for a model with a binary outcome developed with any
#' fitting method, by simulation:
#' * [pmgen_reference()] specifies a reference (data-generating) model,
#'   calibrated to a target prevalence and C-statistic, from data or from
#'   distributional assumptions;
#' * [pmgen_target_population()] generates the population in which the
#'   developed models are evaluated;
#' * [pmgen_simulate()] simulates development datasets at each candidate
#'   sample size, fits the model (built in: logistic, ridge, lasso, elastic
#'   net, Firth; or a user-supplied function) and evaluates its performance;
#' * [pmgen_min_n()] gives the smallest sample size meeting user-specified
#'   criteria, and [pmgen_compare()] compares two sets of results;
#' * [summary()][summary.pmgen_result()] and [plot()][plot.pmgen_result()]
#'   summarise and display the results.
#'
#' See `vignette("pmgen", package = "pmsampsize")` for a worked example.
#'
#' @section Citation:
#' Please cite the package and the methods papers for the functions used; see
#' `citation("pmsampsize")`.
#'
#' @references
#' Riley RD, Snell KIE, Ensor J, Burke DL, Harrell FE Jr, Moons KGM,
#' Collins GS. Minimum sample size for developing a multivariable prediction
#' model: Part I - Continuous outcomes. *Statistics in Medicine*.
#' 2019;38(7):1262-1275. \doi{10.1002/sim.7993}
#'
#' Riley RD, Snell KIE, Ensor J, Burke DL, Harrell FE Jr, Moons KGM,
#' Collins GS. Minimum sample size for developing a multivariable prediction
#' model: Part II - Binary and time-to-event outcomes. *Statistics in
#' Medicine*. 2019;38(7):1276-1296. \doi{10.1002/sim.7992}
#'
#' Riley RD, Ensor J, Snell KIE, Harrell FE Jr, Martin GP, Reitsma JB,
#' Moons KGM, Collins G, van Smeden M. Calculating the sample size required
#' for developing a clinical prediction model. *BMJ*. 2020;368:m441.
#' \doi{10.1136/bmj.m441}
#'
#' Riley RD, Van Calster B, Collins GS. A note on estimating the Cox-Snell
#' R2 from a reported C statistic (AUROC) to inform sample size calculations
#' for developing a prediction model with a binary outcome. *Statistics in
#' Medicine*. 2021;40(4):859-864. \doi{10.1002/sim.8806} (2021a)
#'
#' Riley RD, Debray TPA, Collins GS, Archer L, Ensor J, van Smeden M,
#' Snell KIE. Minimum sample size for external validation of a clinical
#' prediction model with a binary outcome. *Statistics in Medicine*.
#' 2021;40(19):4230-4251. \doi{10.1002/sim.9025} (2021b)
#'
#' Archer L, Snell KIE, Ensor J, Hudda MT, Collins GS, Riley RD. Minimum
#' sample size for external validation of a clinical prediction model with a
#' continuous outcome. *Statistics in Medicine*. 2021;40(1):133-146.
#' \doi{10.1002/sim.8766}
#'
#' Riley RD, Whittle R, Sadatsafavi M, Martin GP, Pate A, Collins GS, Ensor J. A general sample size framework
#' for developing or updating a predictive algorithm: with application to
#' clinical prediction models. *BMC Medical Research Methodology*.
#' 2026;26:121. \doi{10.1186/s12874-026-02856-7}
#'
#' @name pmsampsize-package
#' @aliases pmsampsize-package
#' @keywords internal
"_PACKAGE"
