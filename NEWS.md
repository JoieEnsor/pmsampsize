pmsampsize v2.0.0 
=============================================
pmsampsize now covers the development and external validation of clinical prediction models, in three families of functions: `pmsampsize()` (development), `pmvalsampsize()` (external validation) and the new `pmgen_*()` functions (generalised, simulation-based approach). **Existing `pmsampsize()` calls are unchanged** (same arguments and results). Existing `pmvalsampsize()` calls give the same results as pmvalsampsize 0.1.0, apart from the `lpcstat` option and the type of `simulated_data_cstat` (see below).

* New package title ("Sample Size Calculations for Clinical Prediction Models"), a package help page (`?pmsampsize-package`) and a `CITATION` file: `citation("pmsampsize")` lists the package and the methods papers for each function.
* References in the help pages updated with full details and DOIs.
* `pmvalsampsize()` (sample size for external validation of a prediction model; previously the standalone pmvalsampsize package, v0.1.0, released 2023-11-14) is now part of pmsampsize. Existing `pmvalsampsize()` calls give the same results.
* `pmvalsampsize()` no longer depends on pROC. The c-statistic of the simulated data for the `lpcstat` option is now computed internally (same value), and `simulated_data_cstat` in the returned object is now a plain number rather than a pROC `auc` object.
* `pmvalsampsize()` `lpcstat` option rewritten. The mean of the non-events LP distribution is now calculated exactly from `prevalence` and `cstatistic`, instead of by an iterative simulation search from a user-supplied start value. This fixes cases where the search never finished (e.g. `prevalence = 0.2, cstatistic = 0.75`). Results for `lpcstat` change slightly because the mean is now more precise (e.g. the help-file example: criteria 2 sample size 16,175 -> 16,270).
* `lpcstat` is now logical (`lpcstat = TRUE` selects the option; default `FALSE`) and requires `cstatistic` greater than 0.5. A numeric value, as used in 0.1.0, is still accepted and treated as `TRUE`, with a warning.
* `pmvalsampsize()` arguments `tolerance`, `increment` and `trace` are retired (they only controlled the old `lpcstat` search). They are ignored, with a warning, if supplied, and will be removed in a future major version.
* The `lpcstat` iteration messages are replaced by a single message giving the calculated non-events mean and the simulated vs target prevalence.
* New `pmgen_*()` functions for the general, simulation-based sample size framework of Riley et al. (2026, BMC Medical Research Methodology 26:121) for developing a clinical prediction model (binary outcomes), for any model-fitting method:
  * `pmgen_reference()` specifies a reference model, calibrated to a target prevalence and c-statistic, from data or from distributional assumptions.
  * `pmgen_target_population()` generates the target population on which the developed models are evaluated.
  * `pmgen_simulate()` simulates development datasets at each candidate sample size, fits the model (built in: logistic, ridge, lasso, elastic net, Firth; or a user-supplied fitting function) and evaluates its performance (c-statistic, calibration slope and intercept, O:E ratio, MAPE, RMSPE, Brier score, net benefit, ERVSI).
  * `pmgen_min_n()` gives the smallest sample size meeting user-specified criteria; `pmgen_compare()` compares two sets of results.
  * `summary()`, `print()` and `plot()` methods for the results.
  * Reproducible with `seed`: each N value's simulations use seeds that depend only on `seed`, that N value and `n_sims`. Results are therefore identical sequentially or in parallel (any future plan, any number of workers, using the session's `RNGkind()` on the workers too), an N value gives the same results whatever other values are in `N`, and increasing `n_sims` keeps the earlier simulations.
  * `dvt_data`: an example dataset (simulated DVT data, 9,805 patients, 11 binary predictors).
* New dependencies: MASS, mvtnorm and rlang (Imports). Packages needed only for particular features are suggested, not required: glmnet (ridge, lasso, elastic net), logistf (Firth), ggplot2 and scales (plots), furrr and future (parallel simulation). Results are the same with or without furrr.
* pmsampsize now requires R >= 4.0.0.

pmsampsize v1.1.3 (Release date: 2023-12-05)
=============================================
* added options to specify the cox-snell or nagelkerke's r-squared directly for binary outcome models
* simplified continuous outcome models criteria 3 calculation fixing MMOE=10% so that N=234+p

pmsampsize v1.1.2 (Release date: 2022-02-12)
=============================================
* fixes issue for continuous outcome models in criteria 4 in calculating the CI for the intercept based on final SS
* edit to criteria 2 for continuous outcomes to round up N

pmsampsize v1.1.1 (Release date: 2021-11-21)
=============================================
* fixes issues with rounding leading to small differences in Max_Rsq and subsequently sample size calculations
* additional output including Nagelkerke's R-sq (binary & survival outcomes)

pmsampsize v1.1.0 (Release date: 2021-06-24)
=============================================
* update allowing user to input C-statistic from existing model to approximate Cox-Snell R-squared (for binary outcome model sample size calculations) 
* additional parameter checks	

pmsampsize v1.0.4 (Release date: 2020-07-24)
=============================================
* version issues fixes

pmsampsize v1.0.3 (Release date: 2020-07-23)
=============================================

Continuous outcomes:

* calculation of shrinkage for each criteria based on N needed in specific criteria

* also allowance for the unique scenario when starting N already meets shrinkage requirement set by user (this is more for academic exercise, as in practice a high shrinkage will be desired and so it is very rare that the base N (set to p+2) will provide good shrinkage)

Binary and TTE outcomes:

* more informative error message added for rare scenario where shrinkage is specified lower than r-squared which results in NaNs from shrinkage formula 

pmsampsize v1.0.2 (Release date: 2019-11-15)
=============================================

* update forces the time-unit to be the mean follow-up time, as described in forthcoming BMJ paper which extends the original Statistics in Medicine paper (see reference in help file)

* updated in calculation of max R-sq for survival option, to prevent manipulation based on time units used 


pmsampsize v1.0.1 (Release date: 2019-08-02)
=============================================

* updated in line with correction to original article in eq. 13


pmsampsize v1.0.0 (Release date: 2019-01-02)
=============================================

* this is a new release

