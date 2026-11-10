## Submission

This is a major update (1.1.3 -> 2.0.0). The package now also covers the
external validation of prediction models (`pmvalsampsize()`, previously the
standalone CRAN package pmvalsampsize, by the same maintainer) and a
general, simulation-based sample size framework (`pmgen_*()` functions).
Existing `pmsampsize()` calls are unchanged. Once this version is on CRAN,
pmvalsampsize will be updated to a thin package that re-exports
`pmvalsampsize()` from pmsampsize, so existing code keeps working.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are 2 reverse dependencies on CRAN (predictNMB, R4VN). Both call only
`pmsampsize()`, whose arguments and results are unchanged.
