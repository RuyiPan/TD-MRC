# Optimized running_time_sim

This folder contains an experimental optimized version of `../running_time_sim.R`.
The original file is not modified.

## Files

- `running_time_sim_optimized.R`: wrapped R implementation with cached index sets.
- `mcmc_rcpp.cpp`: Rcpp helpers for the hot likelihood/allocation loops.

## Main Changes

- Wrapped the script in `run_running_time_sim_optimized()`.
- Moved repeated Clayton rotated log-density calculations into Rcpp.
- Moved latent allocation sampling for `Z` into Rcpp.
- Moved component-wise theta log-likelihood summation into Rcpp.
- Cached `eta_subset` and `inv_subset` values once before MCMC.
- Replaced repeated `log(gamma(...))` with `lgamma(...)` in the eta update.

## Run

From the `MixClayton` folder:

```r
Rscript optimized_running_time_sim/running_time_sim_optimized.R running_time_sim_optimized 1 /path/to/output
```

For a quick interactive smoke test:

```r
source("optimized_running_time_sim/running_time_sim_optimized.R")
res <- run_running_time_sim_optimized(
  job_name = "smoke",
  job_num = 2,
  path = tempdir(),
  TT = 3,
  nt = 20,
  B = 2,
  batch.size = 3,
  verbose = FALSE
)
```

The full run writes:

- `<job_name>_job<job_num>_optimized.rds`
