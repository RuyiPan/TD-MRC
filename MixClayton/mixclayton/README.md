# TD-MRC

Reusable R package for MCMC fitting of temporal dependence models via mixtures
of rotated copulas. The package supports `m = 2, 3, 4, 5` dimensions and automatically uses
`2^m` mixture components.

The GitHub repository can be named `TD-MRC`. The R package name is `TDMRC`
because R package names cannot contain hyphens.

This package only contains general MCMC, scoring, simulation, and helper
functions. Paper-specific simulation studies and real-data grid runs are kept
outside the package in `MixClayton/analysis/`.

## Install

From the `TD-MRC` folder:

```sh
R CMD INSTALL MixClayton/mixclayton
```

From RStudio:

```r
setwd("/Users/ruyipan/Desktop/Project/MixCopula/Revision/TD-MRC")
install.packages("MixClayton/mixclayton", repos = NULL, type = "source")
```

After reinstalling in RStudio, restart the R session before running long jobs.

From GitHub, if the repository is `RuyiPan/TD-MRC`:

```r
install.packages("remotes")
remotes::install_github("RuyiPan/TD-MRC", subdir = "MixClayton/mixclayton")
```

## Main Functions

```r
run_mixclayton_mcmc()
run_mixclayton_prediction_mcmc()
run_mixclayton_mdim_mcmc()
run_mixclayton_mdim_prediction_mcmc()
compute_waic()
compute_lpml()
compute_dic()
add_prediction_scores()
simulate_mixclayton()
```

## Example

```r
library(TDMRC)

fit <- run_mixclayton_mdim_prediction_mcmc(
  data = U_train,
  U_test = NULL,
  burn_in = 100,
  B = 200,
  batch.size = 50,
  thin = 10,
  C_tune = 2,
  ek = 1,
  at = 4,
  c0 = 1,
  p = NULL,
  dep_type = list(MA = 4, Season = c(12, 2)),
  seed = 20231213
)

fit$WAIC
fit$LPML
fit$DIC
```

If `U_test = NULL`, LPS is skipped. If `U_test` is supplied, LPS is computed.
