# TD-MRC

Reusable R package for MCMC fitting of temporal dependence models via mixtures
of rotated copulas.

The GitHub repository can be named `TD-MRC`. The R package name is `TDMRC`
because R package names cannot contain hyphens.

This package only contains general MCMC, scoring, simulation, and helper
functions. Paper-specific simulation studies and real-data grid runs are kept
outside the package in `MixClayton/analysis/`.

The fitted dimension `m` is inferred from the number of columns in the data;
the mixture uses `2^m` rotated components, subject only to practical memory and
runtime limits.

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

## Example: fit an m = 2 simulated data set

```r
library(TDMRC)

# Simulate a small m = 2 data set.
# The data object is a list, with one matrix for each time point.
set.seed(20231213)
U_train <- simulate_mixclayton(
  TT = 20,
  nt = 300,
  theta = c(5, 3, 4, 3),
  weights0 = c(0.4, 0.25, 0.1, 0.25)
)

# For m = 2, the model has 2^m = 4 mixture components.
fit <- run_mixclayton_mdim_prediction_mcmc(
  data = U_train,
  U_test = NULL,
  burn_in = 60,
  B = 140,
  batch.size = 50,
  thin = 10,
  C_tune = 2,
  ek = 1,
  at = 4,
  c0 = 1,
  p = NULL,
  dep_type = list(MA = 2, Season = c(12, 0)),
  seed = 20231213,
  verbose = FALSE
)

fit$WAIC
fit$LPML
fit$DIC
```

Increase `burn_in` and `B` for the final analysis. If `U_test = NULL`, LPS is
skipped. If `U_test` is supplied, LPS is computed.
