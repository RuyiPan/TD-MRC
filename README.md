# TD-MRC

Temporal dependence modeling via mixtures of rotated copulas.

This repository contains the R package source for the TD-MRC MCMC sampler plus
separate scripts for simulation and real-data analyses.

## R Package

The package source is in:

```text
MixClayton/mixclayton
```

The GitHub repository is named `TD-MRC`, while the R package name is `TDMRC`
because R package names cannot contain hyphens.

Install from GitHub:

```r
install.packages("remotes")
remotes::install_github("RuyiPan/TD-MRC", subdir = "MixClayton/mixclayton")
library(TDMRC)
```


## Analysis Scripts

Paper-specific scripts are kept outside the package:

```text
MixClayton/analysis/simulation_m2/run_simulation_grid_m2.R
MixClayton/analysis/real_data_m3/run_real_data_grid_m3.R
```

These scripts assume the package has already been installed.

## Other Folders

```text
MixClayton/      Original mixture Clayton scripts and the package source.
SingleClayton/   Single Clayton comparison scripts.
Gaussian/        Gaussian comparison scripts.
```
