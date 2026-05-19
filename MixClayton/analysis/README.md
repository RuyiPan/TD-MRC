# Analysis Scripts

This folder contains paper-specific scripts that use the `TDMRC`
package. These scripts are intentionally kept outside the R package so the
package remains a reusable MCMC implementation.

## Folder Layout

```text
analysis/
  simulation_m2/
    run_simulation_grid_m2.R
  real_data_m3/
    run_real_data_grid_m3.R
```

## Install Package First

Run this from the `TD-MRC` folder:

```r
setwd("/Users/ruyipan/Desktop/Project/MixCopula/Revision/TD-MRC")
install.packages("MixClayton/mixclayton", repos = NULL, type = "source")
```

Restart the R session in RStudio after reinstalling.

## m = 2 Simulation Grid

```r
setwd("/Users/ruyipan/Desktop/Project/MixCopula/Revision/TD-MRC")
source("MixClayton/analysis/simulation_m2/run_simulation_grid_m2.R")
```

## m = 3 Real-Data Grid

This uses all time points in `MCMC_MRC/Data/data_2017_2019.rds` for training.
LPS is skipped because there is no held-out test set.

```r
setwd("/Users/ruyipan/Desktop/Project/MixCopula/Revision/TD-MRC")
source("MixClayton/analysis/real_data_m3/run_real_data_grid_m3.R")
```
