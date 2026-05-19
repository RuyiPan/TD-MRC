# Run the m = 3 real-data grid for LPML and WAIC tables.
#
# This follows the settings in MCMC_MRC/real_dim3_prediction.R:
#   data: Data/data_2017_2019.rds
#   burn_in = 100, B = 200, batch.size = 50, thin = 10
#   ek = 1, c0 = 1, p uniform over 2^3 = 8 components
#
# In RStudio:
#   setwd("/Users/ruyipan/Desktop/Project/MixCopula/Revision/TD-MRC")
#   source("MixClayton/analysis/real_data_m3/run_real_data_grid_m3.R")

suppressPackageStartupMessages({
  library(TDMRC)
})

args <- commandArgs(TRUE)

data_path <- if (length(args) >= 1) args[1] else "MCMC_MRC/Data/data_2017_2019.rds"
output_dir <- if (length(args) >= 2) args[2] else "MixClayton/analysis/real_data_m3/results"

# Use all available real-data time points for training.
use_all_data_for_train <- TRUE

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

burn_in <- 100
B <- 200
batch.size <- 50
thin <- 10
C_tune <- 2

ek <- 1
c0 <- 1
seed <- 20231213
season_lag <- 12

at_grid <- 0:5
MA_grid <- 0:4
Season_grid <- 0:2

# Match the manuscript table: when at = 0, only MA(0) is fitted because the
# dynamic prior contribution is zero and MA does not change the model.
settings <- do.call(
  rbind,
  lapply(Season_grid, function(Season) {
    rbind(
      data.frame(at = 0, MA = 0, Season = Season),
      expand.grid(at = 1:5, MA = MA_grid, Season = Season)
    )
  })
)
settings$setting_id <- seq_len(nrow(settings))

message("Package loaded from: ", find.package("TDMRC"))
message("Package version: ", as.character(packageVersion("TDMRC")))

U_all <- readRDS(data_path)
if (!is.list(U_all)) stop("Expected `data_path` to contain a list of matrices.")
TT_train <- length(U_all)
U_train <- U_all
U_test <- NULL

m <- ncol(U_train[[1]])
K <- 2^m
if (m != 3) stop("Expected m = 3 real data, but ncol(U_train[[1]]) is ", m, ".")

density_component_count <- ncol(TDMRC:::clayton_rot_logdens_cpp(
  U_train[[1]][1:2, , drop = FALSE],
  rep(1, K)
))
message("C++ density returns ", density_component_count, " components.")
if (density_component_count != K) {
  stop(
    "The loaded package is stale: expected ", K, " components but C++ returned ",
    density_component_count, ". Reinstall the package, restart RStudio's R session, ",
    "and run this script again."
  )
}

message("Loaded ", length(U_all), " time points from: ", normalizePath(data_path))
message("Using all data for training: TT_train = ", TT_train)
message("Using m = ", m, " dimensions and K = ", K, " mixture components.")
message("Total grid settings: ", nrow(settings))

saveRDS(
  list(
    data_path = data_path,
    TT_train = TT_train,
    use_all_data_for_train = use_all_data_for_train,
    U_train = U_train,
    U_test = U_test,
    settings = settings,
    burn_in = burn_in,
    B = B,
    batch.size = batch.size,
    thin = thin,
    C_tune = C_tune,
    ek = ek,
    c0 = c0,
    seed = seed
  ),
  file.path(output_dir, "real_data_grid_m3_setup.rds")
)

make_metric_table <- function(results_df, metric, Season) {
  sub <- results_df[results_df$Season == Season, c("at", "MA", metric)]
  names(sub)[names(sub) == metric] <- "value"

  if (nrow(sub) == 0) {
    tab <- data.frame(at = integer())
    for (MA in MA_grid) tab[[paste0("MA", MA)]] <- numeric()
    return(tab)
  }

  tab <- reshape(
    sub,
    idvar = "at",
    timevar = "MA",
    direction = "wide"
  )
  tab <- as.data.frame(tab)
  tab <- tab[order(tab$at), ]

  ma_cols <- paste0("value.", MA_grid)
  missing_cols <- setdiff(ma_cols, names(tab))
  for (col in missing_cols) tab[[col]] <- NA_real_
  tab <- tab[, c("at", ma_cols)]
  names(tab) <- c("at", paste0("MA", MA_grid))
  tab
}

write_metric_tables <- function(results_df) {
  for (Season_s in Season_grid) {
    write.csv(
      make_metric_table(results_df, "LPML", Season_s),
      file.path(output_dir, sprintf("LPML_table_S%d.csv", Season_s)),
      row.names = FALSE
    )
    write.csv(
      make_metric_table(results_df, "WAIC", Season_s),
      file.path(output_dir, sprintf("WAIC_table_S%d.csv", Season_s)),
      row.names = FALSE
    )
  }
}

results <- vector("list", nrow(settings))

for (s in seq_len(nrow(settings))) {
  at_s <- settings$at[s]
  MA_s <- settings$MA[s]
  Season_s <- settings$Season[s]

  message("")
  message(sprintf(
    "Running setting %d/%d: at=%s, MA=%s, S=%s",
    s, nrow(settings), at_s, MA_s, Season_s
  ))

  fit <- run_mixclayton_mdim_prediction_mcmc(
    data = U_train,
    U_test = NULL,
    burn_in = burn_in,
    B = B,
    batch.size = batch.size,
    thin = thin,
    C_tune = C_tune,
    ek = ek,
    at = at_s,
    c0 = c0,
    p = NULL,
    dep_type = list(MA = MA_s, Season = c(season_lag, Season_s)),
    seed = seed,
    verbose = FALSE
  )

  result_s <- data.frame(
    setting_id = settings$setting_id[s],
    m = fit$m,
    K = fit$K,
    TT_train = TT_train,
    ek = ek,
    gk = ek,
    at = at_s,
    MA = MA_s,
    q = MA_s,
    Season = Season_s,
    WAIC = fit$WAIC$WAIC,
    WAIC_fit = fit$WAIC$fit,
    WAIC_penalty = fit$WAIC$penalty,
    LPML = fit$LPML,
    DIC = fit$DIC,
    LPS = if (!is.null(fit$LPS)) fit$LPS else NA_real_
  )

  results[[s]] <- result_s

  saveRDS(
    fit,
    file.path(output_dir, sprintf("fit_setting_%03d_at%s_MA%s_S%s.rds", s, at_s, MA_s, Season_s))
  )

  partial_df <- do.call(rbind, results[seq_len(s)])
  write.csv(
    partial_df,
    file.path(output_dir, "real_data_grid_m3_results_partial.csv"),
    row.names = FALSE
  )
  write_metric_tables(partial_df)

  message(
    "Finished setting ", s,
    ": LPML=", signif(fit$LPML, 6),
    ", WAIC=", signif(fit$WAIC$WAIC, 6),
    ", LPS=NA"
  )
}

results_df <- do.call(rbind, results)
write.csv(results_df, file.path(output_dir, "real_data_grid_m3_results.csv"), row.names = FALSE)
saveRDS(results_df, file.path(output_dir, "real_data_grid_m3_results.rds"))
write_metric_tables(results_df)

message("")
message("Done. Results written to: ", normalizePath(output_dir))
