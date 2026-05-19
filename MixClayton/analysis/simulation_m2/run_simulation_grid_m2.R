#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(TDMRC)
})

args <- commandArgs(TRUE)
output_dir <- if (length(args) >= 1) args[1] else "MixClayton/analysis/simulation_m2/results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# MCMC settings
burn_in <- 60
B <- 140
batch.size <- 50
thin <- 10
C_tune <- 2

# Data settings
data_seed <- 20231213
TT <- 20
nt <- 300

# Table settings: e_k = g_k = ek, a_t = at, q = MA.
# The table has only q = 0 for at = 0 because MA is irrelevant when at = 0.
settings <- rbind(
  data.frame(ek = 1,    at = 0,  MA = 0),
  expand.grid(ek = 1,    at = c(1, 10, 20, 30, 40), MA = 0:7),
  expand.grid(ek = 1000, at = 30,                    MA = 0:7)
)
settings <- settings[order(settings$ek, settings$at, settings$MA), ]
settings$setting_id <- seq_len(nrow(settings))

set.seed(data_seed)
U_train <- simulate_mixclayton(TT = TT, nt = nt)

saveRDS(
  list(
    U_train = U_train,
    data_seed = data_seed,
    settings = settings,
    burn_in = burn_in,
    B = B,
    batch.size = batch.size,
    thin = thin,
    C_tune = C_tune
  ),
  file.path(output_dir, "simulation_grid_setup.rds")
)

results <- vector("list", nrow(settings))

for (s in seq_len(nrow(settings))) {
  ek_s <- settings$ek[s]
  at_s <- settings$at[s]
  MA_s <- settings$MA[s]

  message(sprintf(
    "Running setting %d/%d: ek=gk=%s, at=%s, MA=q=%s",
    s, nrow(settings), ek_s, at_s, MA_s
  ))

  fit <- run_mixclayton_prediction_mcmc(
    data = U_train,
    U_test = NULL,
    burn_in = burn_in,
    B = B,
    batch.size = batch.size,
    thin = thin,
    C_tune = C_tune,
    ek = ek_s,
    at = at_s,
    c0 = 1,
    p = rep(0.25, 4),
    dep_type = list(MA = MA_s, Season = c(12, 0)),
    seed = data_seed,
    verbose = FALSE
  )

  result_s <- data.frame(
    setting_id = settings$setting_id[s],
    ek = ek_s,
    gk = ek_s,
    at = at_s,
    q = MA_s,
    WAIC = fit$WAIC$WAIC,
    WAIC_fit = fit$WAIC$fit,
    WAIC_penalty = fit$WAIC$penalty,
    LPML = fit$LPML,
    DIC = fit$DIC
  )

  results[[s]] <- result_s

  saveRDS(
    fit,
    file.path(output_dir, sprintf("fit_setting_%03d_ek%s_at%s_q%s.rds", s, ek_s, at_s, MA_s))
  )

  write.csv(
    do.call(rbind, results[seq_len(s)]),
    file.path(output_dir, "simulation_grid_results_partial.csv"),
    row.names = FALSE
  )
}

results_df <- do.call(rbind, results)
write.csv(results_df, file.path(output_dir, "simulation_grid_results.csv"), row.names = FALSE)
saveRDS(results_df, file.path(output_dir, "simulation_grid_results.rds"))

lpml_table <- reshape(
  results_df[, c("ek", "gk", "at", "q", "LPML")],
  idvar = c("ek", "gk", "at"),
  timevar = "q",
  direction = "wide"
)
waic_table <- reshape(
  results_df[, c("ek", "gk", "at", "q", "WAIC")],
  idvar = c("ek", "gk", "at"),
  timevar = "q",
  direction = "wide"
)

write.csv(lpml_table, file.path(output_dir, "LPML_table.csv"), row.names = FALSE)
write.csv(waic_table, file.path(output_dir, "WAIC_table.csv"), row.names = FALSE)

message("Done. Results written to: ", normalizePath(output_dir))
