#!/usr/bin/env Rscript

# Example m = 6 simulation and TD-MRC fit.
#
# This script simulates TT = 20 time points with n_t = 300 observations and
# fits the 2^6 = 64 component rotated Clayton mixture.
#
# In RStudio:
#   setwd("/Users/ruyipan/Desktop/Project/MixCopula/Revision/TD-MRC")
#   source("MixClayton/analysis/simulation_m6/run_example_m6.R")

suppressPackageStartupMessages({
  library(TDMRC)
})

args <- commandArgs(TRUE)
output_dir <- if (length(args) >= 1) args[1] else "MixClayton/analysis/simulation_m6/results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

simulate_clayton_copula <- function(n, m, theta) {
  frailty <- rgamma(n, shape = 1 / theta, rate = 1)
  exponential_draws <- matrix(rexp(n * m), nrow = n, ncol = m)
  (1 + exponential_draws / frailty)^(-1 / theta)
}

component_bits <- function(k, m) {
  as.integer(intToBits(k - 1L)[seq_len(m)])
}

simulate_mixclayton_mdim <- function(TT,
                                     nt,
                                     m,
                                     theta = rep(3, 2^m),
                                     seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (length(nt) == 1) nt <- rep(nt, TT)
  if (length(nt) != TT) stop("`nt` must be scalar or length TT.")

  K <- 2^m
  if (length(theta) == 1) theta <- rep(theta, K)
  if (length(theta) != K) stop("`theta` must be scalar or length 2^m.")

  data <- vector("list", TT)
  base_weights <- rep(1 / K, K)

  for (t in seq_len(TT)) {
    n_t <- nt[t]
    active_center <- 1L + ((t - 1L) %% K)
    concentration <- rep(1, K)
    concentration[active_center] <- 12
    weights_t <- as.numeric(rdirichlet_fast(40 * base_weights + concentration))

    z <- sample(seq_len(K), n_t, replace = TRUE, prob = weights_t)
    U_t <- matrix(NA_real_, nrow = n_t, ncol = m)

    for (k in seq_len(K)) {
      idx <- which(z == k)
      if (length(idx) == 0) next
      U_component <- simulate_clayton_copula(length(idx), m, theta[k])
      bits <- component_bits(k, m)
      reflected <- matrix(rep(bits, each = length(idx)), nrow = length(idx), ncol = m)
      U_t[idx, ] <- ifelse(reflected == 1L, 1 - U_component, U_component)
    }

    data[[t]] <- U_t
  }

  data
}

seed <- 20260519
TT <- 20
nt <- 300
m <- 6
K <- 2^m

burn_in <- 60
B <- 140
batch.size <- 50
thin <- 10
C_tune <- 2

ek <- 1
at <- 10
c0 <- 1
MA <- 3

message("Simulating m = ", m, ", K = ", K, ", TT = ", TT, ", nt = ", nt)
U_train <- simulate_mixclayton_mdim(
  TT = TT,
  nt = nt,
  m = m,
  theta = rep(3, K),
  seed = seed
)

message(
  "Fitting TD-MRC with burn_in=", burn_in,
  ", B=", B,
  ", batch.size=", batch.size,
  ", thin=", thin
)

fit_time <- system.time({
  fit <- run_mixclayton_mdim_prediction_mcmc(
    data = U_train,
    U_test = NULL,
    burn_in = burn_in,
    B = B,
    batch.size = batch.size,
    thin = thin,
    C_tune = C_tune,
    ek = ek,
    at = at,
    c0 = c0,
    p = NULL,
    dep_type = list(MA = MA, Season = c(12, 0)),
    seed = seed,
    verbose = FALSE
  )
})

summary_df <- data.frame(
  m = fit$m,
  K = fit$K,
  TT = fit$TT,
  total_n = sum(fit$nts),
  burn_in = burn_in,
  B = B,
  batch.size = batch.size,
  thin = thin,
  ek = ek,
  at = at,
  MA = MA,
  WAIC = fit$WAIC$WAIC,
  WAIC_fit = fit$WAIC$fit,
  WAIC_penalty = fit$WAIC$penalty,
  LPML = fit$LPML,
  DIC = fit$DIC,
  elapsed_seconds = unname(fit_time[["elapsed"]])
)

saveRDS(
  list(
    U_train = U_train,
    fit = fit,
    summary = summary_df,
    seed = seed
  ),
  file.path(output_dir, "example_m6_fit.rds")
)
write.csv(summary_df, file.path(output_dir, "example_m6_summary.csv"), row.names = FALSE)

message("Done. Results written to: ", normalizePath(output_dir))
print(summary_df)
