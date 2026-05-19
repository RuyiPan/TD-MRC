mixture_density_matrix <- function(U, Theta_all_t, Pi_all_t) {
  if (is.null(dim(Theta_all_t))) Theta_all_t <- matrix(Theta_all_t, nrow = 1)
  if (is.null(dim(Pi_all_t))) Pi_all_t <- matrix(Pi_all_t, nrow = 1)

  L <- nrow(Theta_all_t)
  n <- nrow(U)
  out <- matrix(NA_real_, nrow = n, ncol = L)

  for (l in seq_len(L)) {
    logdens <- clayton_rot_logdens_cpp(U, Theta_all_t[l, ])
    if (ncol(logdens) != length(Pi_all_t[l, ])) {
      stop(
        "Density/component dimension mismatch: the C++ density returned ",
        ncol(logdens), " components, but Pi has ", length(Pi_all_t[l, ]),
        " components. Reinstall the package so the R and compiled C++ code ",
        "come from the same package version."
      )
    }
    out[, l] <- exp(logdens) %*% Pi_all_t[l, ]
  }

  out
}

mixture_loglik_one_draw <- function(data, Theta_t, Pi_t) {
  if (is.null(dim(Theta_t))) Theta_t <- matrix(Theta_t, nrow = 1)
  if (is.null(dim(Pi_t))) Pi_t <- matrix(Pi_t, nrow = 1)

  total <- 0
  for (t in seq_along(data)) {
    logdens <- clayton_rot_logdens_cpp(data[[t]], Theta_t[t, ])
    if (ncol(logdens) != length(Pi_t[t, ])) {
      stop(
        "Density/component dimension mismatch: the C++ density returned ",
        ncol(logdens), " components, but Pi has ", length(Pi_t[t, ]),
        " components. Reinstall the package so the R and compiled C++ code ",
        "come from the same package version."
      )
    }
    dens <- as.numeric(exp(logdens) %*% Pi_t[t, ])
    total <- total + sum(log(dens))
  }
  total
}

#' Compute WAIC for posterior draws
#'
#' @param data List of m-column matrices.
#' @param Theta Posterior array with dimensions draws x TT x `2^m`.
#' @param Pi Posterior array with dimensions draws x TT x `2^m`.
#' @export
compute_waic <- function(data, Theta, Pi) {
  fit <- 0
  penalty <- 0

  for (t in seq_along(data)) {
    dist_t <- mixture_density_matrix(data[[t]], Theta[, t, ], Pi[, t, ])
    fit <- fit + sum(log(rowMeans(dist_t)))
    penalty <- penalty + sum(apply(log(dist_t), 1, stats::var))
  }

  list(WAIC = -2 * fit + 2 * penalty, fit = fit, penalty = penalty)
}

#' Compute LPML for posterior draws
#'
#' @param data List of m-column matrices.
#' @param Theta Posterior array with dimensions draws x TT x `2^m`.
#' @param Pi Posterior array with dimensions draws x TT x `2^m`.
#' @export
compute_lpml <- function(data, Theta, Pi) {
  lpml <- 0

  for (t in seq_along(data)) {
    dist_t <- mixture_density_matrix(data[[t]], Theta[, t, ], Pi[, t, ])
    cpo_t <- 1 / rowMeans(1 / dist_t)
    lpml <- lpml + sum(log(cpo_t))
  }

  lpml
}

#' Compute DIC for posterior draws
#'
#' @param data List of m-column matrices.
#' @param Theta Posterior array with dimensions draws x TT x `2^m`.
#' @param Pi Posterior array with dimensions draws x TT x `2^m`.
#' @export
compute_dic <- function(data, Theta, Pi) {
  L <- dim(Theta)[1]
  D_all <- numeric(L)

  for (l in seq_len(L)) {
    D_all[l] <- -2 * mixture_loglik_one_draw(data, Theta[l, , ], Pi[l, , ])
  }

  PMD <- mean(D_all)
  Theta_bar <- apply(Theta, c(2, 3), mean)
  Pi_bar <- apply(Pi, c(2, 3), mean)
  D_at_mean <- -2 * mixture_loglik_one_draw(data, Theta_bar, Pi_bar)

  2 * PMD - D_at_mean
}

posterior_range <- function(burn_in, B, batch.size, thin) {
  seq(burn_in * batch.size, B * batch.size, by = thin)
}

predict_lps <- function(fit, U_test, range, remove_boundary = TRUE) {
  if (is.null(U_test)) return(NULL)

  if (remove_boundary) {
    U_test <- clean_copula_data(U_test)$data
  } else if (is.matrix(U_test)) {
    U_test <- list(U_test)
  }

  U_next <- U_test[[1]]
  predictive_dist <- NULL
  TT <- fit$TT
  K <- fit$K

  for (k in range) {
    theta_t_plus <- pmax(1e-10, rgamma(K, shape = fit$d, rate = fit$beta.all[k, ]))
    at_next <- fit$at[min(length(fit$at), TT)]
    eta_t_plus <- as.vector(rmultinom(1, at_next, fit$omega.all[k, ]))
    eta_set <- eta_subset(TT, TT + 1, fit$dep_type)
    dir_w <- fit$c0 * fit$p + eta_colsum(fit$Eta.all[k, , ], eta_set) + eta_t_plus
    weights_t_plus <- rdirichlet_fast(dir_w)
    temp_dist <- as.numeric(exp(clayton_rot_logdens_cpp(U_next, theta_t_plus)) %*% weights_t_plus)
    predictive_dist <- cbind(predictive_dist, temp_dist)
  }

  predictive_dist_mean <- rowMeans(predictive_dist)
  list(
    predictive_dist = predictive_dist,
    predictive_dist_mean = predictive_dist_mean,
    LPS = sum(log(predictive_dist_mean))
  )
}

#' Add WAIC, LPML, DIC and optional LPS to an MCMC fit
#'
#' @param fit Object returned by `run_mixclayton_mcmc()`.
#' @param thin Keep every `thin`th MCMC draw after burn-in.
#' @param U_test Optional held-out test data for LPS. Matrix or list containing one matrix.
#' @param remove_boundary_test Whether to remove boundary rows from `U_test`.
#' @export
add_prediction_scores <- function(fit,
                                  thin = 10,
                                  U_test = NULL,
                                  remove_boundary_test = TRUE) {
  range <- posterior_range(fit$burn_in, fit$B, fit$batch.size, thin)
  Theta_keep <- fit$Theta[range, , , drop = FALSE]
  Pi_keep <- fit$Pi[range, , , drop = FALSE]

  fit$thin <- thin
  fit$posterior_range <- range
  fit$WAIC <- compute_waic(fit$data, Theta_keep, Pi_keep)
  fit$LPML <- compute_lpml(fit$data, Theta_keep, Pi_keep)
  fit$DIC <- compute_dic(fit$data, Theta_keep, Pi_keep)

  pred <- predict_lps(fit, U_test, range, remove_boundary = remove_boundary_test)
  if (!is.null(pred)) {
    fit$U_test <- if (remove_boundary_test) clean_copula_data(U_test)$data else if (is.matrix(U_test)) list(U_test) else U_test
    fit$predictive_dist <- pred$predictive_dist
    fit$predictive_dist_mean <- pred$predictive_dist_mean
    fit$LPS <- pred$LPS
  }

  fit
}
