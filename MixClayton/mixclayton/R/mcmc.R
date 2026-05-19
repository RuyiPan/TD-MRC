update_eta <- function(Eta, Pi_j, omega, ats, c0, p, eta_sets, inv_sets) {
  TT <- nrow(Eta)
  K <- ncol(Eta)
  ref_k <- K

  for (t in seq_len(TT)) {
    L <- ats[t] / 2
    current_eta <- Eta

    for (k in seq_len(K - 1)) {
      if (ats[t] == 1) {
        prop_etak <- sample(c(1L, 0L), size = 1)
      } else {
        prop_etak <- sample(as.integer(current_eta[t, k] - L):as.integer(current_eta[t, k] + L),
                            size = 1)
      }

      if (prop_etak <= (ats[t] - sum(current_eta[t, seq_len(K - 1)[-k]])) && prop_etak >= 0) {
        current_etak <- current_eta[t, k]
        prop_eta <- current_eta
        prop_eta[t, k] <- prop_etak
        prop_eta[t, ref_k] <- ats[t] - sum(prop_eta[t, seq_len(K - 1)])

        inv_set <- inv_sets[[t]]
        prop_denom_k <- 0
        prop_denom_ref <- 0
        current_denom_k <- 0
        current_denom_ref <- 0

        for (l in inv_set) {
          eta_set_l <- eta_sets[[l]]
          prop_denom_k <- prop_denom_k + lgamma(c0 * p[k] + eta_component_sum(prop_eta, eta_set_l, k))
          prop_denom_ref <- prop_denom_ref + lgamma(c0 * p[ref_k] + eta_component_sum(prop_eta, eta_set_l, ref_k))
          current_denom_k <- current_denom_k + lgamma(c0 * p[k] + eta_component_sum(current_eta, eta_set_l, k))
          current_denom_ref <- current_denom_ref + lgamma(c0 * p[ref_k] + eta_component_sum(current_eta, eta_set_l, ref_k))
        }

        prop_dense <- lgamma(prop_etak + 1) +
          prop_denom_k +
          lgamma(prop_eta[t, ref_k] + 1) +
          prop_denom_ref

        current_dense <- lgamma(current_etak + 1) +
          current_denom_k +
          lgamma(current_eta[t, ref_k] + 1) +
          current_denom_ref

        ratio <- (prop_etak - current_etak) *
          (log(omega[k]) + sum(log(Pi_j[inv_set, k])) - log(omega[ref_k]) - sum(log(Pi_j[inv_set, ref_k]))) +
          current_dense - prop_dense

        if (log(runif(1)) <= ratio) {
          current_eta <- prop_eta
        }
      }
    }

    Eta <- current_eta
  }

  Eta
}

#' Run MCMC for an m-dimensional dynamic mixture of rotated copulas
#'
#' @param data A list of matrices, one matrix per time point. The number of
#'   columns is `m`, and the sampler uses `2^m` rotated components.
#' @param burn_in Number of burn-in batches. Stored in output for downstream use.
#' @param B Number of batches.
#' @param batch.size Number of MCMC iterations per batch.
#' @param C_tune Adaptive proposal tuning multiplier.
#' @param ek Hyperparameter for beta prior, scalar or length-`2^m` vector.
#' @param at Multinomial count for eta, scalar or length `TT`.
#' @param c0 Dirichlet prior concentration multiplier.
#' @param p Baseline mixture probability vector. If `NULL`, a uniform
#'   length-`2^m` vector is used. Must sum to 1.
#' @param dep_type Dependence structure list, e.g. `list(MA = 2, Season = c(12, 0))`.
#' @param d Shape hyperparameter for theta prior, scalar or length-`2^m` vector.
#' @param init_z_prob Initial probabilities for latent allocation at time 1.
#' @param init_theta Optional `TT x 2^m` matrix of initial theta values.
#' @param init_pi Optional `TT x 2^m` matrix of initial mixture weights.
#' @param init_eta Optional `TT x 2^m` matrix of initial eta values.
#' @param init_beta Optional length-`2^m` initial beta vector.
#' @param init_omega Optional length-`2^m` initial omega vector.
#' @param remove_boundary Whether to remove rows containing 0 or 1.
#' @param compute_scores Whether to compute WAIC, LPML and DIC after MCMC.
#' @param thin Thinning interval used when computing posterior scores.
#' @param U_test Optional held-out data for LPS. If `NULL`, LPS is skipped.
#' @param remove_boundary_test Whether to remove boundary rows from `U_test`.
#' @param seed Optional random seed.
#' @param verbose Whether to print batch number.
#' @param ... Backward-compatible aliases `ak`, `ct`, and `cts`.
#' @return A list containing sampled arrays and MCMC metadata.
#' @export
run_mixclayton_mcmc <- function(data,
                                burn_in = 60,
                                B = 140,
                                batch.size = 50,
                                C_tune = 2,
                                ek = 1,
                                at = 2,
                                c0 = 1,
                                p = NULL,
                                dep_type = list(MA = 0, Season = c(12, 0)),
                                d = 1,
                                init_z_prob = NULL,
                                init_theta = NULL,
                                init_pi = NULL,
                                init_eta = NULL,
                                init_beta = NULL,
                                init_omega = NULL,
                                remove_boundary = TRUE,
                                compute_scores = FALSE,
                                thin = 10,
                                U_test = NULL,
                                remove_boundary_test = TRUE,
                                seed = NULL,
                                verbose = FALSE,
                                ...) {
  aliases <- list(...)
  if (!is.null(aliases$ak)) {
    warning("`ak` is deprecated; use `ek` instead.", call. = FALSE)
    ek <- aliases$ak
  }
  if (!is.null(aliases$ct)) {
    warning("`ct` is deprecated; use `at` instead.", call. = FALSE)
    at <- aliases$ct
  }
  if (!is.null(aliases$cts)) {
    warning("`cts` is deprecated and ignored; use `at` instead.", call. = FALSE)
  }
  unknown_aliases <- setdiff(names(aliases), c("ak", "ct", "cts"))
  if (length(unknown_aliases) > 0) {
    stop("Unknown argument(s): ", paste(unknown_aliases, collapse = ", "))
  }

  if (!is.null(seed)) set.seed(seed)

  if (remove_boundary) {
    cleaned <- clean_copula_data(data)
    U_train <- cleaned$data
    nts <- cleaned$nts
  } else {
    U_train <- if (is.matrix(data)) list(data) else data
    nts <- vapply(U_train, nrow, integer(1))
  }

  TT <- length(U_train)
  m <- ncol(U_train[[1]])
  K <- 2^m
  if (is.null(p)) p <- rep(1 / K, K)
  if (is.null(init_z_prob)) init_z_prob <- p

  validate_mcmc_inputs(U_train, burn_in, B, batch.size, C_tune, ek, at, c0, p, dep_type, m, K)

  if (length(at) == 1) at <- rep(at, TT)
  ats <- at
  if (length(ek) == 1) ek <- rep(ek, K)
  if (length(d) == 1) d <- rep(d, K)
  if (length(ek) != K) stop("`ek` must be scalar or length 2^m.")
  if (length(d) != K) stop("`d` must be scalar or length 2^m.")
  if (length(init_z_prob) != K || any(init_z_prob <= 0)) {
    stop("`init_z_prob` must be a positive length-2^m vector.")
  }
  init_z_prob <- init_z_prob / sum(init_z_prob)

  eta_sets <- lapply(seq_len(TT), eta_subset, TT = TT, type = dep_type)
  inv_sets <- lapply(seq_len(TT), inv_subset, TT = TT, type = dep_type)

  M <- B * batch.size
  Z <- vector("list", length = TT)
  for (t in seq_len(TT)) {
    if (t == 1) {
      Z[[t]] <- t(rmultinom(nrow(U_train[[t]]), size = 1, prob = init_z_prob))
    } else {
      Z[[t]] <- matrix(NA_integer_, nrow = nrow(U_train[[t]]), ncol = K)
    }
  }

  a.beta <- ek
  b.beta <- ek

  beta <- if (is.null(init_beta)) rgamma(K, shape = 1, rate = 1) else init_beta
  beta.all <- matrix(0, nrow = M, ncol = K)
  beta.all[1, ] <- beta

  omega <- if (is.null(init_omega)) rdirichlet_fast(c0 * p) else init_omega
  omega.all <- matrix(0, nrow = M, ncol = K)
  omega.all[1, ] <- omega

  Pi <- array(dim = c(M, TT, K))
  Pi[1, , ] <- if (is.null(init_pi)) rdirichlet_matrix_fast(TT, c0 * p) else init_pi

  Theta <- array(dim = c(M, TT, K))
  Theta[1, , ] <- if (is.null(init_theta)) rgamma(TT * K, shape = 1, rate = 1) else init_theta

  Eta <- if (is.null(init_eta)) {
    t(vapply(ats, function(at_t) as.integer(rmultinom(1, size = at_t, prob = rep(1 / K, K))), integer(K)))
  } else {
    init_eta
  }
  Eta.all <- array(dim = c(M, TT, K))
  Eta.all[1, , ] <- Eta

  acc <- matrix(0.3, nrow = TT, ncol = K)
  acc.all <- array(dim = c(B, TT, K))
  ada.shape <- matrix(1, nrow = TT, ncol = K)
  kappa.all <- array(dim = c(B, TT, K))

  for (b in seq_len(B)) {
    if (verbose) print(b)

    ada.shape[acc < 0.3] <- ada.shape[acc < 0.3] * C_tune^(1 / sqrt(b))
    ada.shape[acc > 0.4] <- ada.shape[acc > 0.4] * C_tune^(-1 / sqrt(b))
    kappa.all[b, , ] <- ada.shape
    acc.all[b, , ] <- acc
    count <- matrix(0, nrow = TT, ncol = K)

    for (it in seq_len(batch.size)) {
      j <- it + batch.size * (b - 1)
      if (j == 1) next

      for (t in seq_len(TT)) {
        Z[[t]] <- sample_z_cpp(U_train[[t]], Theta[j - 1, t, ], Pi[j - 1, t, ])
      }

      for (t in seq_len(TT)) {
        current_theta <- Theta[j - 1, t, ]
        prop_theta <- rgamma(K, shape = ada.shape[t, ], rate = ada.shape[t, ] / current_theta)

        ll.new <- theta_component_loglik_cpp(U_train[[t]], prop_theta, Z[[t]]) +
          (d - 1) * log(prop_theta) - beta * prop_theta
        ll.old <- theta_component_loglik_cpp(U_train[[t]], current_theta, Z[[t]]) +
          (d - 1) * log(current_theta) - beta * current_theta

        g.old <- dgamma(current_theta, shape = ada.shape[t, ], rate = ada.shape[t, ] / prop_theta, log = TRUE)
        g.new <- dgamma(prop_theta, shape = ada.shape[t, ], rate = ada.shape[t, ] / current_theta, log = TRUE)
        rate <- ll.new + g.old - ll.old - g.new

        v_theta <- log(runif(K))
        accepted <- v_theta <= rate
        count[t, ] <- count[t, ] + accepted
        prop_theta[!accepted] <- current_theta[!accepted]
        Theta[j, t, ] <- prop_theta
      }

      for (t in seq_len(TT)) {
        tempPar <- c0 * p + eta_colsum(Eta, eta_sets[[t]]) + z_col_sums_cpp(Z[[t]])
        Pi[j, t, ] <- rdirichlet_fast(tempPar)
      }

      Eta <- update_eta(Eta, Pi[j, , ], omega, ats, c0, p, eta_sets, inv_sets)
      Eta.all[j, , ] <- Eta

      omega <- rdirichlet_fast(c0 * p + colSums(Eta))
      omega.all[j, ] <- omega

      shapes <- a.beta + TT * d
      rates <- b.beta + colSums(Theta[j, , ])
      beta <- rgamma(K, shape = shapes, rate = rates)
      beta.all[j, ] <- beta
    }

    acc <- count / batch.size
  }

  fit <- list(
    call = match.call(),
    data = U_train,
    nts = nts,
    m = m,
    K = K,
    TT = TT,
    B = B,
    burn_in = burn_in,
    batch.size = batch.size,
    M = M,
    C_tune = C_tune,
    ek = ek,
    at = at,
    c0 = c0,
    p = p,
    dep_type = dep_type,
    d = d,
    Theta = Theta,
    Pi = Pi,
    Eta.all = Eta.all,
    omega.all = omega.all,
    beta.all = beta.all,
    kappa.all = kappa.all,
    acc.all = acc.all
  )

  class(fit) <- "mixclayton_mcmc"

  if (compute_scores || !is.null(U_test)) {
    fit <- add_prediction_scores(
      fit = fit,
      thin = thin,
      U_test = U_test,
      remove_boundary_test = remove_boundary_test
    )
  }

  fit
}

#' Run MCMC for m-dimensional rotated copula mixtures
#'
#' This is a descriptive alias for `run_mixclayton_mcmc()`. It supports
#' `m = 2, 3, 4, 5`, inferred from the number of data columns.
#'
#' @inheritParams run_mixclayton_mcmc
#' @export
run_mixclayton_mdim_mcmc <- function(...) {
  run_mixclayton_mcmc(...)
}

#' Run MCMC and compute prediction/model-comparison scores
#'
#' Convenience wrapper around `run_mixclayton_mcmc()` with
#' `compute_scores = TRUE`.
#'
#' @inheritParams run_mixclayton_mcmc
#' @export
run_mixclayton_prediction_mcmc <- function(data,
                                           burn_in = 60,
                                           B = 140,
                                           batch.size = 50,
                                           C_tune = 2,
                                           ek = 1,
                                           at = 2,
                                           c0 = 1,
                                           p = NULL,
                                           dep_type = list(MA = 0, Season = c(12, 0)),
                                           d = 1,
                                           thin = 10,
                                           U_test = NULL,
                                           seed = NULL,
                                           verbose = FALSE,
                                           ...) {
  run_mixclayton_mcmc(
    data = data,
    burn_in = burn_in,
    B = B,
    batch.size = batch.size,
    C_tune = C_tune,
    ek = ek,
    at = at,
    c0 = c0,
    p = p,
    dep_type = dep_type,
    d = d,
    compute_scores = TRUE,
    thin = thin,
    U_test = U_test,
    seed = seed,
    verbose = verbose,
    ...
  )
}

#' Run m-dimensional MCMC and compute prediction/model-comparison scores
#'
#' This is a descriptive alias for `run_mixclayton_prediction_mcmc()`. It
#' supports `m = 2, 3, 4, 5`, inferred from the number of data columns.
#'
#' @inheritParams run_mixclayton_prediction_mcmc
#' @export
run_mixclayton_mdim_prediction_mcmc <- function(...) {
  run_mixclayton_prediction_mcmc(...)
}
