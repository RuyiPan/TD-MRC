library(Rcpp)

args <- commandArgs(TRUE)

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file) && nzchar(script_file) && script_file != "-") {
  dirname(normalizePath(script_file))
} else if (dir.exists(file.path(getwd(), "optimized_running_time_sim"))) {
  file.path(getwd(), "optimized_running_time_sim")
} else {
  getwd()
}
source(file.path(dirname(script_dir), "simulation_helper.R"))
Rcpp::sourceCpp(file.path(script_dir, "mcmc_rcpp.cpp"))

eta_subset_fast <- function(TT, t, type) {
  partial_MA <- t:(t - type[["MA"]])
  partial_Season <- t - type[["Season"]][1] * (0:type[["Season"]][2])
  candidate <- union(partial_MA, partial_Season)
  candidate[candidate %in% seq_len(TT)]
}

inv_subset_fast <- function(TT, t, type) {
  partial_MA <- t:(t + type[["MA"]])
  partial_Season <- t + type[["Season"]][1] * (0:type[["Season"]][2])
  candidate <- union(partial_MA, partial_Season)
  candidate[candidate %in% seq_len(TT)]
}

rdirichlet_fast <- function(alpha) {
  x <- rgamma(length(alpha), shape = alpha, rate = 1)
  x / sum(x)
}

rdirichlet_matrix_fast <- function(n, alpha) {
  out <- matrix(rgamma(n * length(alpha), shape = rep(alpha, each = n), rate = 1),
                nrow = n, ncol = length(alpha))
  out / rowSums(out)
}

eta_colsum <- function(Eta, idx) {
  if (length(idx) == 1) {
    Eta[idx, ]
  } else {
    colSums(Eta[idx, , drop = FALSE])
  }
}

eta_component_sum <- function(Eta, idx, k) {
  if (length(idx) == 1) Eta[idx, k] else sum(Eta[idx, k])
}

update_eta_fast <- function(Eta, Pi_j, omega, cts, c0, p, eta_sets, inv_sets) {
  TT <- nrow(Eta)

  for (t in seq_len(TT)) {
    L <- cts[t] / 2
    current_eta <- Eta

    for (k in 1:3) {
      if (cts[t] == 1) {
        prop_etak <- sample(c(1L, 0L), size = 1)
      } else {
        prop_etak <- sample(as.integer(current_eta[t, k] - L):as.integer(current_eta[t, k] + L),
                            size = 1)
      }

      if (prop_etak <= (cts[t] - sum(current_eta[t, 1:3][-k])) && prop_etak >= 0) {
        current_etak <- current_eta[t, k]
        prop_eta <- current_eta
        prop_eta[t, k] <- prop_etak
        prop_eta[t, 4] <- cts[t] - sum(prop_eta[t, 1:3])

        inv_set <- inv_sets[[t]]
        prop_denom_k <- 0
        prop_denom_4 <- 0
        current_denom_k <- 0
        current_denom_4 <- 0

        for (l in inv_set) {
          eta_set_l <- eta_sets[[l]]
          prop_denom_k <- prop_denom_k + lgamma(c0 * p[k] + eta_component_sum(prop_eta, eta_set_l, k))
          prop_denom_4 <- prop_denom_4 + lgamma(c0 * p[4] + eta_component_sum(prop_eta, eta_set_l, 4))
          current_denom_k <- current_denom_k + lgamma(c0 * p[k] + eta_component_sum(current_eta, eta_set_l, k))
          current_denom_4 <- current_denom_4 + lgamma(c0 * p[4] + eta_component_sum(current_eta, eta_set_l, 4))
        }

        prop_dense <- lgamma(prop_etak + 1) +
          prop_denom_k +
          lgamma(prop_eta[t, 4] + 1) +
          prop_denom_4

        current_dense <- lgamma(current_etak + 1) +
          current_denom_k +
          lgamma(current_eta[t, 4] + 1) +
          current_denom_4

        ratio <- (prop_etak - current_etak) *
          (log(omega[k]) + sum(log(Pi_j[inv_set, k])) - log(omega[4]) - sum(log(Pi_j[inv_set, 4]))) +
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

run_running_time_sim_optimized <- function(job_name = "running_time_sim_optimized",
                                           job_num = 1,
                                           path = getwd(),
                                           seed = 20231213,
                                           TT = 20,
                                           nt = 300,
                                           burn_in = 60,
                                           B = 140,
                                           batch.size = 50,
                                           C_tune = 2,
                                           verbose = TRUE) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  set.seed(seed)

  U <- simulate_mixC(TT, nt)
  U_train <- U[seq_len(TT)]

  TT <- length(U_train)
  nts <- seq_len(TT)
  for (day in seq_len(TT)) {
    index <- apply(U_train[[day]], 1, function(row) any(row %in% c(0, 1)))
    U_train[[day]] <- U_train[[day]][!index, , drop = FALSE]
    nts[[day]] <- nrow(U_train[[day]])
  }

  Parameters <- rbind(
    expand.grid(c(1), c(0, 1, 3, 5, 10, 20, 30, 40), c(0, 1, 2, 3, 4, 5, 6, 7)),
    expand.grid(c(1000), c(30), c(0, 1, 2, 3, 4, 5, 6, 7))
  )

  ak <- Parameters[job_num, 1]
  ct <- Parameters[job_num, 2]
  c0 <- 1
  p <- rep(0.25, 4)
  cts <- rep(ct, TT)
  MA <- Parameters[job_num, 3]
  dep_type <- list(MA = MA, Season = c(12, 0))
  eta_sets <- lapply(seq_len(TT), eta_subset_fast, TT = TT, type = dep_type)
  inv_sets <- lapply(seq_len(TT), inv_subset_fast, TT = TT, type = dep_type)

  M <- B * batch.size
  init_z_prob <- c(0.4, 0.25, 0.1, 0.25)

  Z <- vector("list", length = TT)
  for (t in seq_len(TT)) {
    if (t == 1) {
      Z[[t]] <- t(rmultinom(nrow(U_train[[t]]), size = 1, prob = init_z_prob))
    } else {
      Z[[t]] <- matrix(NA_integer_, nrow = nrow(U_train[[t]]), ncol = 4)
    }
  }

  bk <- ak
  d <- rep(1, 4)
  a.beta <- rep(ak, 4)
  b.beta <- rep(bk, 4)
  beta <- rgamma(4, shape = 1, rate = 1)
  beta.all <- matrix(0, nrow = M, ncol = 4)
  beta.all[1, ] <- beta

  omega <- rdirichlet_fast(c0 * p)
  omega.all <- matrix(0, nrow = M, ncol = 4)
  omega.all[1, ] <- omega

  Pi <- array(dim = c(M, TT, 4))
  Pi[1, , ] <- rdirichlet_matrix_fast(TT, c0 * p)

  Theta <- array(dim = c(M, TT, 4))
  Theta[1, , ] <- rgamma(TT * 4, shape = 1, rate = 1)

  Eta <- t(rmultinom(TT, size = ct, c(0.25, 0.25, 0.25, 0.25)))
  Eta.all <- array(dim = c(M, TT, 4))
  Eta.all[1, , ] <- Eta

  acc <- matrix(0.3, nrow = TT, ncol = 4)
  acc.all <- array(dim = c(B, TT, 4))
  ada.shape <- matrix(1, nrow = TT, ncol = 4)
  kappa.all <- array(dim = c(B, TT, 4))

  for (b in seq_len(B)) {
    if (verbose) print(b)

    ada.shape[acc < 0.3] <- ada.shape[acc < 0.3] * C_tune^(1 / sqrt(b))
    ada.shape[acc > 0.4] <- ada.shape[acc > 0.4] * C_tune^(-1 / sqrt(b))
    kappa.all[b, , ] <- ada.shape
    acc.all[b, , ] <- acc
    count <- matrix(0, nrow = TT, ncol = 4)

    for (it in seq_len(batch.size)) {
      j <- it + batch.size * (b - 1)
      if (j == 1) next

      for (t in seq_len(TT)) {
        Z[[t]] <- sample_z_cpp(U_train[[t]], Theta[j - 1, t, ], Pi[j - 1, t, ])
      }

      for (t in seq_len(TT)) {
        current_theta <- Theta[j - 1, t, ]
        prop_theta <- rgamma(4, shape = ada.shape[t, ], rate = ada.shape[t, ] / current_theta)

        ll.new <- theta_component_loglik_cpp(U_train[[t]], prop_theta, Z[[t]]) +
          (d - 1) * log(prop_theta) - beta * prop_theta
        ll.old <- theta_component_loglik_cpp(U_train[[t]], current_theta, Z[[t]]) +
          (d - 1) * log(current_theta) - beta * current_theta

        g.old <- dgamma(current_theta, shape = ada.shape[t, ], rate = ada.shape[t, ] / prop_theta, log = TRUE)
        g.new <- dgamma(prop_theta, shape = ada.shape[t, ], rate = ada.shape[t, ] / current_theta, log = TRUE)
        rate <- ll.new + g.old - ll.old - g.new

        v_theta <- log(runif(4))
        accepted <- v_theta <= rate
        count[t, ] <- count[t, ] + accepted
        prop_theta[!accepted] <- current_theta[!accepted]
        Theta[j, t, ] <- prop_theta
      }

      for (t in seq_len(TT)) {
        tempPar <- c0 * p + eta_colsum(Eta, eta_sets[[t]]) + z_col_sums_cpp(Z[[t]])
        Pi[j, t, ] <- rdirichlet_fast(tempPar)
      }

      Eta <- update_eta_fast(Eta, Pi[j, , ], omega, cts, c0, p, eta_sets, inv_sets)
      Eta.all[j, , ] <- Eta

      omega <- rdirichlet_fast(c0 * p + colSums(Eta))
      omega.all[j, ] <- omega

      shapes <- a.beta + TT * d
      rates <- b.beta + colSums(Theta[j, , ])
      beta <- rgamma(4, shape = shapes, rate = rates)
      beta.all[j, ] <- beta
    }

    acc <- count / batch.size
  }

  filename_base <- paste0(job_name, "_job", job_num, "_optimized")
  res <- list(
    job_name = job_name,
    job_num = job_num,
    Parameters = Parameters[job_num, ],
    B = B,
    burn_in = burn_in,
    batch.size = batch.size,
    M = M,
    TT = TT,
    nts = nts,
    U_train = U_train,
    Theta = Theta,
    Pi = Pi,
    Eta.all = Eta.all,
    omega.all = omega.all,
    beta.all = beta.all,
    kappa.all = kappa.all,
    acc.all = acc.all
  )

  saveRDS(res, file.path(path, paste0(filename_base, ".rds")))

  invisible(res)
}

if (sys.nframe() == 0) {
  job_name <- args[1]
  job_num <- suppressWarnings(as.numeric(args[2]))
  path <- args[3]

  if (is.na(job_name)) job_name <- "running_time_sim_optimized"
  if (is.na(job_num)) job_num <- 1
  if (is.na(path)) path <- getwd()

  run_running_time_sim_optimized(
    job_name = job_name,
    job_num = job_num,
    path = path
  )
}
