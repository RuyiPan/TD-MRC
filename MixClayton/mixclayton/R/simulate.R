vf_clayton <- function(ty, th, p, q) {
  if (ty == 1) {
    (p^(-th) * (q^(-th / (th + 1)) - 1) + 1)^(-1 / th)
  } else if (ty == 4) {
    1 - (p^(-th) * ((1 - q)^(-th / (th + 1)) - 1) + 1)^(-1 / th)
  } else if (ty == 3) {
    1 - ((1 - p)^(-th) * ((1 - q)^(-th / (th + 1)) - 1) + 1)^(-1 / th)
  } else if (ty == 2) {
    ((1 - p)^(-th) * (q^(-th / (th + 1)) - 1) + 1)^(-1 / th)
  } else {
    stop("Unknown rotation type.")
  }
}

#' Simulate from the four-way rotated Clayton mixture used in the paper examples
#'
#' @param TT Number of time points.
#' @param nt Number of observations per time point, or vector of length `TT`.
#' @param theta Length-4 vector of component Clayton parameters.
#' @param weights0 Initial length-4 mixture weights.
#' @param verbose Print evolving weights.
#' @export
simulate_mixclayton <- function(TT,
                                nt,
                                theta = c(5, 3, 4, 3),
                                weights0 = c(0.4, 0.25, 0.1, 0.25),
                                verbose = FALSE) {
  if (length(nt) == 1) nt <- rep(nt, TT)
  if (length(nt) != TT) stop("`nt` must be scalar or length TT.")

  U <- vector("list", length = TT)
  weights <- weights0
  rotations <- 1:4

  for (t in seq_len(TT)) {
    if (t != 1) {
      weights[1] <- weights[1] * 0.95
      weights[2] <- weights[2] * 1.05
      weights[3] <- weights[3]
      weights[4] <- 1 - sum(weights[1:3])
    }
    if (verbose) print(weights)

    n_t <- nt[t]
    p <- matrix(runif(n_t * 4), nrow = n_t, ncol = 4)
    q <- matrix(runif(n_t * 4), nrow = n_t, ncol = 4)
    v <- matrix(NA_real_, nrow = n_t, ncol = 4)

    for (l in rotations) {
      v[, l] <- vf_clayton(rotations[l], theta[l], p[, l], q[, l])
    }

    z <- sample(rotations, n_t, replace = TRUE, prob = weights)
    U[[t]] <- cbind(p[cbind(seq_len(n_t), z)], v[cbind(seq_len(n_t), z)])
  }

  U
}
