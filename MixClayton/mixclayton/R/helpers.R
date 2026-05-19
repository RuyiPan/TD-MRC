#' Remove boundary observations from copula data
#'
#' @param data A list of copula-data matrices, or a single matrix.
#' @return A list with cleaned `data` and sample sizes `nts`.
#' @export
clean_copula_data <- function(data) {
  if (is.matrix(data)) data <- list(data)
  stopifnot(is.list(data))

  cleaned <- vector("list", length(data))
  nts <- integer(length(data))

  for (t in seq_along(data)) {
    U_t <- as.matrix(data[[t]])
    if (ncol(U_t) < 2) stop("Each data element must have at least two columns.")
    index <- apply(U_t, 1, function(row) any(row %in% c(0, 1)))
    cleaned[[t]] <- U_t[!index, , drop = FALSE]
    nts[t] <- nrow(cleaned[[t]])
  }

  list(data = cleaned, nts = nts)
}

#' Find eta indices entering a time-varying weight
#'
#' @param TT Number of time points.
#' @param t Current time index.
#' @param type Dependence structure list, e.g. `list(MA = 2, Season = c(12, 0))`.
#' @export
eta_subset <- function(TT, t, type) {
  partial_MA <- t:(t - type[["MA"]])
  partial_Season <- t - type[["Season"]][1] * (0:type[["Season"]][2])
  candidate <- union(partial_MA, partial_Season)
  candidate[candidate %in% seq_len(TT)]
}

#' Find inverse eta indices affected by a proposed eta update
#'
#' @param TT Number of time points.
#' @param t Current time index.
#' @param type Dependence structure list, e.g. `list(MA = 2, Season = c(12, 0))`.
#' @export
inv_subset <- function(TT, t, type) {
  partial_MA <- t:(t + type[["MA"]])
  partial_Season <- t + type[["Season"]][1] * (0:type[["Season"]][2])
  candidate <- union(partial_MA, partial_Season)
  candidate[candidate %in% seq_len(TT)]
}

#' Draw one Dirichlet vector
#'
#' @param alpha Positive Dirichlet parameters.
#' @export
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

#' Generate the reflection index table for m-dimensional components
#'
#' The row order matches the m-dimensional scripts in `MCMC_MRC`: for `m = 3`
#' the rows are `000, 001, 010, 011, 100, 101, 110, 111`.
#'
#' @param m Copula dimension.
#' @return A `2^m x m` integer matrix. A 1 means use `1 - u_j` for that
#'   component and a 0 means use `u_j`.
#' @export
generate_component_index <- function(m) {
  if (!(m %in% 2:5)) stop("This package currently supports m = 2, 3, 4, or 5.")
  idx <- as.matrix(expand.grid(rep(list(c(0L, 1L)), m)))
  idx <- idx[do.call(order, as.data.frame(idx)), , drop = FALSE]
  storage.mode(idx) <- "integer"
  rownames(idx) <- NULL
  colnames(idx) <- paste0("u", seq_len(m))
  idx
}

#' Precompute all reflected data arrays
#'
#' This helper mirrors `precompute_U_transformed()` in the m-dimensional
#' reference scripts. It is not required by the Rcpp sampler, but is useful for
#' checking component ordering or external density calculations.
#'
#' @param data A matrix or list of matrices.
#' @return A list of arrays with dimensions `n_t x 2^m x m`.
#' @export
precompute_reflected_data <- function(data) {
  if (is.matrix(data)) data <- list(data)
  if (!is.list(data)) stop("`data` must be a matrix or list of matrices.")

  m <- ncol(data[[1]])
  comp_dict <- generate_component_index(m)
  out <- vector("list", length(data))

  for (t in seq_along(data)) {
    U_t <- as.matrix(data[[t]])
    if (ncol(U_t) != m) stop("All data matrices must have the same number of columns.")
    n_t <- nrow(U_t)
    U_transformed <- array(NA_real_, dim = c(n_t, nrow(comp_dict), m))

    for (i in seq_len(n_t)) {
      u <- U_t[i, ]
      U_transformed[i, , ] <- t(apply(comp_dict, 1, function(b) {
        b * (1 - u) + (1 - b) * u
      }))
    }

    out[[t]] <- U_transformed
  }

  out
}

validate_mcmc_inputs <- function(data, burn_in, B, batch.size, C_tune,
                                 ek, at, c0, p, dep_type, m, K) {
  if (!is.list(data)) stop("`data` must be a list of matrices after cleaning.")
  data_dims <- vapply(data, ncol, integer(1))
  if (any(data_dims != m)) stop("All data matrices must have the same number of columns.")
  if (!(m %in% 2:5)) stop("This package currently supports m = 2, 3, 4, or 5.")
  if (length(p) != K || any(p <= 0) || abs(sum(p) - 1) > 1e-8) {
    stop("`p` must be a positive length-2^m vector that sums to 1.")
  }
  if (!is.list(dep_type) || is.null(dep_type[["MA"]]) || is.null(dep_type[["Season"]])) {
    stop("`dep_type` must be a list like list(MA = q, Season = c(s, p)).")
  }
  if (length(dep_type[["Season"]]) != 2) stop("`dep_type$Season` must have length 2.")
  if (any(c(burn_in, B, batch.size) <= 0)) stop("`burn_in`, `B`, and `batch.size` must be positive.")
  if (C_tune <= 0) stop("`C_tune` must be positive.")
  if (any(ek <= 0) || c0 <= 0) stop("`ek` and `c0` must be positive.")
  if (length(at) != 1 && length(at) != length(data)) stop("`at` must be scalar or length TT.")
  invisible(TRUE)
}
