#' Returns the entropy of a multivariate normal distribution with a specified
#' covariance matrix
#'
#' @param x A k x k covariance matrix
#'
#' @return The entropy of the corresponding multivariate normal distribution
entropy <- function(x) {
  k <- ncol(x)
  ent <- k/2 + (k/2)*log(2*pi) + (1/2)*log(det(x))
  return(ent)
}

#' Returns the entropies for a set of covariance matrices
#'
#' @param Omega An n_iter x n_traits x n_traits array of covariance matrices
#'
#' @return A numeric vector containing the entropy associated with each
#' covariance matrix
#'
#' @export
get_entropy <- function(Omega) {
  n_iter <- dim(Omega)[1]
  ent <- numeric(n_iter)
  for (i in 1:n_iter) {
    ent[i] <- entropy(Omega[i, , ])
  }
  return(ent)
}

#' Calculate the entropy for a k x k correlation matrix with a constant
#' (positive) correlation
#'
#' @param r The correlation coefficient between all pairs of variables
#' @param k The dimension of the correlation matrix
#'
#' @return The entropy of the corresponding multivariate normal distribution
#'
#' @export
entropy_omega <- function(r, k) {
  ent_omega <- diag(k)
  for (i in 1:(k-1)) {
    for (j in (i+1):k) {
      ent_omega[i, j] <- r
      ent_omega[j, i] <- r
    }
  }
  return(entropy(ent_omega))
}

#' Calculate the correlation coefficient a multivariate normal with all
#' pairwise correlations equal that has an entropy equivalent to the
#' observed entropy.
#'
#' @param Omega An n_iter x n_traits x n_traits array of correlation matrices
#' @param min_ent The minimum entropy at which to begin the binary search
#' @param max_ent The maximum entropy at which to begin the binary search
#' @param tol Stop the search when the entropy difference is less than tol
#'
#' @return A vector of r equivalents, i.e., the pairwise corrleation in a
#' constant correlation matrix producing entropies equivalent to those in
#' the vector of entropies supplied
#'
#' @export
get_ent_R <- function(Omega, min_ent = 0.01, max_ent = 0.99, tol = 1.0e-14) {
  n_ent <- dim(Omega)[1]
  k <- dim(Omega)[2]
  ent_obs <- get_entropy(Omega)
  ent_R <- numeric(n_ent)
  diff <- 1.0
  for (i in 1:n_ent) {
    min <- min_ent
    max <- max_ent
    mid <- (max + min)/2.0
    while (diff > tol)  {
      min_val <- entropy_omega(min, k)
      mid_val <- entropy_omega(mid, k)
      max_val <- entropy_omega(max, k)
      stopifnot((min_val > ent_obs[i]) & (max_val < ent_obs[i]))
      if (mid_val < ent_obs[i]) {
        max <- mid
      } else {
        min <- mid
      }
      mid <- (max + min)/2.0
      diff <- abs(mid_val - ent_obs[i])
    }
    ent_R[i] <- mid
    diff <- 1.0
  }
  return(ent_R)
}
