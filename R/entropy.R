#' Returns the entropy of a multivariate normal distribution with a specified
#' covariance matrix
#'
#' @param x A k x k covariance matrix
#'
#' @return The entropy of the corresponding multivariate normal distribution
#' @export
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
#' @param tol Stop the search when the entropy difference is less than tol
#'
#' @return A vector of r equivalents, i.e., the pairwise corrleation in a
#' constant correlation matrix producing entropies equivalent to those in
#' the vector of entropies supplied
#'
#' @export
get_ent_R <- function(Omega, tol = 1.0e-14) {
  n_ent <- dim(Omega)[1]
  k <- dim(Omega)[2]
  ent_obs <- get_entropy(Omega)
  ent_R <- numeric(n_ent)
  diff <- 1.0
  min_r <- 1.0e-16
  max_r <- 1 - min_r
  f_min_r <- entropy_omega(min_r, k)
  f_max_r <- entropy_omega(max_r, k)
  for (i in 1:n_ent) {
    min <- min_r
    max <- max_r
    while (diff > tol)  {
      mid <- (max + min)/2.0
      f_mid_r <- entropy_omega(mid, k)
      if (f_mid_r < ent_obs[i]) {
        max <- mid
        f_max_r <- f_mid_r
      } else {
        min <- mid
        f_min_r <- f_mid_r
      }
      diff <- abs(f_mid_r - ent_obs[i])
    }
    ent_R[i] <- mid
    diff <- 1.0
  }
  return(ent_R)
}
