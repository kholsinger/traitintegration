#' Return eigenvalue statistics from a set of posterior correlation
#' matrices
#'
#' This function takes an n_iter x n_traits x n_traits set of correlation
#' matrices and returns a data frame. The first n_trait columns of the data
#' frame are the eigenvalues (prefixed with "EV_"). The remaining two are
#' the lead eigenvalue (Lead_EV), and the variance of eigenvalues (Var_EV).
#'
#' @param Omega An array of correlation matrices, typically from get_Omega().
#'
#' @return A data frame with the eigenvalues, lead eigenvalue, and variance
#' of eigenvalues for each of the correlation matrices
#'
#' @export
get_eigen <- function(Omega){
  eigen <- numeric(0)
  n_iter <- dim(Omega)[1]
  n_traits <- dim(Omega)[2]
  for(i in 1:n_iter){
    A <- eigen(Omega[i,,])
    eigen <- rbind(eigen, A$values)
  }

  eigen_df <- as.data.frame(eigen)
  prefix <- "EV"
  suffix <- seq(1:n_traits)
  colnames(eigen_df) <- paste(prefix, suffix, sep = "_")

  eigen_df <- eigen_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Lead_EV = max(dplyr::c_across(cols = 1:all_of(n_traits))),
                  Var_EV = var(dplyr::c_across(cols = 1:all_of(n_traits))))
  return(eigen_df)
}
