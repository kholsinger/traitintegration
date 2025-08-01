#' Extract array of correlation matrices from output of brm()
#'
#' @description
#' This function takes an object of class brmsfit produced by fitting a
#' multivariate response model, extracts the correlations, checks that all of
#' the correlation matrices are positive definite (printing a warning for any
#' that aren't), and returns them.
#'
#' @param output_model An object of class brmsfit
#'
#' @return An n_iter x trait x trait array
#'
#' @importFrom magrittr %>%
#'
#' @export
get_Omega <- function(output_model) {
  model_df <- as.data.frame(output_model)
  tmp <- model_df %>%
    dplyr::select(starts_with("b_"))
  traits <- gsub("b_(.*)_Intercept", "\\1", colnames(tmp))

  n_trait <- length(traits)
  n_iter <- dim(model_df)[1]

  Omega <- array(dim = c(n_iter, n_trait, n_trait))
  for (i in 1:(n_trait-1)) {
    Omega[ , i, i] <- 1.0
    for (j in (i+1):n_trait) {
      colname <- paste("rescor__", traits[i],"__",traits[j], sep = "")
      tmp_df <- model_df %>% dplyr::select(all_of(colname))
      Omega[ , i, j] <- tmp_df[[colname]]
      Omega[ , j, i] <- tmp_df[[colname]]
    }
  }
  Omega[ , n_trait, n_trait] <- 1.0

  for (i in 1:n_iter) {
    if (!matrixcalc::is.positive.definite(Omega[i, , ])) {
      print("Warning: not positive definite")
      print(Omega[i, , ])
    }
  }
  return(Omega)
}

#' Construct an array of correlation matrices from a bootstrap sample of
#' correlation matrices.
#'
#' @description
#' This function takes an matrix consisting of bootstrap samples of a
#' correlation matrix with each row consisting of a bootstrap sample with
#' for one element of the correlation matrix. Note: The number of rows
#' in the matrix must be a perfect square, and the square root of the number
#' of rows must be equal to the number of traits.
#'
#' @param sample The bootstrap sample
#'
#' @return An n_iter x trait x trait array
#'
#' @importFrom magrittr %>%
#'
#' @export
get_Omega_boot <- function(sample) {
  thetastar <- sample$thetastar
  n_sample <- ncol(thetastar)
  n_traits <- sqrt(nrow(thetastar))
  Omega <- array(dim = c(n_sample, n_traits, n_traits))
  for (i in 1:n_sample) {
    ct <- 0
    for (j in 1:n_traits) {
      for (k in 1:n_traits) {
        ct <- ct + 1
        Omega[i, j, k] <- thetastar[ct, i]
      }
    }
  }
  return(Omega)
}
