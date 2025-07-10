#' Extract array of correlation matrices from output of brm()
#'
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
