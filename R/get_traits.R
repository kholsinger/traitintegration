#' Report the traits found in the model
#'
#' This is a utility function to make it easy to find the traits used in a
#' `brms` analysis and to report them in the order they will appear when
#' the posterior correlation matrix is plotted
#'
#' @param output_model An object of class `brmsfit`
#'
#' @return A vector of parameter names
#'
#' @export
report_trait_names <- function(output_model) {
  Omega <- get_Omega(output_model)
  model_df <- as.data.frame(output_model)
  tmp <- model_df %>%
    dplyr::select(starts_with("b_"))
  traits <- gsub("b_(.*)_Intercept", "\\1", colnames(tmp))
  return(traits)
}
