#' Assemble eigenvalue and entR statistics from posterior samples
#'
#' @description
#' This function encapsulates all of the computations associated with
#' assessing phenotypic integration using a brmsfit object from a
#' multiresponse model. For example, the wild sunflower results distributed
#' with this package (in `mod_wild`) were produced from individual-level data
#' with the following set of commands:
#'
#' \preformatted{
#' wild_CN <- bf(Leaf.C.N.ratio ~ (1|population_id)
#' wild_SLA <- bf(SLA ~ (1|population_id))
#' wild_branch <- bf(Primary.branches ~ (1|population_id))
#' wild_stem <- bf(Stem.diameter.at.flowering ~ (1|population_id))
#' wild_ligule <- bf(Ligule.length ~ (1|population_id))
#' wild_phyllary <- bf(Phyllaries.length ~ (1|population_id))
#'
#' mod_wild <- brm(wild_CN + wild_SLA + wild_branch + wild_stem +
#'                 wild_ligule + wild_phyllary + set_rescor(rescor = TRUE),
#'                 data = wild,
#'                 family = "gaussian")
#'}
#'
#' Once you've loaded this library, `mod_wild` is available for use. See the
#' example below. In addition, if you are familiar with `brms`, you can use
#' any of the usual `brms` methods to learn more about `mod_wild`.
#'
#' @param fit A brmsfit object from a multiresponse model
#'
#' @return A data frame with the eigenvalues, lead eigenvalue, variance
#' of eigenvalues, and r_eff statistic for each of the correlation matrices
#' in the brmsfit object
#'
#' @examples
#' assemble_statistics(mod_wild)
#'
#' @export
assemble_statistics <- function(fit) {
  Omega <- get_Omega(fit)
  eigen_stats <- get_eigen(Omega)
  r_eff <- get_r_eff(Omega)
  df <- cbind(eigen_stats, r_eff)
  return(df)
}

#' Assemble eigenvalue and entR statistics from bootstrap samples
#'
#' @description
#' This function encapsulates all of the computations associated with
#' assessing phenotypic integration using a bootstrap sample of correlation
#' matrices. See `vignette("Bootstrapping", package = "traitintegration")`
#' for a complete example.
#'
#' @param sample A bootstrap sample of correlation matrices. Each row
#' is the bootstrap sample for one element of the correlation matrix.
#'
#' @return A data frame with the eigenvalues, lead eigenvalue, variance
#' of eigenvalues, and r_eff statistic for each of the correlation matrices
#' in the bootstrap sample
#'
#' @export
assemble_statistics_boot <- function(Omega) {
  eigen_stats <- get_eigen(Omega)
  r_eff <- get_r_eff(Omega)
  df <- cbind(eigen_stats, r_eff)
  return(df)
}
