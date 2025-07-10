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
#' of eigenvalues, and ent_R statistic for each of the correlation matrices
#' in the brmsfit object
#'
#' @examples
#' assemble_statistics(mod_wild)
#'
#' @export
assemble_statistics <- function(fit) {
  Omega <- get_Omega(fit)
  eigen_stats <- get_eigen(Omega)
  ent_R <- get_ent_R(Omega)
  df <- cbind(eigen_stats, ent_R)
  return(df)
}
