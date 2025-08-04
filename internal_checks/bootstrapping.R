rm(list = ls())

## retrieve the phenotype matrix
##
phenos <- readr::read_csv("http://www.helianthome.org/rest/population/phenotype_matrix/wild",
                          show_col_types = FALSE)
## retrieve individual id and population id
##
index <- readr::read_csv("http://www.helianthome.org/rest/individual/list/",
                          show_col_types = FALSE)
## merge population id into phenotypic data. drop individuals lacking
## complete data
## relocate puts population_id in the first column, which simplifies
## dropping it for the correlation calculations
##
wild <- merge(phenos, index, by = "individual_id") |>
  dplyr::select(`Leaf C N ratio`, SLA, `Primary branches`,
                `Stem diameter at flowering`, `Ligule length`,
                `Phyllaries length`, population_id) |>
  dplyr::rename(Leaf.C.N.ratio = `Leaf C N ratio`,
                Primary.branches = `Primary branches`,
                Stem.diameter.at.flowering = `Stem diameter at flowering`,
                Ligule.length = `Ligule length`,
                Phyllaries.length = `Phyllaries length`) |>
  dplyr::filter(!is.na(Leaf.C.N.ratio) &
                !is.na(SLA) &
                !is.na(Primary.branches) &
                !is.na(Stem.diameter.at.flowering) &
                !is.na(Ligule.length) &
                !is.na(Phyllaries.length)) |>
  dplyr::relocate("population_id")

scale_within_population <- function(dat, pop) {
  tmp <- subset(dat, population_id == pop)
  tmp <- as.data.frame(scale(tmp[, -1]))
  col_names <- c(colnames(tmp), "population_id")
  tmp <- tibble::add_column(tmp, pop)
  colnames(tmp) <- col_names
  tmp <- tmp %>%
    dplyr::relocate(population_id)
  return(tmp)
}

## scaling may fail either because there is only one row in the data frame
## (hence no standard deviation for the denominator), producing an NaN,
## or because a particular trait has no variation within the data frame,
## producing an NA.
##
scaling_valid <- function(dat) {
  tmp <- as.matrix(dat[, -1])
  result <- sum(is.nan(tmp))
  result <- result & sum(is.na(tmp))
  return(result == 0)
}

scale_within_populations <- function(dat) {
  pops <- unique(dat$population_id)
  tmp <- scale_within_population(dat, pops[1])
  for (i in 2:length(pops)) {
    scaled <- scale_within_population(dat, pops[i])
    ## don't include populations where scaling fails
    ##
    if (scaling_valid(scaled)) {
      tmp <- tibble::add_row(tmp, scaled)
    }
  }
  return(tmp)
}

correlation_within_populations <- function(dat) {
  tmp <- scale_within_populations(dat)
  R <- cor(tmp[, -1], use = "complete.obs")
  return(R)
}

bootstrap <- function(n_sample, dat) {
  n_traits <- ncol(dat) - 1
  Omega <- array(dim = c(n_sample, n_traits, n_traits),
                 dimnames = list(paste("Sample", 1:n_sample, sep = ""),
                                 colnames(dat)[-1],
                                 colnames(dat)[-1]))
  for (i in 1:n_sample) {
    cat("Sample ", i, " of ", n_sample, "\n", sep = "")
    rows <- sample(1:nrow(dat), nrow(dat), replace = TRUE)
    R <- correlation_within_populations(dat[rows, ])
    Omega[i, , ] <- R
  }
  return(Omega)
}

n_sample <- 1000
Omega <- bootstrap(n_sample,  wild)

## check correlation matrices to make sure we can get eigenvalues
##
for (i in 1:n_sample) {
  res <- try(eigen(Omega[i, , ]), silent = TRUE)
  if (inherits(res, "try-error")) {
    cat("Omega[", i, ", , ] not valid\n")
    print(Omega[i, , ])
    stop()
  }
}

results <- assemble_statistics_boot(Omega)
summarize_statistics(results)

p <- bayesplot::mcmc_areas(results[, c("Lead_EV", "Var_EV", "ent_R")])
print(p)

post_corr <- plot_posterior_correlation_boot(Omega)
print(round(post_corr$R, 3))



