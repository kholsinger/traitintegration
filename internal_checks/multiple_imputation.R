## retrieve the phenotype matrix
##
phenos <- readr::read_csv("http://www.helianthome.org/rest/population/phenotype_matrix/wild",
                          show_col_types = FALSE)
# ## retrieve the list of populations (for population_id)
# ##
# populations <- readr::read_csv("http://www.helianthome.org/rest/population/list/",
#                           show_col_types = FALSE)
## retrieve individual id and population id
##
index <- readr::read_csv("http://www.helianthome.org/rest/individual/list/",
                          show_col_types = FALSE)
## merge population id into phenotypic data
##
wild <- merge(phenos, index, by = "individual_id") |>
  dplyr::select(`Leaf C N ratio`, SLA, `Primary branches`,
                `Stem diameter at flowering`, `Ligule length`,
                `Phyllaries length`, population_id) |>
  dplyr::rename(Leaf.C.N.ratio = `Leaf C N ratio`,
                Primary.branches = `Primary branches`,
                Stem.diameter.at.flowering = `Stem diameter at flowering`,
                Ligule.length = `Ligule length`,
                Phyllaries.length = `Phyllaries length`)

options(mc.cores = parallel::detectCores())

wild_CN <- brms::bf(Leaf.C.N.ratio ~ (1|population_id))
wild_SLA <- brms::bf(SLA ~ (1|population_id))
wild_branch <- brms::bf(Primary.branches ~ (1|population_id))
wild_stem <- brms::bf(Stem.diameter.at.flowering ~ (1|population_id))
wild_ligule <- brms::bf(Ligule.length ~ (1|population_id))
wild_phyllary <- brms::bf(Phyllaries.length ~ (1|population_id))

mod_wild <- brms::brm(wild_CN + wild_SLA + wild_branch + wild_stem +
                        wild_ligule + wild_phyllary +
                        brms::set_rescor(rescor = TRUE),
                      data = wild,
                      family = "gaussian",
                      refresh = 0)

m <- 4
imp <- mice::mice(wild, m = m)
mod_imp <- brms::brm_multiple(wild_CN + wild_SLA + wild_branch + wild_stem +
                                wild_ligule + wild_phyllary +
                                brms::set_rescor(rescor = TRUE),
                              data = imp,
                              family = "gaussian",
                              refresh = 0)

cat("Original dataset...\n")
tmp <- posterior::as_draws(mod_wild) |>
  posterior::subset_draws(variable = "^Intercept_", regex = TRUE) |>
  posterior::summarise_draws()
print(tmp)
tmp <- posterior::as_draws(mod_wild) |>
  posterior::subset_draws(variable = "^sd_", regex = TRUE) |>
  posterior::summarise_draws()
print(tmp)
cat("Imputed datasets (", m, ")...\n")
tmp <- posterior::as_draws(mod_imp) |>
  posterior::subset_draws(variable = "^Intercept_", regex = TRUE) |>
  posterior::summarise_draws()
print(tmp)
tmp <- posterior::as_draws(mod_imp) |>
  posterior::subset_draws(variable = "^sd_", regex = TRUE) |>
  posterior::summarise_draws(default_summary_measures())
print(tmp)
draws <- posterior::as_draws_array(mod_imp) |>
  posterior::subset_draws(variable = "^Intercept_", regex = TRUE)
draws_per_dat <- lapply(1:m, \(i) summarize_draws(draws, chain = (i-1)*m + i))
lapply(draws_per_dat, summarize_draws, default_convergence_measures())
draws <- posterior::as_draws_array(mod_imp) |>
  posterior::subset_draws(variable = "^sd_", regex = TRUE)
draws_per_dat <- lapply(1:m, \(i) summarize_draws(draws, chain = (i-1)*m + i))
lapply(draws_per_dat, summarize_draws, default_convergence_measures())


