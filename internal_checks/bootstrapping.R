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
## merge population id into phenotypic data. drop individuals lacking
## complete data
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
  dplyr::filter(population_id != "ANN_62" &
                population_id != "ANN_64" &
                population_id != "ANN_71")

scale_within_population <- function(dat, pop) {
  tmp <- subset(dat, population_id == pop)
  tmp <- as.data.frame(scale(tmp[, 1:6]))
  col_names <- c(colnames(tmp), "population_id")
  tmp <- tibble::add_column(tmp, pop)
  colnames(tmp) <- col_names
  return(tmp)
}

scale_within_populations <- function(dat) {
  pops <- unique(wild$population_id)
  tmp <- scale_within_population(wild, pops[1])
  for (i in 2:length(pops)) {
    tmp <- tibble::add_row(tmp,
                           scale_within_population(wild, pops[i]))
  }
  return(tmp)
}

correlation_within_populations <- function(dat) {
  tmp <- scale_within_populations(dat)
  R <- cor(tmp[, 1:6])
  return(R)
}

theta <- function(x, x_data) {
 R <- correlation_within_populations(dat[x, 1:6])
}

boot_results <- bootstrap::bootstrap(1:nrow(wild), 4000, theta, x_data)

## boot_results is a matrix with 36 columns and 4000 rows. That will
## need to be converted to a format compatible with the array that
## the other functions are expecting
##
## but boot_results isn't right. It's returning the same result every time
## I'll hard code it



