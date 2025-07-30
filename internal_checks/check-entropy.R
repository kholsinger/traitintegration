library(traitintegration)
library(tidyverse)

rm(list = ls())

dat <- tibble(r = NA,
              k = NA,
              brute_force = NA,
              analytical  = NA)
for (r in seq(0.1, 0.9, by = 0.1)) {
  for (k in 2:20) {
    brute_force <- entropy_omega(r, k)
    analytical <- k/2 + (k/2)*log(2*pi) + (1/2)*log(1 + (k - 1)*r) +
      (1/2)*(k - 1)*log(1 - r)
    dat <- add_row(dat,
                   r = r,
                   k = k,
                   brute_force = brute_force,
                   analytical = analytical)
  }
}
dat <- dat %>%
  filter(!is.na(r)) %>%
  mutate(difference = abs(brute_force - analytical))
sink("internal_checks/check-entropy-results.txt", split = TRUE)
cat("r in seq(0.1, 0.9, by = 0.1)\n",
    "k in 2:20\n", sep = "")
cat("Maximum absolute difference: ", max(dat$difference))
sink()

