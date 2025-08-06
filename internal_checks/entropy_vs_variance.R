library(traitintegration)
library(ggplot2)
library(tibble)

rm(list = ls())

my_var <- function(x) {
  x_bar <- mean(x)
  n <- length(x)
  sum_sq <- sum((x - x_bar)^2)
  return(sum_sq/n)
}

r_vals <- seq(0, 0.999, by = 0.001)
k <- 5
R <- matrix(nrow = k, ncol = k)
ent_R <- numeric(length(r_vals))
var_l <- numeric(length(r_vals))
for (i in 1:length(r_vals)) {
  for (j in 1:(k-1)) {
    R[j, j] <- 1.0
    for (k in (j+1):k) {
      R[j, k] <- r_vals[i]
      R[k, j] <- R[j, k]
    }
  }
  R[k, k] <- 1.0
  R_array <- array(dim = c(1, k, k))
  R_array[1, , ] <- R
  ent_R[i] <- get_ent_R(R_array)
  var_l[i] <- my_var(eigen(R)$values)/(k - 1)
}

for_plot <- tibble(r = c(r_vals, r_vals),
                   value = c(ent_R, var_l),
                   Measure = c(rep("entR", length(r_vals)),
                               rep("var(lambda)", length(r_vals))))
p <- ggplot(for_plot, aes(x = r, y = value, color = Measure)) +
  geom_line() +
  theme_bw()
print(p)


