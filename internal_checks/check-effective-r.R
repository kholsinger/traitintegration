library(traitintegration)

rm(list = ls())

get_var <- function(Omega) {
  n_iter <- dim(Omega)[1]
  e_var <- numeric(n_iter)
  for (i in 1:n_iter) {
    values <- eigen(Omega[i, ,])$values
    e_var[i] <- var(values)
  }
  return(e_var)
}

var_omega <- function(r, k) {
  var_omega <- diag(k)
  for (i in 1:(k-1)) {
    for (j in (i+1):k) {
      var_omega[i, j] <- r
      var_omega[j, i] <- r
    }
  }
  values <- eigen(var_omega)$values
  return(var(values))
}

get_r_eff_var <- function(Omega, tol = 1.0e-14) {
  n_var <- dim(Omega)[1]
  k <- dim(Omega)[2]
  var_obs <- get_var(Omega)
  r_eff <- numeric(n_var)
  diff <- 1.0
  min_r <- 1.0e-16
  max_r <- 1 - min_r
  var_min_r <- 0
  var_max_r <- k - 1
  for (i in 1:n_var) {
    min <- min_r
    max <- max_r
    while (diff > tol)  {
      mid <- (max + min)/2.0
      f_mid_r <- var_omega(mid, k)
      if (f_mid_r > var_obs[i]) {
        max <- mid
        f_max_r <- f_mid_r
      } else {
        min <- mid
        f_min_r <- f_mid_r
      }
      diff <- abs(f_mid_r - var_obs[i])
    }
    r_eff[i] <- mid
    diff <- 1.0
  }
  return(r_eff)
}

get_lead <- function(Omega) {
  n_iter <- dim(Omega)[1]
  l_var <- numeric(n_iter)
  for (i in 1:n_iter) {
    values <- eigen(Omega[i, ,])$values
    l_var[i] <- max(values)
  }
  return(l_var)
}

lead_omega <- function(r, k) {
  lead_omega <- diag(k)
  for (i in 1:(k-1)) {
    for (j in (i+1):k) {
      lead_omega[i, j] <- r
      lead_omega[j, i] <- r
    }
  }
  values <- eigen(lead_omega)$values
  return(max(values))
}

get_r_eff_lead <- function(Omega, tol = 1.0e-14) {
  n_lead <- dim(Omega)[1]
  k <- dim(Omega)[2]
  lead_obs <- get_lead(Omega)
  r_eff <- numeric(n_lead)
  diff <- 1.0
  min_r <- 1.0e-16
  max_r <- 1 - min_r
  var_min_r <- 0
  var_max_r <- k - 1
  for (i in 1:n_lead) {
    min <- min_r
    max <- max_r
    while (diff > tol)  {
      mid <- (max + min)/2.0
      f_mid_r <- var_omega(mid, k)
      if (f_mid_r > lead_obs[i]) {
        max <- mid
        f_max_r <- f_mid_r
      } else {
        min <- mid
        f_min_r <- f_mid_r
      }
      diff <- abs(f_mid_r - lead_obs[i])
    }
    r_eff[i] <- mid
    diff <- 1.0
  }
  return(r_eff)
}

Omega <- get_Omega(mod_wild)

cat("effective r, entropy calculation...\n")
entropy <- get_r_eff(Omega)
cat("effective r, variance calculation...\n")
variance <- get_r_eff_var(Omega)
cat("effective r, lead eigenvalue calculation...\n")
lead <- get_r_eff_lead(Omega)

df <- data.frame(entropy = entropy,
                 variance = variance,
                 lead = lead)
summary(df)
