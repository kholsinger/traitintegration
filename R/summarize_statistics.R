#' Summarize results returned by `assemble_statistics()`
#'
#' @description
#' This function provides a quick summary of results returned by
#' `assemble_statistics`. The summary includes the posterior mean and
#' quantiles for the leading eigenvalue, the variance of eigenvalues, and
#' entR. By default only the 2.5% and 97.5% quantiles are reported, but
#' any numeric vector of probabilities in [0, 1] can be specified. See the
#' documentation for `quantile()` for additional details.
#'
#' @param results A data fram returned by assemble statistics
#' @param probs Numeric vector of probabilities for quantiles. Default:
#' c(0.025, 0.975)
#' @param digits Number of digits to report in results. Default: 3
#'
#' @examples
#' tmp <- assemble_statistics(mod_wild)
#' summarize_statistics(tmp)
#'
#' @export
summarize_statistics <- function(results, probs = c(0.025, 0.975),
                                 digits = 3)
{
  na_vec <- rep(NA, 3)
  dat <- data.frame(Statistic = c("Lead_EV", "Var_EV", "r_eff"),
                    Mean = na_vec)

  ## There is probably a better way to add named columns for the quantiles,
  ## but I haven't found it
  ##
  col_names <- paste(probs*100, "%", sep = "")
  tmp <- matrix(ncol = length(probs), data = rep(NA, length(probs)*3))
  colnames(tmp) <- col_names
  dat <- cbind(dat, tmp)

  Lead_EV_stats <- get_summary(results$Lead_EV, probs)
  Var_EV_stats <- get_summary(results$Var_EV, probs)
  r_eff_stats <- get_summary(results$r_eff, probs)

  dat$Mean[1] <- round(Lead_EV_stats$Mean, digits)
  dat$Mean[2] <- round(Var_EV_stats$Mean, digits)
  dat$Mean[3] <- round(r_eff_stats$Mean, digits)

  for (i in 1:length(probs)) {
    dat[1, i+2] <- round(Lead_EV_stats$quantile[i], digits)
    dat[2, i+2] <- round(Var_EV_stats$quantile[i], digits)
    dat[3, i+2] <- round(r_eff_stats$quantile[i], digits)
  }
  return(dat)
}

#' Return posterior mean and specified quantiles for one of the trait
#' integration measures
#'
#' @param result_vector The vector of sampling results for a trait integration
#' measure
#' @param probs Numeric vector of probabilities for the quantiles
#'
#' @return A named list with components Mean (corresponding to the posterior
#' mean) and quantile (a vector with components corresponding to the
#' specified probabilities
#'
get_summary <- function(result_vector, probs) {
  x_bar <- mean(result_vector)
  q_tile <- quantile(result_vector, probs)
  return(list(Mean = x_bar, quantile = q_tile))
}
