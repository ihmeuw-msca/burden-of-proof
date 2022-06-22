egger_regression <- function(residual, residual_sd, one_sided = TRUE) {
  weighted_residual <- residual/residual_sd
  r_mean <- mean(weighted_residual)
  r_sd <- 1/sqrt(length(weighted_residual))
  r_pval <- get_pval(r_mean, r_sd, one_sided = one_sided)
  list(mean = r_mean, sd = r_sd, pval = r_pval)
}

get_pval <- function(beta, beta_sd, one_sided = FALSE) {
  zscore <- abs(beta/beta_sd)
  if (one_sided) {
    pval <- 1 - pnorm(zscore)
  } else {
    pval <- 2*(1 - pnorm(zscore))
  }
  pval
}