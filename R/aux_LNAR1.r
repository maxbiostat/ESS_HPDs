### Latent autoregressive (LAR) stuff; the LAR is a **binary** process
LAR_autocorr_2 <- function(l, a, phi, sigma_e){
  ## computes the autocorrelation function of a LAR of the type
  ### exp(Z_t) <= log(a), where Z_t is an AR(1, m, v, rho).
  sigma_m <- sigma_e/sqrt(1 - phi^2)
  covk <-  phi^l * sigma_m^2
  Sigma <- matrix(c(sigma_m^2, covk, covk, sigma_m^2), nrow = 2)
  P1 <- mvtnorm::pmvnorm(lower = c(-Inf, -Inf),
                         upper = c(a, a),
                         mean = rep(0, 2),
                         sigma = Sigma)
  P2 <- mvtnorm::pmvnorm(lower = c(-Inf, a),
                         upper = c(a, Inf),
                         mean = rep(0, 2),
                         sigma = Sigma)
  ans <- (P1[1]/pnorm(q = a, sd = sigma_m)) -
    (P2[1]/pnorm(q = a, sd = sigma_m, lower.tail = FALSE))
  return(ans)
}
LAR_autocorr_2 <- Vectorize(LAR_autocorr_2)
#################

autocov_expAR1 <- function(k, m, v, rho){
  uy.sq <- exp(2*m + v)
  EW <- exp(2 * m + v * (1 + rho^k))
  ans <- EW - uy.sq
  return(ans)
}
autocov_expAR1 <- Vectorize(autocov_expAR1)

autocorr_expAR1 <- function(k, m, v, rho){
  covar <- autocov_expAR1(k, m, v, rho)
  vy <- (exp(v) - 1) * exp(2*m + v)
  ans <- covar/vy
  return(ans)
}
autocorr_expAR1 <- Vectorize(autocorr_expAR1)

LNAR1_eff <- function(mean, variance, rho, cap = 1e5){
  corrk <- function(k) autocorr_expAR1(k = k,
                                       m = mean, v = variance, rho = rho)
  theoEffy <-  1 / (1 + 2 * sum(corrk(1:cap)))
  return(theoEffy)
}

find_rho <- function(eff, mean, variance){
  opt_eff <- function(cand){
    theo <- LNAR1_eff(mean = mean, variance = variance, rho = cand)
    return(
      (eff - theo)^2
    )
  }
  opt_eff <- Vectorize(opt_eff)
  Opt <- optimise(opt_eff, interval = c(-1, 1))
  return(Opt$minimum)
}

LNAR_quantile_autocorr <- function(k, xi_p, m, v, rho, K){
  
  ## Depends on the LAR
  sgE <- sqrt(v) * sqrt(1 - rho ^ 2)
  
  g_j <- function(j){
    LAR_autocorr_2(l = j,
                   a = log(xi_p) - m,
                   sigma_e = sgE,
                   phi = rho)
  }
  g_j <- Vectorize(g_j)
  
  f_j <- function(j){
    return((1 - j/K) * g_j(j))
  }
  f_j <- Vectorize(f_j)
  
  return(f_j(k))
}

compute_LT_variance <- function(p, x, m, v, rho, K = 5000){
  
  autcorrs <- LNAR_quantile_autocorr(k = 1:(K - 1),
                                     xi_p = x,
                                     m = m,
                                     v = v,
                                     rho = rho,
                                     K = K)
  
  SigmaSq_p <- p * (1 - p) * (1 + 2 * sum(autcorrs))
  
  return(list(asymp_var = SigmaSq_p,
              iid_var = p * (1 - p)))
}
