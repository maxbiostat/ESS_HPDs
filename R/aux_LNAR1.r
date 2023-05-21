autocorr_expAR1 <- function(k, m, v, rho){
  uy.sq <- exp(2*m + v)
  vw <- 2*v*(1 + rho^k) 
  EW <- exp(2*m + vw/2)
  covar <- EW - uy.sq
  vy <- (exp(v)-1)*exp(2*m + v)
  ans <- covar/vy
  return(ans)
}
autocorr_expAR1 <- Vectorize(autocorr_expAR1)

LNAR1_eff <- function(mean, variance, rho){
  covk <- function(k) autocorr_expAR1(k = k,
                                      m = mean, v = variance, rho = rho)
  theoEffy <-  1/ (1 + 2 * sum(covk(1:1E5)))
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
