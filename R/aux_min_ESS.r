ESS_bound_Liu_known_lognormal <- function(m,
                                v,
                                delta = 0.01,
                                p,
                                asymp_var, ## if iid, please supply p*(1-p)
                                ci.level = 0.95) {
  Sigma <- sqrt(v)
  true_Finv <- function(p)
    qlnorm(p, m, Sigma)
  true_f <- function(x)
    dlnorm(x = x, m, Sigma)
  
  
  true.csi <- true_Finv(p) ## true quantile
  
  Z <- qnorm(p = (1 + ci.level) / 2)
  R <- (Z ^ 2 * asymp_var) /
    (delta * true.csi * true_f(true.csi)) ^ 2
  S <- floor(R) + 1
  
  return(S)
}


min_samples_iid_CLT <- function(q, dens, Qf,
                            r_precision = 0.01, level = 0.95){
  Z <- qnorm(p = (1 + level)/2)
  num <- Z^2 * q * (1 - q)
  xi_q <- Qf(q)
  denom <- (r_precision * xi_q * dens(xi_q))^2
  N <- floor(num/denom) + 1
  return(N)
}

cdf_order_2 <- function(x, k, n, m, s){
  thelogF <- function(x) plnorm(x, m, s, log.p = TRUE)
  theP <- exp(thelogF(x))
  res <- pbinom(q = n, size = n, prob = theP) - pbinom(q = k -1, size = n, prob = theP)
  return(res)
}

min_samples_iid_order <- function(q, m, s,
                                  r_precision = 0.01, level = 0.95, Nmax = 1e7){
  iniM <- round(1/q)
  M <- iniM
  xi_q <- qlnorm(q, meanlog = m, sdlog = s)
  A <- (1 - r_precision) * xi_q
  B <- (1 + r_precision) * xi_q
  get_prob <- function(M){
    theK <- round(q * M)
    prob <- cdf_order_2(B, k = theK, n = M, m = m, s = s) - cdf_order_2(A, k = theK, n = M, m = m, s = s)
    return(prob)
  }
  obj <- function(M) return((get_prob(round(M)) - level)^2)
  Opt <- optimise(obj, interval = c(iniM, Nmax))
  resN <- round(Opt$minimum)
  out <- list(
    N = resN,
    error = Opt$objective,
    attained_prob = get_prob(resN)
  )
  return(out)
}
