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
