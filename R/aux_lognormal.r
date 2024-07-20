lognormal_bci <- function(alpha, lmean, lsd){
  qlnorm(p = c(1 - alpha, 1 + alpha)/2,
         meanlog = lmean,
         sdlog = lsd)
}

lognormal_hpd <- function(alpha, lmean, lsd){
  opt_int <- function(q1, a){
    q2 <- q1 + a
    cand <- qlnorm(p = c(q1, q2),
                   meanlog = lmean, sdlog = lsd)
    return(
      cand[2] - cand[1]
    )
  } 
  Opt <- optimise(
    opt_int, a = alpha,
    lower = 1E-10, upper = 1-alpha
  )
  q1.opt <- Opt$minimum
  q2.opt <- q1.opt + alpha
  return(
    qlnorm(p = c(q1.opt, q2.opt),
           meanlog = lmean, sdlog = lsd)
  )
}
lognormal_hpd_percentiles <- function(alpha, lmean, lsd){
  hpd <- lognormal_hpd(alpha = alpha, lmean, lsd)
  return(
    plnorm(q = hpd, meanlog = lmean, sdlog = lsd)
  )
}

lognormal_skewness <- function(lsd){
  v <- lsd^2
  skw <- (exp(v) + 2) * sqrt(exp(v) - 1)
  return(skw)
}
