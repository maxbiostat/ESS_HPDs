is_in <- function(x, l, u){
  below <- x >= l
  above <- x <= u
  result <- as.logical(below * above)
  return(result)
}
get_intervals <- function(X, alpha = 0.95){
  BCI <- stats::quantile(X, probs = c(1-alpha, 1 + alpha)/2)
  HPD <- HDInterval::hdi(X, credMass = alpha)
  
  mm <- mean(X)
  vv <- var(X)
  
  Nint <- qnorm(p = c(1-alpha, 1 + alpha)/2, mean = mm, sd = sqrt(vv))
  
  out <- tibble::tibble(
    mean = mm,
    median = median(X),
    variance = vv,
    n = length(X),
    bci_lwr = as.numeric(BCI[1]),
    bci_upr = as.numeric(BCI[2]),
    bci_width = as.numeric(BCI[2]) - as.numeric(BCI[1]),
    hpd_lwr = as.numeric(HPD[1]),
    hpd_upr = as.numeric(HPD[2]),
    hpd_width = as.numeric(HPD[2]) - as.numeric(HPD[1]),
    napp_lwr =  Nint[1],
    napp_upr =  Nint[2],
    napp_width =  Nint[2] - Nint[1]
  )
  return(out)
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
    c(q1.opt, q2.opt)
  )
}
compute_MRAE <- function(hats, theta0){
  mean(abs(1-(hats/theta0)))
}
compute_MSE <- function(hats, theta0){
  VV <- var(hats)
  BB <- mean(hats) - theta0
  MSE <- mean((hats-theta0)^2)
  return(c(
    est_variance = VV,
    est_bias = BB,
    mse = MSE
  ))
}
estimate_p <- function(q, X){
  mean(X <= q)
}