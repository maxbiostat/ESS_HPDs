#' Estimate empirical CDF 
#'
#' @param q A value on the domain for which one wants p = F(X <= q)
#' @param X A sample from F
#'
#' @return an estimate of p 
#' @export estimate_p
#'
#' @examples
#' estimate_p(0, rnorm(100))
estimate_p <- function(q, X) {
  mean(X <= q)
}

#' Compute the mean squared error (MSE)
#'
#' @param Xs values of the estimator
#' @param true.qoi true value
#'
#' @return mean absolute relative error (MRAE)
#' @export mse
#'
mse <- function(Xs, true.qoi){
  Ys <- (Xs - true.qoi)^2
  return(mean(Ys, na.rm = TRUE))
}

#' Compute the mean absolute relative error (MRAE)
#'
#' @param Xs values of the estimator
#' @param true.qoi true value
#'
#' @return mean absolute relative error (MRAE)
#' @export mrae
#'
mrae <- function(Xs, true.qoi){
  Ys <- abs(Xs - true.qoi)/abs(true.qoi)
  return(mean(Ys, na.rm = TRUE))
}

#' Check if a value is inside a given interval
#'
#' @param x value to be tested
#' @param l lower endpoint of the interval
#' @param u  upper endpoint of the interval
#'
#' @return whether x is in (l, u)
#' @export is_in
#'
is_in <- function(x, l, u){
  below <- x >= l
  above <- x <= u
  result <- as.logical(below * above)
  return(result)
}

#' Compute BCI and HPD on a set of observations/draws
#'
#' @param X a set of draws from some process
#' @param alpha the level of confidence (between 0 and 1)
#'
#' @return data.frame containing BCI, HPD and approximate normal BCI/HPD
#' @export get_intervals
#'
get_intervals <- function(X, alpha = 0.95){
  
  BCI <- stats::quantile(X, probs = c(1-alpha, 1 + alpha)/2)
  HPD <- HDInterval::hdi(X, credMass = alpha)
  
  mm <- mean(X)
  vv <- var(X)
  
  Nint <- qnorm(p = c(1 - alpha, 1 + alpha)/2,
                mean = mm, sd = sqrt(vv)) ## normal interval
  
  out <- data.frame(
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