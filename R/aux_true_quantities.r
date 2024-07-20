#' Compute the (true) mean of a log-normal distribution
#'
#' @param lmean log mean
#' @param lsd log standard deviation
#'
#' @return the mean
#' @export lognormal_mean
#'
#' @examples 
#' lognormal_mean(0, 1)
#' mean(exp(rnorm(1e6)))
lognormal_mean <- function(lmean, lsd){
 exp(lmean + lsd^2/2) 
}
#' Compute the (true) median of a log-normal distribution
#'
#' @param lmean log mean
#' @param lsd log standard deviation
#'
#' @return the median
#' @export lognormal_median
#'
#' @examples
#' lognormal_median(0, 1)
#' median(exp(rnorm(1e6)))
lognormal_median <- function(lmean, lsd){
  exp(lmean) 
}
#' Compute the (true) variance of a log-normal distribution
#'
#' @param lmean log mean
#' @param lsd log standard deviation
#'
#' @return the variance
#' @export lognormal_variance
#'
#' @examples 
#' lognormal_variance(0, 1)
#' var(exp(rnorm(1e6)))
lognormal_variance <- function(lmean, lsd){
   (exp(lsd^2) - 1) * exp(2*lmean + lsd^2) 
}
#' Compute the (true) HPD of a log-normal distribution
#'
#' @param alpha the HPD level  (default is 0.95)
#' @param lmean log mean
#' @param lsd log standard deviation
#'
#' @return a vector of size two with the HPD endpoints
#' @export lognormal_hpd
#'
#' @examples
#' lognormal_hpd(0.95, 0, 1)
#' lognormal_hpd(0.90, 0, 1)
lognormal_hpd <- function(alpha = 0.95, lmean, lsd) {
  opt_int <- function(q1, a) {
    q2 <- q1 + a
    cand <- qlnorm(p = c(q1, q2),
                   meanlog = lmean,
                   sdlog = lsd)
    return(cand[2] - cand[1])
  }
  Opt <- optimise(opt_int,
                  a = alpha,
                  lower = 1E-10,
                  upper = 1 - alpha)
  q1.opt <- Opt$minimum
  q2.opt <- q1.opt + alpha
  return(qlnorm(
    p = c(q1.opt, q2.opt),
    meanlog = lmean,
    sdlog = lsd
  ))
}
#' Compute the (true) central BCI for a log-normal distribution
#'
#' @param alpha the HPD level  (default is 0.95)
#' @param lmean log mean
#' @param lsd log standard deviation
#'
#' @return a vector of size two with the BCI endpoints
#' @export lognormal_bci
#'
#' @examples
#' lognormal_bci(0.95, 0, 1)
#' lognormal_bci(0.90, 0, 1)
lognormal_bci <- function(alpha, lmean, lsd) {
  ps <- c((1 - alpha), (1 + alpha))/2
  return(qlnorm(p = ps, meanlog = lmean, sdlog = lsd))
}

#'  Compute the (true) HPD of an exponential distribution
#'
#' @param alpha the HPD level  (default is 0.95)
#' @param theta the rate parameter.
#'
#' @return a vector of size two with the HPD endpoints
#' @export exponential_hpd
#'
#' @examples
#' exponential_hpd(0.95, 1)
#' exponential_hpd(0.90, 1)
#' exponential_hpd(0.95, 10)
#' exponential_hpd(0.90, 10)
exponential_hpd <- function(alpha = 0.95, theta) {
  opt_int <- function(q1, a) {
    q2 <- q1 + a
    cand <- qexp(p = c(q1, q2), rate = theta)
    return(cand[2] - cand[1])
  }
  Opt <- optimise(opt_int,
                  a = alpha,
                  lower = 1E-10,
                  upper = 1 - alpha)
  q1.opt <- Opt$minimum
  q2.opt <- q1.opt + alpha
  return(qexp(
    p = c(q1.opt, q2.opt), rate = theta))
}
#'  Compute the (true) central BCI of an exponential distribution
#'
#' @param alpha the BCI level  (default is 0.95)
#' @param theta the rate parameter.
#'
#' @return a vector of size two with the HPD endpoints
#' @export exponential_bci
#'
#' @examples
#' exponential_bci(0.95, 1)
#' exponential_bci(0.90, 1)
#' exponential_bci(0.95, 10)
#' exponential_bci(0.90, 10)
exponential_bci <- function(alpha = 0.95, theta) {
  ps <- c((1 - alpha), (1 + alpha))/2
  return(qexp(p = ps, rate = theta))
}

