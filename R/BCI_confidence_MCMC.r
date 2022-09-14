library(tidyverse)
source("aux.r")
##############################
############# Functions
quantile_CI <- function(n, q, alpha = 0.95) {
  # Copied from https://stats.stackexchange.com/questions/99829/how-to-obtain-a-confidence-interval-for-a-percentile
  # Search over a small range of upper and lower order statistics for the 
  # closest coverage to 1-alpha (but not less than it, if possible).
  gg <- 1-alpha
  u <- qbinom(1-gg/2, n, q) + (-2:2) + 1
  l <- qbinom(gg/2, n, q) + (-2:2)
  u[u > n] <- Inf
  l[l < 0] <- -Inf
  coverage <- outer(l, u, function(a, b) pbinom(b-1,n,q) - pbinom(a-1,n,q))
  if (max(coverage) < 1-gg) i <- which(coverage == max(coverage)) else
    i <- which(coverage == min(coverage[coverage >= 1-gg]))
  i <- i[1]
  # Return the order statistics and the actual coverage.
  u <- rep(u, each=5)[i]
  l <- rep(l, 5)[i]
  return(list(Interval=c(l,u), Coverage=coverage[i]))
}
##
get_bci_intervals_Doss <- function(X,
                                   b_level = 0.95,
                                   ci_level = 0.95){
  
  Z <- qnorm(p = (1 + ci_level)/2)
  
  BL.ests <- mcmcse::mcse.q(x = X,
                            q = (1-b_level)/2)
  BU.ests <- mcmcse::mcse.q(x = X,
                            q = (1+b_level)/2)
  
  L.stuff <- BL.ests$est[1] + Z*c(-1, 1)*BL.ests$se[1] 
  U.stuff <- BU.ests$est[1] + Z*c(-1, 1)*BU.ests$se[1]  
  
  res <- tibble::tibble(
    L_lwr1 = ifelse(is.na(L.stuff[1]), -Inf, L.stuff[1]),
    L_mean1 = BL.ests$est,
    L_upr1 = L.stuff[2],
    U_lwr1 = U.stuff[1],
    U_mean1 = BU.ests$est,
    U_upr1 = ifelse(is.na(U.stuff[2]), Inf, U.stuff[2])
  )
  return(res)
}
##
get_bci_intervals <- function(X,
                              b_level = 0.95,
                              ci_level = 0.95){
  
  BCI <- stats::quantile(X, probs = c(1-b_level, 1 + b_level)/2)
  
  phat.LB <- estimate_p(q = BCI[1], X = X)
  phat.UB <- estimate_p(q = BCI[2], X = X)  
  sX <- sort(X) 
  L.stuff <- sX[quantile_CI(n = length(X), q = phat.LB,
                            alpha = ci_level)$Interval]
  U.stuff <- sX[quantile_CI(n = length(X), q = phat.UB,
                            alpha = ci_level)$Interval]
  
  res <- tibble::tibble(
    hatp_L = phat.LB,
    L_lwr2 = ifelse(is.na(L.stuff[1]), -Inf, L.stuff[1]),
    L_mean2 = BCI[1],
    L_upr2 = L.stuff[2],
    hatp_U = phat.UB,
    U_lwr2 = U.stuff[1],
    U_mean2 = BCI[2],
    U_upr2 = ifelse(is.na(U.stuff[2]), Inf, U.stuff[2])
  )
  return(res)
}
##
generate_fake_MCMC <- function(N, mu, sigma, phi, LN = FALSE){
  sigma_e <- sigma*sqrt(1-phi^2) 
  rawSimu <- arima.sim(n = N, list(ar = c(phi)), sd = sigma_e) + mu 
  if(LN){
    exampleSimu <- exp(rawSimu)
  }else{
    exampleSimu <- rawSimu
  }  
  return(exampleSimu) 
}
##
run_once <- function(i, Mu, Sigma, Phi, Alpha, Omega, M, logn){
  
  Samples <- generate_fake_MCMC(N = M,
                                mu = Mu, sigma = Sigma, phi = Phi,
                                LN = logn)
  
  Meeker.ests <- get_bci_intervals(X = Samples, b_level = Alpha,
                                    ci_level = Omega)
  Doss.ests <- get_bci_intervals_Doss(X = Samples, b_level = Alpha,
                                      ci_level = Omega)
  
  out <- tibble::tibble(
    min = min(Samples),
    max = max(Samples),
    Doss.ests,
    Meeker.ests,
    replicate = i
  )
  return(out)
}
##########################################
### Asymmetric : (2.58, 0.254)
### Symmetric : (4.01, 0.0352)
logNormal <- TRUE ## should the target be log-normal (TRUE) or normal (FALSE)?

Mu <- 2.58
Sigma <- .254
eff <- .1
Phi <- (1-eff)/(eff + 1) ## if Phi = 0, independent samples

Alpha <- 0.95
TargetCoverage <- .95
M <- 500
Nrep <- 1000

if(logNormal){
  true.BCI <- qlnorm(p = c(1-Alpha, 1 + Alpha)/2,
                     mean = Mu, sd = Sigma)
}else{
  true.BCI <- qnorm(p = c(1-Alpha, 1 + Alpha)/2,
                     mean = Mu, sd = Sigma)
}

simus <- do.call(rbind,
                 parallel::mclapply(seq_len(Nrep),
                                    function(i){
                                      run_once(i, M = M,
                                               Mu = Mu,
                                               Sigma = Sigma,
                                               Phi = Phi,
                                               Alpha = Alpha,
                                               Omega = TargetCoverage,
                                               logn = logNormal)
                                    }, mc.cores = 6)
)

simus <- simus %>%
  add_column(l_coverage_1 = is_in(
    x = true.BCI[1],
    l = simus$L_lwr1,
    u = simus$L_upr1
  ), .after = "L_upr1")

simus <- simus %>%
  add_column(u_coverage_1 = is_in(
    x = true.BCI[2],
    l = simus$U_lwr1,
    u = simus$U_upr1
  ), .after = "U_upr1")

simus <- simus %>%
  add_column(l_coverage_2 = is_in(
    x = true.BCI[1],
    l = simus$L_lwr2,
    u = simus$L_upr2
  ), .after = "L_upr2")

simus <- simus %>%
  add_column(u_coverage_2 = is_in(
    x = true.BCI[2],
    l = simus$U_lwr2,
    u = simus$U_upr2
  ), .after = "U_upr2")


widths.L.1 <- simus$L_upr1-simus$L_lwr1
widths.U.1 <- simus$U_upr1-simus$U_lwr1
widths.L.2 <- simus$L_upr2-simus$L_lwr2
widths.U.2 <- simus$U_upr2-simus$U_lwr2

par(mfrow = c(2, 2))
hist(widths.L.1, probability = TRUE,
     xlab = expression(W[L]^1), main = "Width lwr Doss")
hist(widths.L.2, probability = TRUE,
     xlab = expression(W[L]^2), main = "Width lwr Meeker")
hist(widths.U.1, probability = TRUE,
     xlab = expression(W[U]^1), main = "Width upr Doss")
hist(widths.U.2, probability = TRUE,
     xlab = expression(W[U]^2), main = "Width upr Meeker")

results.1 <- tibble::tibble(
  method = "Doss",
  quantity = c("lower_BCI", "upper_BCI"),
  coverage = c(mean(simus$l_coverage_1), mean(simus$u_coverage_1)),
  mean_interval_width = c(mean(widths.L.1), mean(widths.U.1))
)
results.2 <- tibble::tibble(
  method = "Meeker",
  quantity = c("lower_BCI", "upper_BCI"),
  coverage = c(mean(simus$l_coverage_2), mean(simus$u_coverage_2)),
  mean_interval_width = c(mean(widths.L.2), mean(widths.U.2))
)

results <- rbind(results.1, results.2)
results

eff*M
##
compute_MRAE(hats = simus$L_mean1, theta0 = true.BCI[1])
compute_MRAE(hats = simus$L_mean2, theta0 = true.BCI[1])
compute_MRAE(hats = simus$U_mean1, theta0 = true.BCI[2])
compute_MRAE(hats = simus$U_mean2, theta0 = true.BCI[2])
## 
compute_MSE(hats = simus$L_mean1, theta0 = true.BCI[1])
compute_MSE(hats = simus$L_mean2, theta0 = true.BCI[1])
compute_MSE(hats = simus$U_mean1, theta0 = true.BCI[2])
compute_MSE(hats = simus$U_mean2, theta0 = true.BCI[2])